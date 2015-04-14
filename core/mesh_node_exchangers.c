// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/mesh.h"
#include "core/unordered_set.h"
#include "core/hilbert.h"
#include "core/kd_tree.h"
#include "core/array_utils.h"

exchanger_t* mesh_1v_node_exchanger_new(mesh_t* mesh)
{
  exchanger_t* ex = exchanger_new(mesh->comm);
#if POLYMEC_HAVE_MPI

  int nprocs;
  MPI_Comm_size(mesh->comm, &nprocs);
  if (nprocs == 1)
    return ex;

  // First construct the n-valued face exchanger and the associated offsets.
  int offsets[mesh->num_nodes+1];
  exchanger_t* noden_ex = mesh_nv_node_exchanger_new(mesh, offsets);
exchanger_fprintf(noden_ex, stdout);

  // Determine the owner of each node. We assign a face to the process on 
  // which it is present, having the lowest rank.
  int rank;
  MPI_Comm_rank(mesh->comm, &rank);
printf("%d: node offsets = {", rank);
for (int n = 0; n < mesh->num_nodes; ++n)
  printf("%d ", offsets[n]);
printf("}\n");
  int node_procs[offsets[mesh->num_nodes]];
  for (int n = 0; n < mesh->num_nodes; ++n)
  {
    for (int i = offsets[n]; i < offsets[n+1]; ++i)
      node_procs[i] = rank;
  }
  exchanger_exchange(noden_ex, node_procs, 1, 0, MPI_INT);
printf("%d: node owners = {", rank);
for (int n = 0; n < offsets[mesh->num_nodes]; ++n)
  printf("%d ", node_procs[n]);
printf("}\n");
#ifndef NDEBUG
  for (int n = 0; n < offsets[mesh->num_nodes]; ++n)
  {
    ASSERT(node_procs[n] < nprocs);
  }
#endif

  // Now set up our single-valued exchanger.
  int_ptr_unordered_map_t* send_map = int_ptr_unordered_map_new();
  int_ptr_unordered_map_t* receive_map = int_ptr_unordered_map_new();
  for (int n = 0; n < mesh->num_nodes; ++n)
  {
    int node_sender = INT_MAX;
    for (int i = offsets[n]; i < offsets[n+1]; ++i)
      node_sender = MIN(node_sender, node_procs[i]);
    if (node_sender == rank)
    {
      for (int i = offsets[n]; i < offsets[n+1]; ++i)
      {
        int node_receiver = node_procs[i];
        if (node_receiver != node_sender)
        {
          int_array_t** send_p = (int_array_t**)int_ptr_unordered_map_get(send_map, node_receiver);
          int_array_t* send = NULL;
          if (send_p == NULL)
          {
            send = int_array_new();
            int_ptr_unordered_map_insert_with_v_dtor(send_map, node_receiver, send, DTOR(int_array_free));
          }
          else
            send = *send_p;
          int_array_append(send, n);
        }
      }
    }
    else
    {
      ASSERT(node_sender < rank);
      int_array_t** receive_p = (int_array_t**)int_ptr_unordered_map_get(receive_map, node_sender);
      int_array_t* receive = NULL;
      if (receive_p == NULL)
      {
        receive = int_array_new();
        int_ptr_unordered_map_insert_with_v_dtor(receive_map, node_sender, receive, DTOR(int_array_free));
      }
      else
        receive = *receive_p;
      int_array_append(receive, n);
    }
  }
  exchanger_set_sends(ex, send_map);
  exchanger_set_receives(ex, receive_map);
  int_ptr_unordered_map_free(send_map);
  int_ptr_unordered_map_free(receive_map);
  exchanger_free(noden_ex);
#endif
  return ex;
}

exchanger_t* mesh_nv_node_exchanger_new(mesh_t* mesh, int* node_offsets)
{
  exchanger_t* ex = exchanger_new(mesh->comm);

  // Initialize the node offsets array.
  node_offsets[0] = 0;
  for (int n = 0; n <= mesh->num_nodes; ++n)
    node_offsets[n+1] = node_offsets[n] + 1;

#if POLYMEC_HAVE_MPI

  int nprocs;
  MPI_Comm_size(mesh->comm, &nprocs);
  if (nprocs == 1)
    return ex;

  int rank;
  MPI_Comm_rank(mesh->comm, &rank);

  // Create a 2-valued face exchanger.
  exchanger_t* face_ex = mesh_2v_face_exchanger_new(mesh);

  // Generate a list of all the neighbors of our set of neighbors.
  int_array_t* all_neighbors_of_neighbors = int_array_new();
  int num_neighbors = exchanger_num_sends(face_ex);
  int neighbors[num_neighbors];
  {
    MPI_Request requests[2*num_neighbors];
    MPI_Status statuses[2*num_neighbors];

    // Identify our own neighbors.
    int pos = 0, proc, *indices, size, p = 0;
    while (exchanger_next_send(face_ex, &pos, &proc, &indices, &size))
      neighbors[p++] = proc;

    // Get the number of neighbors for our pth neighbor.
    int num_neighbors_of_neighbors[num_neighbors];
    for (int p = 0; p < num_neighbors; ++p)
    {
      MPI_Irecv(&num_neighbors_of_neighbors[p], 1, MPI_INT, neighbors[p], 0, 
                mesh->comm, &requests[p]);
      MPI_Isend(&num_neighbors, 1, MPI_INT, neighbors[p], 0, 
          mesh->comm, &requests[p + num_neighbors]);
    }
    MPI_Waitall(2 * num_neighbors, requests, statuses);

    // Get the ranks of the neighbors for our pth neighbor.
    int* neighbors_of_neighbors[num_neighbors];
    for (int p = 0; p < num_neighbors; ++p)
    {
      neighbors_of_neighbors[p] = polymec_malloc(sizeof(int) * num_neighbors_of_neighbors[p]);
      MPI_Irecv(neighbors_of_neighbors[p], num_neighbors_of_neighbors[p], 
          MPI_INT, neighbors[p], 0, mesh->comm, &requests[p]);
      MPI_Isend(neighbors, num_neighbors, MPI_INT, neighbors[p], 0, 
          mesh->comm, &requests[p + num_neighbors]);
    }
    MPI_Waitall(2 * num_neighbors, requests, statuses);

    // Mash them all together.
    int_unordered_set_t* neighbor_neighbor_set = int_unordered_set_new();
    for (int p = 0; p < num_neighbors; ++p)
    {
      int_unordered_set_insert(neighbor_neighbor_set, neighbors[p]);
      for (int pp = 0; pp < num_neighbors_of_neighbors[p]; ++pp)
      {
        if (neighbors_of_neighbors[p][pp] != rank)
        {
          int_unordered_set_insert(neighbor_neighbor_set, 
                                   neighbors_of_neighbors[p][pp]);
        }
      }
      polymec_free(neighbors_of_neighbors[p]);
    }
    int npos = 0, neighbor;
    while (int_unordered_set_next(neighbor_neighbor_set, &npos, &neighbor))
      int_array_append(all_neighbors_of_neighbors, neighbor);
    int_unordered_set_free(neighbor_neighbor_set);
  }

  // Make a list of all the nodes that can be communicated to neighboring 
  // processes.
  int_array_t* my_node_indices = int_array_new();
  point_array_t* my_nodes = point_array_new();
  {
    int_unordered_set_t* node_set = int_unordered_set_new();
    int pos = 0, proc, *indices, size;
    while (exchanger_next_send(face_ex, &pos, &proc, &indices, &size))
    {
      for (int f = 0; f < size; ++f)
      {
        int face = indices[f]/2; // see docs on face exchanger!
        int npos = 0, node;
        while (mesh_face_next_node(mesh, face, &npos, &node))
          int_unordered_set_insert(node_set, node);
      }
    }
    int npos = 0, node;
    while (int_unordered_set_next(node_set, &npos, &node))
    {
      int_array_append(my_node_indices, node);
      point_array_append(my_nodes, mesh->nodes[node]);
    }
    int_unordered_set_free(node_set);
  }
//printf("%d: nodes = {", rank);
//for (int n = 0; n < my_nodes->size; ++n)
//  printf("(%g, %g, %g) ", my_nodes->data[n].x, my_nodes->data[n].y, my_nodes->data[n].z);
//printf("}\n");

  // Sort our nodes so that they are in Hilbert order.
  {
    bbox_t bbox = {.x1 = FLT_MAX, .x2 = -FLT_MAX,
                   .y1 = FLT_MAX, .y2 = -FLT_MAX,
                   .z1 = FLT_MAX, .z2 = -FLT_MAX};
    for (int i = 0; i < my_nodes->size; ++i)
      bbox_grow(&bbox, &(my_nodes->data[i]));
    hilbert_t* curve = hilbert_new(&bbox);
    hilbert_sort_points(curve, my_nodes->data, my_node_indices->data, 
                        my_nodes->size);
  }

  // Now send/receive the positions of all nodes that can interact with 
  // neighbors of our neighbors. 
  int num_neighbor_neighbors = all_neighbors_of_neighbors->size;
  MPI_Request requests[2*num_neighbor_neighbors];
  MPI_Status statuses[2*num_neighbor_neighbors];
  int_ptr_unordered_map_t* neighbor_neighbor_nodes = int_ptr_unordered_map_new();
  {
    // How many nodes does each neighbor neighbor have?
    int num_neighbor_neighbor_nodes[num_neighbor_neighbors];
    for (int p = 0; p < num_neighbor_neighbors; ++p)
    {
      int proc = all_neighbors_of_neighbors->data[p];
      MPI_Irecv(&num_neighbor_neighbor_nodes[p], 1, MPI_INT, proc, 0, 
                mesh->comm, &requests[p]);
    }
    for (int p = 0; p < num_neighbor_neighbors; ++p)
    {
      int proc = all_neighbors_of_neighbors->data[p];
      MPI_Isend(&(my_nodes->size), 1, MPI_INT, proc, 0, 
                mesh->comm, &requests[p + num_neighbor_neighbors]);
    }
    MPI_Waitall(2 * num_neighbor_neighbors, requests, statuses);

    // Get 'em!
    for (int p = 0; p < num_neighbor_neighbors; ++p)
    {
      int proc = all_neighbors_of_neighbors->data[p];
      point_array_t* nn_nodes = point_array_new();
      point_array_resize(nn_nodes, num_neighbor_neighbor_nodes[p]);
      int_ptr_unordered_map_insert_with_v_dtor(neighbor_neighbor_nodes, proc, nn_nodes, DTOR(point_array_free));
      MPI_Irecv(nn_nodes->data, 3*num_neighbor_neighbor_nodes[p], MPI_REAL_T, proc, 0, 
                mesh->comm, &requests[p]);
    }

    // Send 'em!
    for (int p = 0; p < num_neighbor_neighbors; ++p)
    {
      int proc = all_neighbors_of_neighbors->data[p];
      MPI_Isend(my_nodes->data, 3*my_nodes->size, MPI_REAL_T, proc, 0, 
                mesh->comm, &requests[p + num_neighbor_neighbors]);
    }
    MPI_Waitall(2 * num_neighbor_neighbors, requests, statuses);
  }

  // At this point, neighbor_neighbor_nodes maps the ranks of all processes
  // we can possibly interact with to the positions of the nodes on those 
  // processes. 

  // Now we cut out the nodes that don't match up.
  {
    // Post the receives for the numbers of culled nodes.
    int num_culled_nodes[num_neighbor_neighbors];
    for (int p = 0; p < num_neighbor_neighbors; ++p)
    {
      int proc = all_neighbors_of_neighbors->data[p];
      MPI_Irecv(&(num_culled_nodes[p]), 1, MPI_INT, proc, 0, 
                mesh->comm, &requests[p]);
    }

    // Figure out which nodes we will cull and send the number.
    // culled_nodes[p] contains a list of the nodes on neighbor p 
    // that DON'T correspond to any of my_nodes.
    int_array_t* culled_nodes[num_neighbor_neighbors];
    real_t tolerance = 1e-8; // FIXME: This is bad.
    for (int p = 0; p < num_neighbor_neighbors; ++p)
    {
      int proc = all_neighbors_of_neighbors->data[p];
      point_array_t* nn_nodes = *int_ptr_unordered_map_get(neighbor_neighbor_nodes, proc);
      culled_nodes[p] = int_array_new();
      int_unordered_set_t* kept_nodes = int_unordered_set_new();
      for (int i = 0; i < my_nodes->size; ++i)
      {
        point_t* xi = &my_nodes->data[i];
        for (int j = 0; j < nn_nodes->size; ++j)
        {
          point_t* xj = &nn_nodes->data[j];
          if (point_distance(xi, xj) <= tolerance)
            int_unordered_set_insert(kept_nodes, j);
        }
      }
      for (int i = 0; i < nn_nodes->size; ++i)
      {
        if (!int_unordered_set_contains(kept_nodes, i))
          int_array_append(culled_nodes[p], i);
      }
      int_unordered_set_free(kept_nodes);
      MPI_Isend(&culled_nodes[p]->size, 1, MPI_INT, proc, 0, 
                mesh->comm, &requests[p + num_neighbor_neighbors]);
    }
    MPI_Waitall(2 * num_neighbor_neighbors, requests, statuses);

    // Now receive/send the culled nodes.
    // my_culled_nodes[p] is the list of indices of nodes in my_nodes
    // that DON'T correspond to any of the nodes on neighbor p.
    int_array_t* my_culled_nodes[num_neighbor_neighbors];
    for (int p = 0; p < num_neighbor_neighbors; ++p)
    {
      my_culled_nodes[p] = int_array_new();
      int_array_resize(my_culled_nodes[p], num_culled_nodes[p]);
      int proc = all_neighbors_of_neighbors->data[p];
      MPI_Irecv(my_culled_nodes[p]->data, num_culled_nodes[p], MPI_INT, proc, 0, 
                mesh->comm, &requests[p]);
    }
    for (int p = 0; p < num_neighbor_neighbors; ++p)
    {
      int proc = all_neighbors_of_neighbors->data[p];
      MPI_Isend(culled_nodes[p]->data, culled_nodes[p]->size, MPI_INT, proc, 0, 
                mesh->comm, &requests[p + num_neighbor_neighbors]);
    }
    MPI_Waitall(2 * num_neighbor_neighbors, requests, statuses);

//for (int p = 0; p < num_neighbor_neighbors; ++p)
//{
//printf("%d: culled nodes for %d = {", rank, all_neighbors_of_neighbors->data[p]);
//for (int n = 0; n < my_culled_nodes[p]->size; ++n)
//  printf("%d ", my_culled_nodes[p]->data[n]);
//printf("}\n");
//}
    // Organized the culled nodes into sets for querying.
    int_unordered_set_t* my_culled_node_sets[num_neighbor_neighbors];
    int_unordered_set_t* their_culled_node_sets[num_neighbor_neighbors];
    for (int p = 0; p < num_neighbor_neighbors; ++p)
    {
      my_culled_node_sets[p] = int_unordered_set_new();
      for (int i = 0; i < my_culled_nodes[p]->size; ++i)
        int_unordered_set_insert(my_culled_node_sets[p], my_culled_nodes[p]->data[i]);
      int_array_free(my_culled_nodes[p]);
      their_culled_node_sets[p] = int_unordered_set_new();
      for (int i = 0; i < culled_nodes[p]->size; ++i)
        int_unordered_set_insert(their_culled_node_sets[p], culled_nodes[p]->data[i]);
      int_array_free(culled_nodes[p]);
    }

    // First, figure out the node offsets so that we can consistently 
    // create receive nodes.
    kd_tree_t* node_tree = kd_tree_new(mesh->nodes, mesh->num_nodes);
    for (int p = 0; p < num_neighbor_neighbors; ++p)
    {
      int proc = all_neighbors_of_neighbors->data[p];
      point_array_t* their_nodes = *int_ptr_unordered_map_get(neighbor_neighbor_nodes, proc);
      int sorted_indices[their_nodes->size];
      for (int i = 0; i < their_nodes->size; ++i)
        sorted_indices[i] = i;
      int_qsort(sorted_indices, their_nodes->size);

      for (int i = 0; i < their_nodes->size; ++i)
      {
        int j = sorted_indices[i];
        if (!int_unordered_set_contains(their_culled_node_sets[p], j))
        {
          // Find the node in "their_nodes" that matches our local node.
          int node = kd_tree_nearest(node_tree, &their_nodes->data[j]);

          // Bump back all the nodes behind it.
          for (int n = node+1; n <= mesh->num_nodes; ++n)
            ++node_offsets[n];
        }
      }
    }

    // Now set up the exchanger send/receive maps and the node offsets.
    int_ptr_unordered_map_t* send_map = int_ptr_unordered_map_new();
    int_ptr_unordered_map_t* receive_map = int_ptr_unordered_map_new();
    int_int_unordered_map_t* receive_offsets = int_int_unordered_map_new();
    for (int p = 0; p < num_neighbor_neighbors; ++p)
    {
      int proc = all_neighbors_of_neighbors->data[p];

      // Send nodes:
      for (int i = 0; i < my_nodes->size; ++i)
      {
        int_array_t* send_nodes = NULL;
        if (!int_unordered_set_contains(my_culled_node_sets[p], i))
        {
          if (send_nodes == NULL)
          {
            int_array_t** send_nodes_p = (int_array_t**)int_ptr_unordered_map_get(send_map, proc);
            if (send_nodes_p != NULL)
              send_nodes = *send_nodes_p;
            else
            {
              send_nodes = int_array_new();
              int_ptr_unordered_map_insert_with_v_dtor(send_map, proc, send_nodes, DTOR(int_array_free));
            }
          }
          int node = my_node_indices->data[i];
          int_array_append(send_nodes, node_offsets[node]);
        }
      }
      int_unordered_set_free(my_culled_node_sets[p]);

      // Received nodes: received values from nodes go into the slots just 
      // above the local nodal value in an exchanged array. Of course, we 
      // have to sort these nodes by their indices so that we don't mess up 
      // the ordering of offsets.
      point_array_t* their_nodes = *int_ptr_unordered_map_get(neighbor_neighbor_nodes, proc);
      int sorted_indices[their_nodes->size];
      for (int i = 0; i < their_nodes->size; ++i)
        sorted_indices[i] = i;
      int_qsort(sorted_indices, their_nodes->size);

      for (int i = 0; i < their_nodes->size; ++i)
      {
        int j = sorted_indices[i];
        int_array_t* receive_nodes = NULL;
        if (!int_unordered_set_contains(their_culled_node_sets[p], j))
        {
          if (receive_nodes == NULL)
          {
            int_array_t** receive_nodes_p = (int_array_t**)int_ptr_unordered_map_get(receive_map, proc);
            if (receive_nodes_p != NULL)
              receive_nodes = *receive_nodes_p;
            else
            {
              receive_nodes = int_array_new();
              int_ptr_unordered_map_insert_with_v_dtor(receive_map, proc, receive_nodes, DTOR(int_array_free));
            }
          }

          // Find the node in "their_nodes" that matches our local node and 
          // append it to our list of receive nodes.
          int node = kd_tree_nearest(node_tree, &their_nodes->data[j]);
          int* offset_p = int_int_unordered_map_get(receive_offsets, node);
          int receive_index;
          if (offset_p == NULL)
            receive_index = node_offsets[node] + 1;
          else
            receive_index = *offset_p + 1;
          ASSERT(receive_index < node_offsets[node+1]);
          int_int_unordered_map_insert(receive_offsets, node, receive_index);
          int_array_append(receive_nodes, receive_index);
        }
      }

      int_unordered_set_free(their_culled_node_sets[p]);
    }

#ifndef NDEBUG
    // Check the node offsets.
    for (int n = 0; n < mesh->num_nodes; ++n)
    {
      int* offset_p = int_int_unordered_map_get(receive_offsets, n);
      if (offset_p != NULL)
      {
        ASSERT(*offset_p == node_offsets[n+1] - 1);
      }
    }
#endif

    int_int_unordered_map_free(receive_offsets);

    kd_tree_free(node_tree);
    int_array_free(my_node_indices);
    point_array_free(my_nodes);

    // Add the maps to the exchanger.
    exchanger_set_sends(ex, send_map);
    exchanger_set_receives(ex, receive_map);
    int_ptr_unordered_map_free(send_map);
    int_ptr_unordered_map_free(receive_map);
  }

  // Clean up.
  int_ptr_unordered_map_free(neighbor_neighbor_nodes);
  int_array_free(all_neighbors_of_neighbors);

#endif
  return ex;
}

