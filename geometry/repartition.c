// Copyright (c) 2012-2014, Jeffrey N. Johnson
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without 
// modification, are permitted provided that the following conditions are met:
// 
// 1. Redistributions of source code must retain the above copyright notice, this 
// list of conditions and the following disclaimer.
// 
// 2. Redistributions in binary form must reproduce the above copyright notice, 
// this list of conditions and the following disclaimer in the documentation 
// and/or other materials provided with the distribution.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE 
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include "geometry/repartition.h"

#if POLYMEC_HAVE_MPI
#include "core/unordered_set.h"
#include "ptscotch.h"

// This helper partitions a local graph, creating and returning a 
// local partition vector.
static SCOTCH_Num* repartition_graph(adj_graph_t* local_graph, 
                                     int num_ghost_vertices,
                                     int* weights,
                                     real_t imbalance_tol,
                                     exchanger_t* local_graph_ex)
{
  int nprocs, rank;
  MPI_Comm comm = adj_graph_comm(local_graph);
  MPI_Comm_size(comm, &nprocs);
  MPI_Comm_rank(comm, &rank);
  int num_vertices = adj_graph_num_vertices(local_graph);

  // Extract the adjacency information.
  SCOTCH_Num* xadj = malloc(sizeof(SCOTCH_Num) * (num_vertices+1));
  int* edge_offsets = adj_graph_edge_offsets(local_graph);
  for (int i = 0; i <= num_vertices; ++i)
    xadj[i] = (SCOTCH_Num)edge_offsets[i];
  SCOTCH_Num num_arcs = xadj[num_vertices];
  SCOTCH_Num* adj = malloc(sizeof(SCOTCH_Num) * num_arcs);
  int* edges = adj_graph_adjacency(local_graph);
  for (int i = 0; i < xadj[num_vertices]; ++i)
    adj[i] = (SCOTCH_Num)edges[i];

  // Replace the ghost entries in adj with global indices.
  index_t* vtx_dist = adj_graph_vertex_dist(local_graph);
  {
    index_t* global_indices = polymec_malloc(sizeof(SCOTCH_Num) * (num_vertices + num_ghost_vertices));
    for (int i = 0; i < num_vertices; ++i)
      global_indices[i] = (index_t)(vtx_dist[rank] + i);
    exchanger_exchange(local_graph_ex, global_indices, 1, 0, MPI_UINT64_T);
    for (int i = 0; i < num_arcs; ++i)
      adj[i] = global_indices[adj[i]];
    polymec_free(global_indices);
  }

  // Build a distributed graph.
  SCOTCH_Num* vtx_weights = NULL;
  if (weights != NULL)
  {
    vtx_weights = polymec_malloc(sizeof(SCOTCH_Num) * num_vertices);
    for (int i = 0; i < num_vertices; ++i)
      vtx_weights[i] = (SCOTCH_Num)weights[i];
  }
  SCOTCH_Dgraph dist_graph;
  SCOTCH_dgraphInit(&dist_graph, comm);
  SCOTCH_dgraphBuild(&dist_graph, 0, (SCOTCH_Num)num_vertices, (SCOTCH_Num)num_vertices,
                     xadj, NULL, vtx_weights, NULL, num_arcs, num_arcs,
                     adj, NULL, NULL);

  // Generate the local partition vector by mapping the distributed graph.
  SCOTCH_Num* local_partition = polymec_malloc(sizeof(SCOTCH_Num) * num_vertices);
  {
    SCOTCH_Arch arch;
    SCOTCH_archInit(&arch);
    if (vtx_weights == NULL)
      SCOTCH_archCmplt(&arch, num_vertices);
    else
      SCOTCH_archCmpltw(&arch, num_vertices, vtx_weights);
    SCOTCH_Strat strategy;
    SCOTCH_stratInit(&strategy);
    SCOTCH_Num strat_flags = SCOTCH_STRATDEFAULT;
    SCOTCH_stratDgraphMapBuild(&strategy, strat_flags, nprocs, nprocs, (double)imbalance_tol);
    int result = SCOTCH_dgraphMap(&dist_graph, &arch, &strategy, local_partition);
    if (result != 0)
      polymec_error("Repartitioning failed.");
    SCOTCH_stratExit(&strategy);
    SCOTCH_archExit(&arch);
  }
  if (vtx_weights != NULL)
    polymec_free(vtx_weights);

//  // Produce a redistributed graph using the partition vector.
//  SCOTCH_Dgraph redist_graph;
//  SCOTCH_dgraphInit(&redist_graph, comm);
//  result = SCOTCH_dgraphRedist(&dist_graph, partition, NULL, 0, 0, &redist_graph); 

  // Return the local partition vector.
  return local_partition;

#if 0
  // Assemble a global partition vector.
  index_t global_num_points = vtx_dist[nprocs];
  SCOTCH_Num* global_partition = polymec_malloc(sizeof(SCOTCH_Num) * global_num_points);
  int recv_counts[nprocs], offsets[nprocs];
  for (int p = 0; p < nprocs; ++p)
  {
    recv_counts[p] = (int)(vtx_dist[p+1] - vtx_dist[p]);
    offsets[p] = (int)vtx_dist[p];
  }
  MPI_Allgatherv(local_partition, num_vertices, MPI_UINT64_T, global_partition, 
                 recv_counts, offsets, MPI_UINT64_T, comm);
  polymec_free(local_partition);

  // Clean up.
  SCOTCH_dgraphExit(&dist_graph);

printf("P = [");
for (int i = 0; i < global_num_points; ++i)
printf("%d ", global_partition[i]);
printf("]\n");

  polymec_free(global_partition);
#endif
}

// This helper sets up the exchanger ex so that it can migrate data from the 
// current process according to the local partition vector.
static exchanger_t* create_migrator(MPI_Comm comm,
                                    SCOTCH_Num* local_partition,
                                    int num_vertices)
{
  exchanger_t* migrator = exchanger_new(comm);
  
  // Tally up the number of vertices we're going to send to every other process, 
  // including those that are staying on our own.
  int nprocs, rank;
  MPI_Comm_size(comm, &nprocs);
  MPI_Comm_rank(comm, &rank);
  int num_vertices_to_send[nprocs];
  memset(num_vertices_to_send, 0, sizeof(int) * nprocs);
  for (int v = 0; v < num_vertices; ++v)
    num_vertices_to_send[local_partition[v]]++;

  // Get the number of vertices we're going to receive from every other process.
  int num_vertices_to_receive[nprocs];
  MPI_Alltoall(num_vertices_to_send, 1, MPI_INT, 
               num_vertices_to_receive, 1, MPI_INT, comm);

  // Send and receive the actual vertices.
  int num_requests = 0;
  for (int p = 0; p < nprocs; ++p)
  {
    if ((rank != p) && (num_vertices_to_receive[p] > 0)) ++num_requests;
    if ((rank != p) && (num_vertices_to_send[p] > 0)) ++num_requests;
  }
  MPI_Request requests[num_requests];
  int r = 0;
  int** receive_vertices = polymec_malloc(sizeof(int*) * nprocs);
  memset(receive_vertices, 0, sizeof(int*) * nprocs);
  for (int p = 0; p < nprocs; ++p)
  {
    if (num_vertices_to_receive[p] > 0)
    {
      receive_vertices[p] = polymec_malloc(sizeof(int) * num_vertices_to_receive[p]);
      if (p != rank)
        MPI_Irecv(receive_vertices, num_vertices_to_receive[p], MPI_INT, p, p, comm, &requests[r++]);
    }
    if (num_vertices_to_send[p] > 0)
    {
      int send_vertices[num_vertices_to_send[p]], s = 0;
      for (int v = 0; v < num_vertices; ++v)
      {
        if (local_partition[v] == p)
          send_vertices[s++] = v;
      }
      exchanger_set_send(migrator, p, send_vertices, num_vertices_to_send[p], true);
      if (p != rank)
        MPI_Isend(send_vertices, num_vertices_to_receive[p], MPI_INT, p, p, comm, &requests[r++]);
    }
  }
  // Do any local copies.
  if (num_vertices_to_receive[rank] > 0)
  {
    ASSERT(num_vertices_to_receive[rank] == num_vertices_to_send[rank]);
    int send_vertices[num_vertices_to_send[rank]], s = 0;
    for (int v = 0; v < num_vertices; ++v)
    {
      if (local_partition[v] == rank)
        send_vertices[s++] = v;
    }
    memcpy(receive_vertices[rank], send_vertices, sizeof(int) * num_vertices_to_receive[rank]);
  }

  // Wait for exchanges to finish.
  MPI_Status statuses[num_requests];
  MPI_Waitall(num_requests, requests, statuses);

  // Now register all the vertices we're receiving with the migrator.
  r = 0;
  for (int p = 0; p < nprocs; ++p)
  {
    if (num_vertices_to_receive[p] > 0)
    {
      ASSERT(receive_vertices[p] != NULL)
      exchanger_set_receive(migrator, p, receive_vertices[p], num_vertices_to_receive[p], true);
      polymec_free(receive_vertices[p]);
    }
  }
  polymec_free(receive_vertices);
  exchanger_fprintf(migrator, stdout);

  return migrator;
}

static void point_cloud_migrate(point_cloud_t** cloud, 
                                adj_graph_t* local_graph, 
                                exchanger_t* migrator)
{
}

static void mesh_migrate(mesh_t** mesh, 
                         adj_graph_t* local_graph, 
                         exchanger_t* migrator)
{
  mesh_t* m = *mesh;

  // Post receives for buffer sizes.
  int num_receives = exchanger_num_receives(migrator);
  int num_sends = exchanger_num_sends(migrator);
  int receive_buffer_sizes[num_receives], receive_procs[num_receives];
  int pos = 0, proc, num_indices, *indices, i_req = 0;
  MPI_Request requests[num_receives + num_sends];
  while (exchanger_next_send(migrator, &pos, &proc, &indices, &num_indices))
  {
    receive_procs[i_req] = proc;
    MPI_Irecv(&receive_buffer_sizes[i_req], 1, MPI_INT, proc, 0, m->comm, &requests[i_req]);
    ++i_req;
  }

  // Build meshes to send to other processes.
  serializer_t* ser = mesh_serializer();
  int_unordered_set_t* cells = int_unordered_set_new();
  int_unordered_set_t* faces = int_unordered_set_new();
  int_unordered_set_t* nodes = int_unordered_set_new();
  byte_array_t* send_buffers[num_sends];
  int send_procs[num_sends];
  pos = 0;
  while (exchanger_next_send(migrator, &pos, &proc, &indices, &num_indices))
  {
    send_procs[i_req-num_receives] = proc;
    byte_array_t* bytes = byte_array_new();

    int_unordered_set_clear(cells);
    int_unordered_set_clear(faces);
    int_unordered_set_clear(nodes);

    // Make a set of these indices and count up the various thingies.
    int num_ghost_cells = 0; 
    for (int j = 0; j < num_indices; ++j)
    {
      int cell = indices[j];
      int_unordered_set_insert(cells, cell);
      int pos = 0, face;
      while (mesh_cell_next_face(m, cell, &pos, &face))
      {
        int_unordered_set_insert(faces, face);
        int pos1 = 0, node;
        while (mesh_face_next_node(m, face, &pos1, &node))
          int_unordered_set_insert(nodes, node);
      }
    }
    int num_cells = cells->size;
    int num_faces = faces->size;
    int num_nodes = nodes->size;

    // Create the mesh to send. We encode the indices as follows: 
    // All internal cells have their original local indices. A ghost cell,
    // on the other hand, has index -G - 2, where G is its global index.
    // Face and node indices are local and used only to express connectivity.
    mesh_t* submesh = mesh_new(m->comm, num_cells, num_ghost_cells, num_faces, num_nodes);
    // FIXME

    // Serialize and send the buffer size.
    size_t offset = 0;
    serializer_write(ser, submesh, bytes, &offset);
    MPI_Isend(&bytes->size, 1, MPI_INT, proc, 0, m->comm, &requests[i_req]);

    // Clean up.
    mesh_free(submesh);
    send_buffers[i_req - num_receives] = bytes;
    ++i_req;
  }
  ASSERT(i_req == num_sends + num_receives);

  // Preliminary clean up.
  int_unordered_set_free(cells);
  int_unordered_set_free(faces);
  int_unordered_set_free(nodes);

  // Wait for the buffer sizes to be transmitted.
  MPI_Status statuses[num_receives + num_sends];
  MPI_Waitall(num_receives + num_sends, requests, statuses);

  // Post receives for the actual messages.
  byte_array_t* receive_buffers[num_receives];
  for (int i = 0; i < num_receives; ++i)
  {
    receive_buffers[i] = byte_array_new();
    byte_array_resize(receive_buffers[i], receive_buffer_sizes[i]);
    MPI_Irecv(receive_buffers[i]->data, receive_buffer_sizes[i], MPI_BYTE, receive_procs[i], 0, m->comm, &requests[i]);
  }

  // Send the actual meshes and wait for receipt.
  for (int i = 0; i < num_sends; ++i)
    MPI_Isend(send_buffers[i]->data, send_buffers[i]->size, MPI_BYTE, send_procs[i], 0, m->comm, &requests[num_receives + i]);
  MPI_Waitall(num_receives + num_sends, requests, statuses);

  // Unpack the meshes.
  mesh_t* submeshes[num_receives];
  for (int i = 0; i < num_receives; ++i)
  {
    size_t offset = 0;
    submeshes[i] = serializer_read(ser, receive_buffers[i], &offset);
  }

  // Fuse 'em into a single mesh.
  // FIXME

  // Clean up all the stuff from the exchange.
  ser = NULL;
  for (int i = 0; i < num_receives; ++i)
    byte_array_free(receive_buffers[i]);
  polymec_free(receive_buffers);
  for (int i = 0; i < num_sends; ++i)
    byte_array_free(send_buffers[i]);
  polymec_free(send_buffers);
}

#endif

exchanger_t* repartition_point_cloud(point_cloud_t** cloud, int* weights, real_t imbalance_tol)
{
  ASSERT(imbalance_tol > 0.0);
  ASSERT(imbalance_tol <= 1.0);
  point_cloud_t* cl = *cloud;

#if POLYMEC_HAVE_MPI
  _Static_assert(sizeof(SCOTCH_Num) == sizeof(index_t), "SCOTCH_Num must be 64-bit.");

  int nprocs, rank;
  MPI_Comm_size(cl->comm, &nprocs);
  MPI_Comm_rank(cl->comm, &rank);

  // On a single process, repartitioning has no meaning.
  if (nprocs == 1)
    return exchanger_new(cl->comm);

  // Generate a local adjacency graph for the point cloud.
  adj_graph_t* local_graph = graph_from_point_cloud(cl);

  // Get the exchanger for the point cloud.
  exchanger_t* cloud_ex = point_cloud_exchanger(cl);

  // Map the graph to the different domains, producing a local partition vector.
  SCOTCH_Num* local_partition = repartition_graph(local_graph, cl->num_ghost_points, weights, imbalance_tol, cloud_ex);

  // Set up an exchanger to migrate field data.
  int num_vertices = adj_graph_num_vertices(local_graph);
  exchanger_t* migrator = create_migrator(cl->comm, local_partition, num_vertices);

  // Migrate the point cloud.
  point_cloud_migrate(cloud, local_graph, migrator);

  // Clean up.
  adj_graph_free(local_graph);
  polymec_free(local_partition);

  // Return the migrator.
  return (migrator == NULL) ? exchanger_new(cl->comm) : migrator;
#else
  return exchanger_new(cl->comm);
#endif
}

exchanger_t* repartition_mesh(mesh_t** mesh, int* weights, real_t imbalance_tol)
{
  ASSERT(imbalance_tol > 0.0);
  ASSERT(imbalance_tol <= 1.0);
  mesh_t* m = *mesh;

#if POLYMEC_HAVE_MPI
  _Static_assert(sizeof(SCOTCH_Num) == sizeof(index_t), "SCOTCH_Num must be 64-bit.");

  int nprocs, rank;
  MPI_Comm_size(m->comm, &nprocs);
  MPI_Comm_rank(m->comm, &rank);

  // On a single process, repartitioning has no meaning.
  if (nprocs == 1)
    return exchanger_new(m->comm);

  // Generate a local adjacency graph for the mesh.
  adj_graph_t* local_graph = graph_from_mesh_cells(m);

  // Get the exchanger for the mesh.
  exchanger_t* mesh_ex = mesh_exchanger(m);

  // Map the graph to the different domains, producing a local partition vector.
  SCOTCH_Num* local_partition = repartition_graph(local_graph, m->num_ghost_cells, weights, imbalance_tol, mesh_ex);

  // Set up an exchanger to migrate field data.
  int num_vertices = adj_graph_num_vertices(local_graph);
  exchanger_t* migrator = create_migrator(m->comm, local_partition, num_vertices);

  // Migrate the mesh.
  mesh_migrate(mesh, local_graph, migrator);

  // Clean up.
  adj_graph_free(local_graph);
  polymec_free(local_partition);

  // Return the migrator.
  return migrator;
#else
  return exchanger_new(m->comm);
#endif
}

