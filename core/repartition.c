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

#include "core/repartition.h"

#if POLYMEC_HAVE_MPI
#include "core/unordered_set.h"
#include "core/kd_tree.h"
#include "ptscotch.h"

// This helper partitions a (serial) global graph, creating and returning a 
// global partition vector. It should only be called on rank 0.
static SCOTCH_Num* partition_graph(adj_graph_t* global_graph, 
                                   MPI_Comm comm,
                                   int* weights,
                                   real_t imbalance_tol)
{
  ASSERT(adj_graph_comm(global_graph) == MPI_COMM_SELF);

  int nprocs, rank;
  MPI_Comm_size(comm, &nprocs);
  MPI_Comm_rank(comm, &rank);
  ASSERT(rank == 0);

  SCOTCH_Dgraph dist_graph;
  SCOTCH_Num *vtx_weights = NULL, *xadj = NULL, *adj = NULL;
  int num_global_vertices = 0;
  if (rank == 0)
  {
    num_global_vertices = adj_graph_num_vertices(global_graph);

    // Extract the adjacency information.
    xadj = polymec_malloc(sizeof(SCOTCH_Num) * (num_global_vertices+1));
    int* edge_offsets = adj_graph_edge_offsets(global_graph);
    for (int i = 0; i <= num_global_vertices; ++i)
      xadj[i] = (SCOTCH_Num)edge_offsets[i];
    SCOTCH_Num num_arcs = xadj[num_global_vertices];
    adj = polymec_malloc(sizeof(SCOTCH_Num) * num_arcs);
    int* edges = adj_graph_adjacency(global_graph);
    for (int i = 0; i < xadj[num_global_vertices]; ++i)
      adj[i] = (SCOTCH_Num)edges[i];

    // Build a graph on rank 0.
    SCOTCH_dgraphInit(&dist_graph, MPI_COMM_SELF);
    if (weights != NULL)
    {
      vtx_weights = polymec_malloc(sizeof(SCOTCH_Num) * num_global_vertices);
      for (int i = 0; i < num_global_vertices; ++i)
        vtx_weights[i] = (SCOTCH_Num)weights[i];
    }
    SCOTCH_dgraphBuild(&dist_graph, 0, (SCOTCH_Num)num_global_vertices, (SCOTCH_Num)num_global_vertices,
        xadj, NULL, vtx_weights, NULL, num_arcs, num_arcs,
        adj, NULL, NULL);
  }

  // Generate the global partition vector by scattering the global graph.
  SCOTCH_Num* global_partition = NULL;
  if (rank == 0)
  {
    global_partition = polymec_malloc(sizeof(SCOTCH_Num) * num_global_vertices);
    SCOTCH_Strat strategy;
    SCOTCH_stratInit(&strategy);
    SCOTCH_Num strat_flags = SCOTCH_STRATDEFAULT;
    SCOTCH_stratDgraphMapBuild(&strategy, strat_flags, nprocs, nprocs, (double)imbalance_tol);
    int result = SCOTCH_dgraphPart(&dist_graph, nprocs, &strategy, global_partition);
    if (result != 0)
      polymec_error("Partitioning failed.");
    SCOTCH_dgraphExit(&dist_graph);
    SCOTCH_stratExit(&strategy);

    if (vtx_weights != NULL)
      polymec_free(vtx_weights);
    polymec_free(xadj);
    polymec_free(adj);
  }

  // Return the global partition vector.
  return global_partition;
}

// This helper repartitions a local graph, creating and returning a 
// local partition vector with destination ranks included for ghost cells.
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
  SCOTCH_Num* xadj = polymec_malloc(sizeof(SCOTCH_Num) * (num_vertices+1));
  int* edge_offsets = adj_graph_edge_offsets(local_graph);
  for (int i = 0; i <= num_vertices; ++i)
    xadj[i] = (SCOTCH_Num)edge_offsets[i];
  SCOTCH_Num num_arcs = xadj[num_vertices];
  SCOTCH_Num* adj = polymec_malloc(sizeof(SCOTCH_Num) * num_arcs);
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
  SCOTCH_Num* local_partition = polymec_malloc(sizeof(SCOTCH_Num) * (num_vertices + num_ghost_vertices));
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
    SCOTCH_dgraphExit(&dist_graph);
    SCOTCH_stratExit(&strategy);
    SCOTCH_archExit(&arch);
    polymec_free(adj);
    polymec_free(xadj);
  }
  if (vtx_weights != NULL)
    polymec_free(vtx_weights);

  // At this point, the local partition vector contains only the destination 
  // ranks of the local vertices. Now we talk amongst ourselves to figure 
  // out the destination ranks of the ghost vertices. This operation should 
  // preserve the ordering of the ghost vertices, too.
  exchanger_exchange(local_graph_ex, local_partition, 1, 0, MPI_UINT64_T);

//printf("%d: P = [", rank);
//for (int i = 0; i < num_vertices + num_ghost_vertices; ++i)
//printf("%d ", local_partition[i]);
//printf("]\n");

  // Return the local partition vector.
  return local_partition;
}

// This helper sets up the exchanger ex so that it can distribute data from the 
// root process (0) according to the global partition vector. The partition 
// vector is NULL on all nonzero ranks and defined on rank 0.
static exchanger_t* create_distributor(MPI_Comm comm, 
                                       SCOTCH_Num* global_partition,
                                       int num_global_vertices)
{
  int nprocs, rank;
  MPI_Comm_size(comm, &nprocs);
  MPI_Comm_rank(comm, &rank);
  ASSERT((rank == 0) || (global_partition == NULL));
  ASSERT((rank == 0) || (num_global_vertices == 0));

  exchanger_t* distributor = exchanger_new(comm);
  
  if (rank == 0)
  {
    // Get the number of vertices we're going to send to each other process.
    int num_vertices_to_send[nprocs];
    memset(num_vertices_to_send, 0, sizeof(int) * nprocs);
    for (int v = 0; v < num_global_vertices; ++v)
      num_vertices_to_send[global_partition[v]]++;

    // Send this number.
    int n;
    MPI_Scatter(num_vertices_to_send, 1, MPI_INT, &n, 1, MPI_INT, 0, comm);

    // Now send the vertices to each process and register these sends
    // with the distributor.
    for (int p = 1; p < nprocs; ++p)
    {
      ASSERT(num_vertices_to_send[p] > 0);
      int vertices[num_vertices_to_send[p]], k = 0, p_tag = p;
      for (int i = 0; i < num_global_vertices; ++i)
      {
        if (global_partition[i] == p)
          vertices[k++] = i;
      }
      MPI_Send(vertices, num_vertices_to_send[p], MPI_INT, p, p_tag, comm);
      exchanger_set_send(distributor, p, vertices, num_vertices_to_send[p], true);
    }

    // Figure out the local vertices for rank 0.
    int num_local_vertices = num_vertices_to_send[0];
    int local_vertices[num_local_vertices], k = 0;
    for (int i = 0; i < num_global_vertices; ++i)
    {
      if (global_partition[i] == 0)
        local_vertices[k++] = i;
    }
  }
  else
  {
    // Get the number of vertices we will receive from rank 0.
    int num_local_vertices;
    MPI_Scatter(NULL, 1, MPI_INT, &num_local_vertices, 1, MPI_INT, 0, comm);

    // Now get the vertices.
    int local_vertices[num_local_vertices], p_tag = rank;
    MPI_Status status;
    MPI_Recv(local_vertices, num_local_vertices, MPI_INT, 0, p_tag, comm, &status);

    // Now register all the vertices we're receiving with the migrator.
    ASSERT(num_local_vertices > 0);
    exchanger_set_receive(distributor, 0, local_vertices, num_local_vertices, true);
  }

  return distributor;
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
      {
        MPI_Irecv(receive_vertices, num_vertices_to_receive[p], MPI_INT, p, p, comm, &requests[r++]);
      }
    }
    if (num_vertices_to_send[p] > 0)
    {
      int send_vertices[num_vertices_to_send[p]], s = 0;
      for (int v = 0; v < num_vertices; ++v)
      {
        if (local_partition[v] == p)
          send_vertices[s++] = v;
      }
      if (p != rank)
      {
        exchanger_set_send(migrator, p, send_vertices, num_vertices_to_send[p], true);
        MPI_Isend(send_vertices, num_vertices_to_receive[p], MPI_INT, p, p, comm, &requests[r++]);
      }
      else
        memcpy(receive_vertices[rank], send_vertices, sizeof(int) * num_vertices_to_receive[rank]);
    }
  }
  ASSERT(r == num_requests);

  // Wait for exchanges to finish.
  MPI_Status statuses[num_requests];
  MPI_Waitall(num_requests, requests, statuses);

  // Now register all the vertices we're receiving with the migrator.
  r = 0;
  for (int p = 0; p < nprocs; ++p)
  {
    if ((rank != p) && (num_vertices_to_receive[p] > 0))
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

// This helper constructs and returns a point cloud from the given points in 
// the given point cloud.
static point_cloud_t* create_subcloud(point_cloud_t* cloud, int* indices, int num_indices)
{
  return NULL;
}

static void point_cloud_distribute(point_cloud_t** cloud, 
                                   MPI_Comm comm,
                                   adj_graph_t* global_graph, 
                                   SCOTCH_Num* global_partition)
{
}

static void point_cloud_migrate(point_cloud_t** cloud, 
                                adj_graph_t* local_graph, 
                                exchanger_t* migrator)
{
}

// This helper constructs and returns a mesh from the cells with the given
// indices in the given mesh. The submesh is valid with the following 
// exceptions:
// 1. The indices of ghost cells referenced in submesh->face_cells are 
//    replaced with destination process ranks.
// 2. The geometry for the submesh is not computed.
static mesh_t* create_submesh(MPI_Comm comm, mesh_t* mesh, 
                              SCOTCH_Num* partition, index_t* vtx_dist, 
                              int* indices, int num_indices)
{
  // Make a set of cells for querying membership in this submesh.
  int_unordered_set_t* cell_set = int_unordered_set_new();
  for (int i = 0; i < num_indices; ++i)
    int_unordered_set_insert(cell_set, indices[i]);

  // Count unique mesh elements.
  int num_cells = num_indices, num_ghost_cells = 0;
  int_unordered_set_t* face_indices = int_unordered_set_new();
  int_unordered_set_t* node_indices = int_unordered_set_new();
  for (int i = 0; i < num_indices; ++i)
  {
    int cell = indices[i], fpos = 0, face;
    while (mesh_cell_next_face(mesh, cell, &fpos, &face))
    {
      int opp_cell = mesh_face_opp_cell(mesh, face, cell);
      //if (opp_cell >= mesh->num_cells)
      if ((opp_cell != -1) && !int_unordered_set_contains(cell_set, opp_cell))
        ++num_ghost_cells;

      int_unordered_set_insert(face_indices, face);
      int npos = 0, node;
      while (mesh_face_next_node(mesh, face, &npos, &node))
        int_unordered_set_insert(node_indices, node);
    }
  }
  int_unordered_set_free(cell_set);

  // Create the submesh container.
  int num_faces = face_indices->size, num_nodes = node_indices->size;
  mesh_t* submesh = mesh_new(comm, num_cells, num_ghost_cells, num_faces, num_nodes);

  // Create mappings for faces and nodes.
  int* face_map = polymec_malloc(sizeof(int) * num_faces);
  {
    int fpos = 0, face, f = 0;
    while (int_unordered_set_next(face_indices, &fpos, &face))
      face_map[f++] = face;
  }
  int_unordered_set_free(face_indices);
  int* node_map = polymec_malloc(sizeof(int) * num_nodes);
  {
    int npos = 0, node, n = 0;
    while (int_unordered_set_next(node_indices, &npos, &node))
      node_map[n++] = node;
  }
  int_unordered_set_free(node_indices);

  // Allocate cell faces and face nodes (and generate inverse maps).
  int_int_unordered_map_t* inverse_cell_map = int_int_unordered_map_new();
  submesh->cell_face_offsets[0] = 0;
  for (int c = 0; c < submesh->num_cells; ++c)
  {
    int_int_unordered_map_insert(inverse_cell_map, indices[c], c);
    int num_cell_faces = mesh->cell_face_offsets[indices[c]+1] - mesh->cell_face_offsets[indices[c]];
    submesh->cell_face_offsets[c+1] = submesh->cell_face_offsets[c] + num_cell_faces;
  }
  int_int_unordered_map_t* inverse_face_map = int_int_unordered_map_new();
  submesh->face_node_offsets[0] = 0;
  for (int f = 0; f < submesh->num_faces; ++f)
  {
    int_int_unordered_map_insert(inverse_face_map, face_map[f], f);
    int num_face_nodes = mesh->face_node_offsets[face_map[f]+1] - mesh->face_node_offsets[face_map[f]];
    submesh->face_node_offsets[f+1] = submesh->face_node_offsets[f] + num_face_nodes;
  }
  int_int_unordered_map_t* inverse_node_map = int_int_unordered_map_new();
  for (int n = 0; n < submesh->num_nodes; ++n)
    int_int_unordered_map_insert(inverse_node_map, node_map[n], n);
  mesh_reserve_connectivity_storage(submesh);

  // Copy cell faces.
  for (int c = 0; c < submesh->num_cells; ++c)
  {
    int num_cell_faces = submesh->cell_face_offsets[c+1] - submesh->cell_face_offsets[c];
    for (int f = 0; f < num_cell_faces; ++f)
    {
      int face_index = mesh->cell_faces[mesh->cell_face_offsets[indices[c]]+f];
      bool flipped = false;
      if (face_index < 0)
      {
        flipped = true;
        face_index = ~face_index;
      }
      int subface_index = *int_int_unordered_map_get(inverse_face_map, face_index);
      if (flipped)
        subface_index = ~subface_index;
      submesh->cell_faces[submesh->cell_face_offsets[c]+f] = subface_index;
    }
  }

  // Copy face cells and face nodes.
  for (int f = 0; f < submesh->num_faces; ++f)
  {
    // Identify the cells attached to the face and construct them within the submesh.
    int orig_mesh_face = face_map[f];
    int* cell_p = int_int_unordered_map_get(inverse_cell_map, mesh->face_cells[2*orig_mesh_face]);
    int* opp_cell_p = int_int_unordered_map_get(inverse_cell_map, mesh->face_cells[2*orig_mesh_face+1]);
    ASSERT((cell_p != NULL) || (opp_cell_p != NULL));
    int this_cell = -1, that_cell = -1;
    if ((cell_p != NULL) && (opp_cell_p != NULL))
    {
      this_cell = *cell_p;
      that_cell = *opp_cell_p;
    }
    else
    {
      // One of the cells attached to this face is a ghost cell.
      int ghost_cell;
      if (cell_p != NULL)
      {
        this_cell = *cell_p;
        ghost_cell = mesh->face_cells[2*orig_mesh_face+1];
      }
      else
      {
        this_cell = *opp_cell_p;
        ghost_cell = mesh->face_cells[2*orig_mesh_face];
      }

      // We encode the destination process rank in the ghost cell.
      if (ghost_cell != -1)
      {
        that_cell = -partition[ghost_cell] - 2;
      }
    }
    submesh->face_cells[2*f] = this_cell;
    submesh->face_cells[2*f+1] = that_cell;

    // Copy over the nodes for this face.
    int num_face_nodes = submesh->face_node_offsets[f+1] - submesh->face_node_offsets[f];
    for (int n = 0; n < num_face_nodes; ++n)
    {
      int subnode = *int_int_unordered_map_get(inverse_node_map, mesh->face_nodes[mesh->face_node_offsets[orig_mesh_face]+n]);
      submesh->face_nodes[submesh->face_node_offsets[f]+n] = subnode;
    }
  }

  // Copy node positions.
  for (int n = 0; n < submesh->num_nodes; ++n)
  {
    int orig_mesh_node = node_map[n];
    submesh->nodes[n] = mesh->nodes[orig_mesh_node];
  }

  // NOTE: we don't need to construct edges or compute geometry, since 
  // NOTE: this can be done on the "far end."

  // Clean up.
  int_int_unordered_map_free(inverse_cell_map);
  int_int_unordered_map_free(inverse_face_map);
  int_int_unordered_map_free(inverse_node_map);
  polymec_free(face_map);
  polymec_free(node_map);

  return submesh;
}

static void mesh_distribute(mesh_t** mesh, 
                            MPI_Comm comm,
                            adj_graph_t* global_graph, 
                            SCOTCH_Num* global_partition)
{
  int nprocs, rank;
  MPI_Comm_size(comm, &nprocs);
  MPI_Comm_rank(comm, &rank);

  // Make sure we're all here.
  MPI_Barrier(comm);

  mesh_t* global_mesh = *mesh;
  mesh_t* local_mesh = NULL;
  uint64_t vtx_dist[nprocs+1];
  if (rank == 0)
  {
    // Take stock of how many cells we'll have per process.
    int num_cells[nprocs];
    memset(num_cells, 0, sizeof(int) * nprocs);
    for (int i = 0; i < global_mesh->num_cells; ++i)
      num_cells[global_partition[i]]++;

    // Construct the distribution of vertices for the partitioning.
    vtx_dist[0] = 0;
    for (int p = 0; p < nprocs; ++p)
      vtx_dist[p+1] = vtx_dist[p] + num_cells[p];

    // Carve out the portion of the mesh that will stick around on process 0.
    {
      int indices[num_cells[0]], k = 0;
      for (int i = 0; i < global_mesh->num_cells; ++i)
      {
        if (global_partition[i] == rank)
          indices[k++] = i;
      }
      local_mesh = create_submesh(comm, global_mesh, global_partition, NULL, indices, num_cells[0]);

      // Construct edges, geometry, since this hasn't been done yet.
      mesh_construct_edges(local_mesh);
      mesh_compute_geometry(local_mesh);
    }

    // Now do the other processes.
    serializer_t* ser = mesh_serializer();
    byte_array_t* bytes = byte_array_new();
    for (int p = 1; p < nprocs; ++p)
    {
      // Share the vertex distribution.
      MPI_Send(vtx_dist, nprocs+1, MPI_UINT64_T, p, p, comm);

      // Create the pth submesh.
      int indices[num_cells[p]], k = 0;
      for (int i = 0; i < global_mesh->num_cells; ++i)
      {
        if (global_partition[i] == p)
          indices[k++] = i;
      }
      mesh_t* p_mesh = create_submesh(comm, global_mesh, global_partition, NULL, indices, num_cells[p]);

      // Serialize it and send its size (and it) to process p.
      size_t offset = 0;
      serializer_write(ser, p_mesh, bytes, &offset);
      MPI_Send(&bytes->size, 1, MPI_INT, p, p, comm);
      MPI_Send(bytes->data, bytes->size, MPI_BYTE, p, p, comm);

      // Clean up.
      byte_array_clear(bytes);
      mesh_free(p_mesh);
    }
    ser = NULL;
    byte_array_free(bytes);
  }
  else
  {
    // Receive the vertex distribution of the incoming mesh.
    MPI_Status status;
    MPI_Recv(vtx_dist, nprocs+1, MPI_UINT64_T, 0, rank, comm, &status);

    // Receive the size of the incoming mesh.
    int mesh_size;
    MPI_Recv(&mesh_size, 1, MPI_INT, 0, rank, comm, &status);

    // Now receive the mesh.
    byte_array_t* bytes = byte_array_new();
    byte_array_resize(bytes, mesh_size);

    MPI_Recv(bytes->data, mesh_size, MPI_BYTE, 0, rank, comm, &status);
    serializer_t* ser = mesh_serializer();
    size_t offset = 0;
    local_mesh = serializer_read(ser, bytes, &offset);
    
    byte_array_free(bytes);
    ser = NULL;
  }

  *mesh = local_mesh;

  // Clean up.
  if (global_mesh != NULL)
    mesh_free(global_mesh);

  // Now we create the exchanger using the encoded destination ranks.
  int_ptr_unordered_map_t* ghost_cell_indices = int_ptr_unordered_map_new();
  int num_ghosts = 0;
  for (int f = 0; f < local_mesh->num_faces; ++f)
  {
    if (local_mesh->face_cells[2*f+1] < -1) 
    {
      // Get the destination process for this ghost cell.
      int proc = -local_mesh->face_cells[2*f+1] - 2;
      ASSERT(proc >= 0);
      ASSERT(proc < nprocs);
      local_mesh->face_cells[2*f+1] = num_ghosts;

      if (!int_ptr_unordered_map_contains(ghost_cell_indices, proc))
        int_ptr_unordered_map_insert_with_v_dtor(ghost_cell_indices, proc, int_array_new(), DTOR(int_array_free));
      int_array_t* indices = *int_ptr_unordered_map_get(ghost_cell_indices, proc);
      int_array_append(indices, local_mesh->face_cells[2*f]);
      int_array_append(indices, num_ghosts++);
    }
  }
  ASSERT(num_ghosts == local_mesh->num_ghost_cells);
  exchanger_t* ex = mesh_exchanger(local_mesh);
  int pos = 0, proc;
  index_array_t* indices;
  while (int_ptr_unordered_map_next(ghost_cell_indices, &pos, &proc, (void**)&indices))
  {
    if (proc != rank) 
    {
      int send_indices[indices->size/2], recv_indices[indices->size/2];
      for (int i = 0; i < indices->size/2; ++i)
      {
        send_indices[i] = indices->data[2*i];
        recv_indices[i] = indices->data[2*i+1];
      }
      exchanger_set_send(ex, proc, send_indices, indices->size/2, true);
      exchanger_set_receive(ex, proc, recv_indices, indices->size/2, true);
    }
  }

  // Clean up again.
  int_ptr_unordered_map_free(ghost_cell_indices);
}

// This helper takes an array of submeshes and stitches them all together into 
// one single mesh (contiguous or not) on the current domain. The submeshes are 
// consumed in the process.
static mesh_t* fuse_submeshes(mesh_t** submeshes, 
                              int num_submeshes)
{
  ASSERT(num_submeshes > 0);

  // If we're only given 1 submesh, simply copy it.
  if (num_submeshes == 1)
    return submeshes[0];

  int rank;
  MPI_Comm_rank(submeshes[0]->comm, &rank);

  // First we traverse each of the submeshes and count up all the internal 
  // and ghost cells, faces, nodes.
  int num_cells = 0, num_ghost_cells = 0, num_faces = 0, num_nodes = 0, 
      submesh_cell_offsets[num_submeshes+1], submesh_face_offsets[num_submeshes+1], 
      submesh_node_offsets[num_submeshes+1];
  submesh_cell_offsets[0] = submesh_face_offsets[0] = submesh_node_offsets[0] = 0;
  for (int m = 0; m < num_submeshes; ++m)
  {
    mesh_t* submesh = submeshes[m];
    num_cells += submesh->num_cells;
    num_ghost_cells += submesh->num_ghost_cells;
    submesh_cell_offsets[m+1] = num_cells;
    num_faces += submesh->num_faces;
    submesh_face_offsets[m+1] = num_faces;
    num_nodes += submesh->num_nodes;
    submesh_node_offsets[m+1] = num_nodes;
  }
  printf("%d: Fusing %d submeshes totaling %d cells.", rank, num_submeshes, num_cells);

  // Next, we traverse these submeshes and construct sets of all the 
  // faces/nodes that make up the "seams" of the fused mesh--those whose 
  // would-be-ghost cells belong to other submeshes on this process. 
  int_unordered_set_t* seam_faces = int_unordered_set_new();
  int_unordered_set_t* seam_nodes = int_unordered_set_new();
  for (int m = 0; m < num_submeshes; ++m)
  {
    mesh_t* submesh = submeshes[m];
    for (int f = 0; f < submesh->num_faces; ++f)
    {
      if (submesh->face_cells[2*f+1] == -rank - 2) // belongs to this domain
      {
        int_unordered_set_insert(seam_faces, submesh_face_offsets[m] + f);
        int pos = 0, node;
        while (mesh_face_next_node(submesh, f, &pos, &node))
          int_unordered_set_insert(seam_nodes, submesh_node_offsets[m] + node);
      }
    }
  }
  num_ghost_cells -= seam_faces->size;

  // Construct mappings to remove duplicate faces by mapping redundant 
  // ones to "originals."
  int_int_unordered_map_t* dup_face_map = int_int_unordered_map_new();
  {
    // Make a kd-tree of face centers.
    point_t* xf = polymec_malloc(sizeof(point_t) * seam_faces->size);
    int pos = 0, i = 0, j;
    while (int_unordered_set_next(seam_faces, &pos, &j))
    {
      // Find the submesh that contains this face.
      int m = 0;
      while (j < submesh_face_offsets[m]) ++m;

      // Put the face center into this array.
      int f = j - submesh_face_offsets[m];
      xf[i] = submeshes[m]->face_centers[f];
      ++i;
    }
    kd_tree_t* face_tree = kd_tree_new(xf, seam_faces->size);
    polymec_free(xf);

    // Now examine each seam face and find the 2 nearest faces to it.
    // One will be the face itself, and the other will be a duplicate. 
    // Map the one with the higher index to the lower index.
    pos = 0;
    while (int_unordered_set_next(seam_faces, &pos, &j))
    {
      int nearest[2];

      // Get the face center we're considering (in the same way we did above).
      int m = 0;
      while (j < submesh_face_offsets[m]) ++m;
      int f = j - submesh_face_offsets[m];
      point_t* x = &submeshes[m]->face_centers[f];

      // Find the 2 nearest faces.
      kd_tree_nearest_n(face_tree, x, 2, nearest);

      // Merge the faces by mapping the one with the higher index to the 
      // lower index. 
      int index0 = submesh_face_offsets[m] + nearest[0];
      int index1 = submesh_face_offsets[m] + nearest[1];
      ASSERT(index0 != index1);
      int_int_unordered_map_insert(dup_face_map, MAX(index0, index1), MIN(index0, index1));

      // ALSO: Alter the submesh so that its interior ghost cells actually 
      // refer to the correct (flattened) interior cells. This may seem strange,
      // but it helps us to get the face->cell connectivity right later on.
      int j1;
      if (j == MIN(index0, index1))
        j1 = MAX(index0, index1);
      else
        j1 = MIN(index0, index1);
      int m1 = 0;
      while (j1 < submesh_face_offsets[m1]) ++m1;
      int f1 = j1 - submesh_face_offsets[m1];
      int flattened_cell1 = submesh_cell_offsets[m1] + submeshes[m1]->face_cells[2*f1];
      submeshes[m]->face_cells[2*f+1] = flattened_cell1;
    }

    // Clean up.
    kd_tree_free(face_tree);
  }

  // Do the same for duplicate nodes.
  int_int_unordered_map_t* dup_node_map = int_int_unordered_map_new();
  {
    // Make a kd-tree of node positions.
    point_t* xn = polymec_malloc(sizeof(point_t) * seam_nodes->size);
    int pos = 0, i = 0, j;
    while (int_unordered_set_next(seam_nodes, &pos, &j))
    {
      // Find the submesh that contains this face.
      int m = 0;
      while (j < submesh_node_offsets[m]) ++m;

      // Put the node position into this array.
      int n = j - submesh_node_offsets[m];
      xn[i] = submeshes[m]->nodes[n];
      ++i;
    }
    kd_tree_t* node_tree = kd_tree_new(xn, seam_nodes->size);
    polymec_free(xn);

    // Now examine each seam node and find the 2 nearest node to it.
    // One will be the node itself, and the other will be a duplicate. 
    // Map the one with the higher index to the lower index.
    pos = 0;
    while (int_unordered_set_next(seam_nodes, &pos, &j))
    {
      int nearest[2];

      // Get the node position we're considering (in the same way we did above).
      int m = 0;
      while (j < submesh_node_offsets[m]) ++m;
      int n = j - submesh_node_offsets[m];
      point_t* x = &submeshes[m]->nodes[n];

      // Find the 2 nearest nodes.
      kd_tree_nearest_n(node_tree, x, 2, nearest);

      // Merge the nodes by mapping the one with the higher index to the 
      // lower index.
      int index0 = submesh_node_offsets[m] + nearest[0];
      int index1 = submesh_node_offsets[m] + nearest[1];
      int_int_unordered_map_insert(dup_node_map, MAX(index0, index1), MIN(index0, index1));
    }

    // Clean up.
    kd_tree_free(node_tree);
  }

  // Reduce the number of faces and nodes in the fused mesh by the ones that 
  // have been merged to others.
  num_faces -= seam_faces->size;
  num_nodes -= seam_nodes->size;

  // We're through with the seam faces/nodes.
  int_unordered_set_free(seam_faces);
  int_unordered_set_free(seam_nodes);

  // Now we create the fused mesh and fill it with the contents of the submeshes.
  mesh_t* fused_mesh = mesh_new(submeshes[0]->comm, num_cells, num_ghost_cells,
                                num_faces, num_nodes);

  // Allocate storage for cell faces and face nodes.
  fused_mesh->cell_face_offsets[0] = 0;
  int cell = 0, face = 0;
  for (int m = 0; m < num_submeshes; ++m)
  {
    mesh_t* submesh = submeshes[m];
    for (int c = 0; c < submesh->num_cells; ++c, ++cell)
    {
      int num_cell_faces = submesh->cell_face_offsets[c+1] - submesh->cell_face_offsets[c];
      fused_mesh->cell_face_offsets[cell+1] = fused_mesh->cell_face_offsets[cell] + num_cell_faces;
    }
    for (int f = 0; f < submesh->num_faces; ++f)
    {
      int flattened_face = submesh_face_offsets[m] + f;
      if (!int_int_unordered_map_contains(dup_face_map, flattened_face))
      {
        int num_face_nodes = submesh->face_node_offsets[f+1] - submesh->face_node_offsets[f];
        fused_mesh->face_node_offsets[face+1] = fused_mesh->face_node_offsets[face] + num_face_nodes;
        ++face;
      }
    }
  }
  mesh_reserve_connectivity_storage(fused_mesh);

  // Copy cell faces.
  cell = 0;
  for (int m = 0; m < num_submeshes; ++m)
  {
    mesh_t* submesh = submeshes[m];
    for (int c = 0; c < submesh->num_cells; ++c, ++cell)
    {
      int num_cell_faces = submesh->cell_face_offsets[c+1] - submesh->cell_face_offsets[c];
      for (int f = 0; f < num_cell_faces; ++f)
      {
        int subface_index = submesh->cell_faces[submesh->cell_face_offsets[c]] + f;
        bool flipped = false;
        if (subface_index < 0)
        {
          flipped = true;
          subface_index = ~subface_index;
        }
        int flattened_face = submesh_face_offsets[m] + subface_index;
        int* face_p = int_int_unordered_map_get(dup_face_map, flattened_face);
        int face_index = (face_p != NULL) ? *face_p : flattened_face;
        if (flipped)
          face_index = ~face_index;
        fused_mesh->cell_faces[fused_mesh->cell_face_offsets[cell]+f] = face_index;
      }
    }
  }

  // Copy face cells and face nodes.
  int node = 0;
  face = 0;
  for (int m = 0; m < num_submeshes; ++m)
  {
    mesh_t* submesh = submeshes[m];
    for (int f = 0; f < submesh->num_faces; ++f)
    {
      int flattened_face = submesh_face_offsets[m] + f;
      if (!int_int_unordered_map_contains(dup_face_map, flattened_face))
      {
        // Set up the cells of the face. We can do this expediently now 
        // because of the way we altered the submeshes with interior ghost cells.
        int subcell_index = submesh->face_cells[2*f];
        int flattened_cell = submesh_cell_offsets[m] + subcell_index;
        int flattened_cell1 = submesh->face_cells[2*f+1]; // we stuck this here earlier!
        fused_mesh->face_cells[2*face] = flattened_cell;
        fused_mesh->face_cells[2*face+1] = flattened_cell1;

        // Set up the nodes of the face.
        int num_face_nodes = submesh->face_node_offsets[f+1] - submesh->face_node_offsets[f];
        for (int n = 0; n < num_face_nodes; ++n)
        {
          int subnode_index = submesh->face_nodes[submesh->face_node_offsets[f]] + n;
          int flattened_node = submesh_node_offsets[m] + subnode_index;
          int* node_p = int_int_unordered_map_get(dup_node_map, flattened_node);
          int node_index = (node_p != NULL) ? *node_p : flattened_node;
          fused_mesh->face_nodes[fused_mesh->face_node_offsets[face]+f] = node_index;
        }
        ++face;
      }
    }
  }

  // Copy node positions.
  node = 0;
  for (int m = 0; m < num_submeshes; ++m)
  {
    mesh_t* submesh = submeshes[m];
    for (int n = 0; n < submesh->num_nodes; ++n)
    {
      int flattened_node = submesh_node_offsets[m] + n;
      if (!int_int_unordered_map_contains(dup_node_map, flattened_node))
      {
        fused_mesh->nodes[node] = submesh->nodes[n];
        ++node;
      }
    }
  }

  // Consume the submeshes.
  for (int m = 0; m < num_submeshes; ++m)
    mesh_free(submeshes[m]);

  // Now fill the exchanger for the fused mesh with data.
  int_ptr_unordered_map_t* send_map = int_ptr_unordered_map_new(); 
  int_ptr_unordered_map_t* recv_map = int_ptr_unordered_map_new(); 
  int ghost_cell = fused_mesh->num_cells;
  for (int f = 0; f < fused_mesh->num_faces; ++f)
  {
    if (fused_mesh->face_cells[2*f+1] < -1)
    {
      // Found a ghost cell. Get its owning process.
      int proc = -fused_mesh->face_cells[2*f+1] + 2;
      ASSERT(proc != rank);
      if (!int_ptr_unordered_map_contains(send_map, proc))
        int_ptr_unordered_map_insert_with_v_dtor(send_map, proc, int_array_new(), DTOR(int_array_free));
      int_array_t* send_indices = *int_ptr_unordered_map_get(send_map, proc);
      int_array_append(send_indices, fused_mesh->face_cells[2*f]);
      if (!int_ptr_unordered_map_contains(recv_map, proc))
        int_ptr_unordered_map_insert_with_v_dtor(recv_map, proc, int_array_new(), DTOR(int_array_free));
      int_array_t* recv_indices = *int_ptr_unordered_map_get(recv_map, proc);
      int_array_append(recv_indices, ghost_cell);
      ++ghost_cell;
    }
  }
  ASSERT(ghost_cell == (fused_mesh->num_cells + fused_mesh->num_ghost_cells));
  exchanger_t* fused_ex = mesh_exchanger(fused_mesh);
  exchanger_set_sends(fused_ex, send_map);
  int_ptr_unordered_map_free(send_map);
  exchanger_set_receives(fused_ex, recv_map);
  int_ptr_unordered_map_free(recv_map);

  // Return the final fused mesh.
  return fused_mesh;
}

static void mesh_migrate(mesh_t** mesh, 
                         adj_graph_t* local_graph, 
                         SCOTCH_Num* local_partition,
                         exchanger_t* migrator)
{
  mesh_t* m = *mesh;
  index_t* vtx_dist = adj_graph_vertex_dist(local_graph);

  // Post receives for buffer sizes.
  int num_receives = exchanger_num_receives(migrator);
  int num_sends = exchanger_num_sends(migrator);
  int receive_buffer_sizes[num_receives], receive_procs[num_receives];
  int pos = 0, proc, num_indices, *indices, i_req = 0;
  MPI_Request requests[num_receives + num_sends];
  while (exchanger_next_receive(migrator, &pos, &proc, &indices, &num_indices))
  {
    receive_procs[i_req] = proc;
    MPI_Irecv(&receive_buffer_sizes[i_req], 1, MPI_INT, proc, 0, m->comm, &requests[i_req]);
    ++i_req;
  }

  // Build meshes to send to other processes.
  int_unordered_set_t* sent_cells = int_unordered_set_new();
  serializer_t* ser = mesh_serializer();
  byte_array_t* send_buffers[num_sends];
  int send_procs[num_sends];
  pos = 0;
  while (exchanger_next_send(migrator, &pos, &proc, &indices, &num_indices))
  {
    send_procs[i_req-num_receives] = proc;
    byte_array_t* bytes = byte_array_new();

    // Add the indices of the cells we are sending.
    for (int i = 0; i < num_indices; ++i)
      int_unordered_set_insert(sent_cells, indices[i]);

    // Create the mesh to send. Recall that the submesh encodes the destination
    // process rank in the ghost cells referenced in submesh->face_cells.
    mesh_t* submesh = create_submesh(m->comm, m, local_partition, vtx_dist, 
                                     indices, num_indices);

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
  mesh_t* submeshes[1+num_receives];
  for (int i = 0; i < num_receives; ++i)
  {
    size_t offset = 0;
    submeshes[i+1] = serializer_read(ser, receive_buffers[i], &offset);
  }

  // Clean up all the stuff from the exchange.
  ser = NULL;
  for (int i = 0; i < num_receives; ++i)
    byte_array_free(receive_buffers[i]);
  for (int i = 0; i < num_sends; ++i)
    byte_array_free(send_buffers[i]);

  // Construct a local submesh and store it in submeshes[0]. This submesh
  // consists of all cells not sent to other processes.
  {
    int num_cells = adj_graph_num_vertices(local_graph);
    int num_local_cells = num_cells - sent_cells->size;
    int local_cells[num_local_cells], j = 0;
    for (int i = 0; i < num_cells; ++i)
    {
      if (!int_unordered_set_contains(sent_cells, i))
        local_cells[j++] = i;
    }
    submeshes[0] = create_submesh(m->comm, m, local_partition, vtx_dist, 
                                  local_cells, num_local_cells);
  }

  // Fuse all the submeshes into a single mesh.
  mesh_free(m);
  *mesh = fuse_submeshes(submeshes, 1+num_receives);
}

#endif

exchanger_t* partition_point_cloud(point_cloud_t** cloud, MPI_Comm comm, int* weights, real_t imbalance_tol)
{
  ASSERT(imbalance_tol > 0.0);
  ASSERT(imbalance_tol <= 1.0);
  point_cloud_t* cl = *cloud;

#if POLYMEC_HAVE_MPI
  _Static_assert(sizeof(SCOTCH_Num) == sizeof(index_t), "SCOTCH_Num must be 64-bit.");

  int nprocs, rank;
  MPI_Comm_size(cl->comm, &nprocs);
  MPI_Comm_rank(cl->comm, &rank);

  // On a single process, partitioning has no meaning.
  if (nprocs == 1)
    return exchanger_new(cl->comm);

  // Clouds on rank != 0 must be NULL.
  ASSERT((rank == 0) || (*cloud == NULL));

  // Generate a global adjacency graph for the point cloud.
  adj_graph_t* global_graph = (cl != NULL) ? graph_from_point_cloud(cl) : NULL;

  // Map the graph to the different domains, producing a local partition vector.
  SCOTCH_Num* global_partition = (rank == 0) ? partition_graph(global_graph, comm, weights, imbalance_tol) : NULL;

  // Distribute the point cloud.
  point_cloud_distribute(cloud, comm, global_graph, global_partition);

  // Set up an exchanger to distribute field data.
  int num_vertices = (cl != NULL) ? adj_graph_num_vertices(global_graph) : 0;
  exchanger_t* distributor = create_distributor(comm, global_partition, num_vertices);

  // Clean up.
  adj_graph_free(global_graph);
  polymec_free(global_partition);

  // Return the migrator.
  return (distributor == NULL) ? exchanger_new(cl->comm) : distributor;
#else
  return exchanger_new(cl->comm);
#endif
}

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

  // Map the graph to the different domains, producing a local partition vector
  // (with values included for ghost cells).
  SCOTCH_Num* local_partition = repartition_graph(local_graph, cl->num_ghost_points, 
                                                  weights, imbalance_tol, cloud_ex);

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

exchanger_t* partition_mesh(mesh_t** mesh, MPI_Comm comm, int* weights, real_t imbalance_tol)
{
  ASSERT(imbalance_tol > 0.0);
  ASSERT(imbalance_tol <= 1.0);

#if POLYMEC_HAVE_MPI
  ASSERT((*mesh == NULL) || ((*mesh)->comm == MPI_COMM_SELF));
  _Static_assert(sizeof(SCOTCH_Num) == sizeof(index_t), "SCOTCH_Num must be 64-bit.");

  int nprocs, rank;
  MPI_Comm_size(comm, &nprocs);
  MPI_Comm_rank(comm, &rank);

  // On a single process, partitioning has no meaning.
  if (nprocs == 1)
    return exchanger_new(comm);

  // If meshes on rank != 0 are not NULL, we delete them.
  mesh_t* m = *mesh;
  if ((rank != 0) && (m != NULL))
  {
    mesh_free(m);
    *mesh = m = NULL; 
  }

  // Generate a global adjacency graph for the mesh.
  adj_graph_t* global_graph = (m != NULL) ? graph_from_mesh_cells(m) : NULL;

  // Map the graph to the different domains, producing a local partition vector.
  SCOTCH_Num* global_partition = (rank == 0) ? partition_graph(global_graph, comm, weights, imbalance_tol): NULL;

  // Distribute the mesh.
  mesh_distribute(mesh, comm, global_graph, global_partition);

  // Set up an exchanger to distribute field data.
  int num_vertices = (m != NULL) ? adj_graph_num_vertices(global_graph) : 0;
  exchanger_t* distributor = create_distributor(comm, global_partition, num_vertices);

  // Clean up.
  if (global_graph != NULL)
    adj_graph_free(global_graph);
  if (global_partition != NULL)
    polymec_free(global_partition);

  // Return the migrator.
  return distributor;
#else
  return exchanger_new(comm);
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
  mesh_migrate(mesh, local_graph, local_partition, migrator);

  // Clean up.
  adj_graph_free(local_graph);
  polymec_free(local_partition);

  // Return the migrator.
  return migrator;
#else
  return exchanger_new(m->comm);
#endif
}

