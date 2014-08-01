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

// This helper partitions a (serial) global graph, creating and returning a 
// global partition vector on rank 0.
static SCOTCH_Num* partition_graph(adj_graph_t* global_graph, 
                                   int* weights,
                                   real_t imbalance_tol)
{
  ASSERT(adj_graph_comm(global_graph) == MPI_COMM_WORLD);

  int nprocs, rank;
  MPI_Comm comm = adj_graph_comm(global_graph);
  MPI_Comm_size(comm, &nprocs);
  MPI_Comm_rank(comm, &rank);

  ASSERT((rank == 0) || (global_graph == NULL));

  SCOTCH_Dgraph dist_graph;
  SCOTCH_Num* vtx_weights = NULL;
  int num_global_vertices = 0;
  if (rank == 0)
  {
    num_global_vertices = adj_graph_num_vertices(global_graph);

    // Extract the adjacency information.
    SCOTCH_Num* xadj = malloc(sizeof(SCOTCH_Num) * (num_global_vertices+1));
    int* edge_offsets = adj_graph_edge_offsets(global_graph);
    for (int i = 0; i <= num_global_vertices; ++i)
      xadj[i] = (SCOTCH_Num)edge_offsets[i];
    SCOTCH_Num num_arcs = xadj[num_global_vertices];
    SCOTCH_Num* adj = malloc(sizeof(SCOTCH_Num) * num_arcs);
    int* edges = adj_graph_adjacency(global_graph);
    for (int i = 0; i < xadj[num_global_vertices]; ++i)
      adj[i] = (SCOTCH_Num)edges[i];

    // Build a graph on rank 0.
    SCOTCH_dgraphInit(&dist_graph, comm);
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
    SCOTCH_stratExit(&strategy);

    if (vtx_weights != NULL)
      polymec_free(vtx_weights);
  }

  // Return the global partition vector.
  return global_partition;
}

// This helper repartitions a local graph, creating and returning a 
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

// This helper sets up the exchanger ex so that it can distribute data from the 
// root process (0) according to the global partition vector. The partition 
// vector is NULL on all nonzero ranks and defined on rank 0.
static exchanger_t* create_distributor(SCOTCH_Num* global_partition,
                                       int num_global_vertices)
{
  MPI_Comm comm = MPI_COMM_WORLD;

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

// This helper constructs and returns a point cloud from the given points in 
// the given point cloud.
static point_cloud_t* create_subcloud(point_cloud_t* cloud, int* indices, int num_indices)
{
  return NULL;
}

static void point_cloud_distribute(point_cloud_t** cloud, 
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
// indices in the given mesh.
static mesh_t* create_submesh(mesh_t* mesh, index_t* vtx_dist, int* indices, int num_indices)
{
  int rank;
  MPI_Comm_rank(mesh->comm, &rank);

  // Construct a set of global cell indices for ghosts on this mesh.
  uint64_t* global_cell_indices = polymec_malloc(sizeof(uint64_t) * (mesh->num_cells + mesh->num_ghost_cells));
  for (int i = 0; i < mesh->num_cells; ++i)
    global_cell_indices[i] = vtx_dist[rank] + i;
  exchanger_t* ex = mesh_exchanger(mesh);
  exchanger_exchange(ex, global_cell_indices, 1, 0, MPI_UINT64_T);

  // Count unique mesh elements.
  int num_cells = num_indices, num_ghost_cells = 0;
  int_unordered_set_t* face_indices = int_unordered_set_new();
  int_unordered_set_t* node_indices = int_unordered_set_new();
  for (int i = 0; i < num_indices; ++i)
  {
    int cell = indices[i], fpos = 0, face;
    while (mesh_cell_next_face(mesh, cell, &fpos, &face))
    {
      if (mesh_face_opp_cell(mesh, face, cell) >= mesh->num_cells)
        ++num_ghost_cells;
      int_unordered_set_insert(face_indices, face);
      int npos = 0, node;
      while (mesh_face_next_node(mesh, face, &npos, &node))
        int_unordered_set_insert(node_indices, node);
    }
  }

  // Create the submesh container.
  int num_faces = face_indices->size, num_nodes = node_indices->size;
  mesh_t* submesh = mesh_new(mesh->comm, num_cells, num_ghost_cells, num_faces, num_nodes);

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
  for (int c = 0; c < num_cells; ++c)
  {
    int_int_unordered_map_insert(inverse_cell_map, indices[c], c);
    int num_cell_faces = mesh->cell_face_offsets[indices[c]+1] - mesh->cell_face_offsets[indices[c]];
    submesh->cell_face_offsets[c+1] = submesh->cell_face_offsets[c] + num_cell_faces;
  }
  int_int_unordered_map_t* inverse_face_map = int_int_unordered_map_new();
  submesh->face_node_offsets[0] = 0;
  for (int f = 0; f < num_faces; ++f)
  {
    int_int_unordered_map_insert(inverse_face_map, face_map[f], f);
    int num_face_nodes = mesh->face_node_offsets[face_map[f]+1] - mesh->face_node_offsets[face_map[f]];
    submesh->face_node_offsets[f+1] = submesh->face_node_offsets[f] + num_face_nodes;
  }
  int_int_unordered_map_t* inverse_node_map = int_int_unordered_map_new();
  for (int n = 0; n < num_nodes; ++n)
    int_int_unordered_map_insert(inverse_node_map, node_map[n], n);
  mesh_reserve_connectivity_storage(submesh);

  // Copy cell faces.
  for (int c = 0; c < num_cells; ++c)
  {
    int num_cell_faces = submesh->cell_face_offsets[c+1] - submesh->cell_face_offsets[c];
    for (int f = 0; f < num_cell_faces; ++f)
      submesh->cell_faces[submesh->cell_face_offsets[c]+f] = *int_int_unordered_map_get(inverse_face_map, mesh->cell_faces[indices[c]+f]);
  }

  // Copy face cells and face nodes.
  for (int f = 0; f < num_faces; ++f)
  {
    submesh->face_cells[2*f] = *int_int_unordered_map_get(inverse_cell_map, mesh->face_cells[2*face_map[f]]);
    int* opp_cell_p = int_int_unordered_map_get(inverse_cell_map, mesh->face_cells[2*face_map[f]+1]);
    if (opp_cell_p != NULL)
      submesh->face_cells[2*f+1] = *opp_cell_p;
    else
    {
      // This is a ghost cell. Either it belongs to the mesh we're carving up, 
      // or we can find it in the exchanger.
      int ghost_cell = mesh->face_cells[2*face_map[f]+1];
      int global_index;
      if (ghost_cell < mesh->num_cells) 
      {
        // It belongs to the mesh. 
        global_index = ghost_cell + vtx_dist[rank];
      }
      else
      {
        // It's actually a ghost cell of the mesh. Fetch its global cell index 
        // from the exchanger.
        ASSERT(ghost_cell >= mesh->num_cells);
        global_index = global_cell_indices[ghost_cell];
      }

      // Encode the global cell index in this entry of the submesh.
      submesh->face_cells[2*f+1] = -global_index - 2;
    }
    int num_face_nodes = submesh->face_node_offsets[f+1] - submesh->face_node_offsets[f];
    for (int n = 0; n < num_face_nodes; ++n)
      submesh->face_nodes[submesh->face_node_offsets[f]+n] = *int_int_unordered_map_get(inverse_node_map, mesh->face_nodes[face_map[f]+n]);
  }

  // Copy node positions.
  for (int n = 0; n < num_nodes; ++n)
  {
    submesh->nodes[n].x = mesh->nodes[node_map[n]].x;
    submesh->nodes[n].y = mesh->nodes[node_map[n]].y;
    submesh->nodes[n].z = mesh->nodes[node_map[n]].z;
  }

  // Finish things up.
  // FIXME: Is this needed?
  mesh_construct_edges(submesh);
  mesh_compute_geometry(submesh);

  // Clean up.
  polymec_free(global_cell_indices);
  int_int_unordered_map_free(inverse_cell_map);
  int_int_unordered_map_free(inverse_face_map);
  int_int_unordered_map_free(inverse_node_map);
  polymec_free(face_map);
  polymec_free(node_map);

  return submesh;
}

static void mesh_distribute(mesh_t** mesh, 
                            adj_graph_t* global_graph, 
                            SCOTCH_Num* global_partition)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  int nprocs, rank;
  MPI_Comm_size(comm, &nprocs);
  MPI_Comm_rank(comm, &rank);

  mesh_t* global_mesh = *mesh;
  mesh_t* local_mesh = NULL;
  index_t* vtx_dist = adj_graph_vertex_dist(global_graph);
  if (rank == 0)
  {
    // Take stock of how many cells we'll have per process.
    int num_cells[nprocs];
    memset(num_cells, 0, sizeof(int) * nprocs);
    for (int i = 0; i < global_mesh->num_cells; ++i)
      num_cells[global_partition[i]]++;

    // Carve out the portion of the mesh that will stick around on process 0.
    {
      int indices[num_cells[0]], k = 0;
      for (int i = 0; i < global_mesh->num_cells; ++i)
      {
        if (global_partition[i] == rank)
          indices[k++] = i;
      }
      local_mesh = create_submesh(global_mesh, vtx_dist, indices, num_cells[0]);
    }

    // Now do the other processes.
    serializer_t* ser = mesh_serializer();
    byte_array_t* bytes = byte_array_new();
    for (int p = 1; p < nprocs; ++p)
    {
      // Create the pth submesh.
      int indices[num_cells[p]], k = 0;
      for (int i = 0; i < global_mesh->num_cells; ++i)
      {
        if (global_partition[i] == p)
          indices[k++] = i;
      }
      mesh_t* p_mesh = create_submesh(global_mesh, vtx_dist, indices, num_cells[p]);

      // Serialize it and send its size (and it) to process p.
      size_t offset = 0;
      serializer_write(ser, p_mesh, bytes, &offset);
      MPI_Send(&bytes->size, 1, MPI_INT, p, p, comm);
      MPI_Send(&bytes->data, bytes->size, MPI_BYTE, p, p, comm);

      // Clean up.
      byte_array_clear(bytes);
      mesh_free(p_mesh);
    }
    ser = NULL;
    byte_array_free(bytes);
  }
  else
  {
    // Receive the size of the incoming mesh.
    int mesh_size;
    MPI_Status status;
    MPI_Recv(&mesh_size, 1, MPI_INT, 0, rank, comm, &status);

    // Now receive the mesh.
    byte_array_t* bytes = byte_array_new();
    byte_array_resize(bytes, mesh_size);
    MPI_Recv(bytes->data, mesh_size, MPI_BYTE, 0, rank, comm, &status);
    size_t offset = 0;
    serializer_t* ser = mesh_serializer();
    local_mesh = serializer_read(ser, bytes, &offset);
    
    byte_array_free(bytes);
    ser = NULL;
  }

  *mesh = local_mesh;

  // Clean up.
  mesh_free(global_mesh);

  // Now we create the exchanger using the encoded ghost cell indices.
  int_ptr_unordered_map_t* ghost_cell_indices = int_ptr_unordered_map_new();
  for (int f = 0; f < local_mesh->num_faces; ++f)
  {
    if (local_mesh->face_cells[2*f+1] < -1) 
    {
      int global_ghost_index = -(local_mesh->face_cells[2*f+1] + 2);

      // Find the process and local index for this ghost cell.
      int proc = 0;
      while (vtx_dist[proc+1] < global_ghost_index) ++proc;
      int local_ghost_index = global_ghost_index - vtx_dist[proc];

      if (!int_ptr_unordered_map_contains(ghost_cell_indices, proc))
        int_ptr_unordered_map_insert_with_v_dtor(ghost_cell_indices, proc, int_array_new(), DTOR(int_array_free));
      int_array_t* indices = *int_ptr_unordered_map_get(ghost_cell_indices, proc);
      int_array_append(indices, local_mesh->face_cells[2*f]);
      int_array_append(indices, local_ghost_index);
    }
  }
  exchanger_t* ex = mesh_exchanger(local_mesh);
  int pos = 0, proc;
  index_array_t* indices;
  while (int_ptr_unordered_map_next(ghost_cell_indices, &pos, &proc, (void**)&indices))
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

  // Clean up again.
  int_ptr_unordered_map_free(ghost_cell_indices);
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

exchanger_t* partition_point_cloud(point_cloud_t** cloud, int* weights, real_t imbalance_tol)
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
  SCOTCH_Num* global_partition = partition_graph(global_graph, weights, imbalance_tol);

  // Distribute the point cloud.
  point_cloud_distribute(cloud, global_graph, global_partition);

  // Set up an exchanger to distribute field data.
  int num_vertices = (cl != NULL) ? adj_graph_num_vertices(global_graph) : 0;
  exchanger_t* distributor = create_distributor(global_partition, num_vertices);

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

exchanger_t* partition_mesh(mesh_t** mesh, int* weights, real_t imbalance_tol)
{
  ASSERT(imbalance_tol > 0.0);
  ASSERT(imbalance_tol <= 1.0);

  MPI_Comm comm = MPI_COMM_WORLD;

#if POLYMEC_HAVE_MPI
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
    m = NULL; 
  }

  // Generate a global adjacency graph for the mesh.
  adj_graph_t* global_graph = (m != NULL) ? graph_from_mesh_cells(m) : NULL;

  // Map the graph to the different domains, producing a local partition vector.
  SCOTCH_Num* global_partition = partition_graph(global_graph, weights, imbalance_tol);

  // Distribute the mesh.
  mesh_distribute(mesh, global_graph, global_partition);

  // Set up an exchanger to distribute field data.
  int num_vertices = (m != NULL) ? adj_graph_num_vertices(global_graph) : 0;
  exchanger_t* distributor = create_distributor(global_partition, num_vertices);

  // Clean up.
  adj_graph_free(global_graph);
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

