// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/partition_mesh.h"

#if POLYMEC_HAVE_MPI
#include "core/array_utils.h"
#include "core/unordered_set.h"
#include "core/kd_tree.h"
#include "core/hilbert.h"
#include "core/timer.h"
#include "ptscotch.h"

// This helper partitions a (serial) global graph, creating and returning a 
// global partition vector. It should only be called on rank 0. This is not 
// declared static because it may be used by other parts of polymec, but is 
// not part of the public API.
int64_t* partition_graph(adj_graph_t* global_graph, 
                         MPI_Comm comm,
                         int* weights,
                         real_t imbalance_tol)
{
  START_FUNCTION_TIMER();
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
  int64_t* global_partition = NULL;
  if (rank == 0)
  {
    global_partition = polymec_malloc(sizeof(int64_t) * num_global_vertices);
    SCOTCH_Strat strategy;
    SCOTCH_stratInit(&strategy);
    SCOTCH_Num strat_flags = SCOTCH_STRATDEFAULT;
    int result = SCOTCH_stratDgraphMapBuild(&strategy, strat_flags, nprocs, nprocs, (double)imbalance_tol);
    if (result != 0)
      polymec_error("Partitioning strategy could not be constructed.");
    result = SCOTCH_dgraphPart(&dist_graph, nprocs, &strategy, global_partition);
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
  STOP_FUNCTION_TIMER();
  return global_partition;
}

// This helper repartitions a local graph, creating and returning a 
// local partition vector with destination ranks included for ghost cells.
static int64_t* repartition_graph(adj_graph_t* local_graph, 
                                  int num_ghost_vertices,
                                  int* weights,
                                  real_t imbalance_tol,
                                  exchanger_t* local_graph_ex)
{
  START_FUNCTION_TIMER();
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
    index_t* global_indices = polymec_malloc(sizeof(index_t) * (num_vertices + num_ghost_vertices));
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
  int64_t* local_partition = polymec_malloc(sizeof(int64_t) * (num_vertices + num_ghost_vertices));
  {
    SCOTCH_Strat strategy;
    SCOTCH_stratInit(&strategy);
    SCOTCH_Num strat_flags = SCOTCH_STRATBALANCE;
    SCOTCH_stratDgraphMapBuild(&strategy, strat_flags, nprocs, nprocs, (double)imbalance_tol);
    int result = SCOTCH_dgraphPart(&dist_graph, nprocs, &strategy, local_partition);
    if (result != 0)
      polymec_error("Repartitioning failed.");
    SCOTCH_dgraphExit(&dist_graph);
    SCOTCH_stratExit(&strategy);
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

  // Return the local partition vector.
  STOP_FUNCTION_TIMER();
  return local_partition;
}

// This creates tags for a submesh whose elements belong to the given set.
static void create_submesh_tags(tagger_t* tagger, 
                                int_int_unordered_map_t* element_map)
{
  int pos = 0, *indices, size;
  char* tag_name;
  int_array_t* sub_indices = int_array_new();
  while (tagger_next_tag(tagger, &pos, &tag_name, &indices, &size))
  {
    int* tag = tagger_create_tag(tagger, tag_name, sub_indices->size);
    if (tag != NULL) // Skip tags that already exist!
    {
      for (int i = 0; i < size; ++i)
      {
        int* elem_p = int_int_unordered_map_get(element_map, indices[i]);
        if (elem_p != NULL)
          int_array_append(sub_indices, *elem_p);
      }
      memcpy(tag, sub_indices->data, sizeof(int) * sub_indices->size);
    }
  }
  int_array_free(sub_indices);
}

// This helper constructs and returns a mesh from the cells with the given
// indices in the given mesh. The submesh is valid with the following 
// exception:
// The indices of ghost cells referenced in submesh->face_cells are 
// replaced with destination process ranks.
static mesh_t* create_submesh(MPI_Comm comm, mesh_t* mesh, 
                              int64_t* partition, index_t* vtx_dist, 
                              int* indices, int num_indices)
{
  START_FUNCTION_TIMER();
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
      if ((opp_cell != -1) && !int_unordered_set_contains(cell_set, opp_cell))
      {
        int_unordered_set_insert(cell_set, opp_cell);
        ++num_ghost_cells;
      }

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

  // Create a mapping of parallel boundary faces to global "ghost" cells.
  // We will use this to create an annotated "parallel_boundary_faces" tag 
  // in the submesh.
  int_int_unordered_map_t* parallel_bface_map = int_int_unordered_map_new();

  // Copy face cells and face nodes.
  for (int f = 0; f < submesh->num_faces; ++f)
  {
    // Identify the cells attached to the face and construct them within the submesh.
    int orig_mesh_face = face_map[f];
    int* cell_p = int_int_unordered_map_get(inverse_cell_map, mesh->face_cells[2*orig_mesh_face]);
    int* opp_cell_p = int_int_unordered_map_get(inverse_cell_map, mesh->face_cells[2*orig_mesh_face+1]);
    ASSERT((cell_p != NULL) || (opp_cell_p != NULL)); // There must be at least one internal cell attached to this face.
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

      // We encode the destination process rank in the ghost cell and stash 
      // the original ghost index in our parallel boundary face map.
      if (ghost_cell != -1)
      {
        that_cell = -partition[ghost_cell] - 2;
        int_int_unordered_map_insert(parallel_bface_map, f, ghost_cell);
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

  // Create a tag for boundary faces and an associated property for 
  // global ghost cell indices. 
  {
    int* pbf_tag = mesh_create_tag(submesh->face_tags, "parallel_boundary_faces", 
                                   parallel_bface_map->size);
    int_array_t* pbf_ghost_cells = int_array_new();
    int pos = 0, face, gcell, i = 0;
    while (int_int_unordered_map_next(parallel_bface_map, &pos, &face, &gcell))
    {
      pbf_tag[i] = face;
      int_array_append(pbf_ghost_cells, gcell);
      ++i;
    }
    serializer_t* s = int_array_serializer();
    mesh_tag_set_property(submesh->face_tags, "parallel_boundary_faces", 
                          "ghost_cell_indices", pbf_ghost_cells, s);

    // Also set up a copy of the  array of global cell indices.
    int_array_t* global_cell_indices = int_array_new();
    int_array_resize(global_cell_indices, num_indices);
    memcpy(global_cell_indices->data, indices, sizeof(int) * num_indices);
    mesh_set_property(submesh, "global_cell_indices", global_cell_indices, s);
  }
  int_int_unordered_map_free(parallel_bface_map);

  // Copy node positions.
  for (int n = 0; n < submesh->num_nodes; ++n)
  {
    int orig_mesh_node = node_map[n];
    submesh->nodes[n] = mesh->nodes[orig_mesh_node];
  }

  // Construct edges.
  mesh_construct_edges(submesh);

  // Verify the submesh's topological correctness.
  ASSERT(mesh_verify_topology(submesh, polymec_error));

  // Do geometry.
  mesh_compute_geometry(submesh);

  // Dump any tags in the existing global mesh in, too.
  // FIXME: Edge tags are not supported!
  create_submesh_tags(submesh->cell_tags, inverse_cell_map);
  create_submesh_tags(submesh->face_tags, inverse_face_map);
  create_submesh_tags(submesh->node_tags, inverse_node_map);

  // Clean up.
  int_int_unordered_map_free(inverse_cell_map);
  int_int_unordered_map_free(inverse_face_map);
  int_int_unordered_map_free(inverse_node_map);
  polymec_free(face_map);
  polymec_free(node_map);

  // Copy properties using their serializers.
  int pos = 0;
  char* prop_name;
  void* prop_data;
  serializer_t* prop_ser;
  while (mesh_next_property(mesh, &pos, &prop_name, &prop_data, &prop_ser))
  {
    if (prop_ser != NULL)
    {
      void* prop_data_copy = serializer_clone_object(prop_ser, prop_data);
      mesh_set_property(submesh, prop_name, prop_data_copy, prop_ser);
    }
    else
      log_debug("create_submesh: property '%s' has no serializer.", prop_name);
  }

  STOP_FUNCTION_TIMER();
  return submesh;
}

static void sort_global_cell_pairs(int* indices, int num_pairs)
{
  START_FUNCTION_TIMER();
  int data[3*num_pairs]; // (i, j, swapped) for each pair
  for (int i = 0; i < num_pairs; ++i)
  {
    if (indices[2*i] > indices[2*i+1])
    {
      data[3*i]   = indices[2*i+1];
      data[3*i+1] = indices[2*i];
      data[3*i+2] = 1;
    }
    else
    {
      data[3*i]   = indices[2*i];
      data[3*i+1] = indices[2*i+1];
      data[3*i+2] = 0;
    }
  }
  // We can sort using the ordinary integer pair comparator, since the 
  // swapped flag doesn't affect the ordering.
  qsort(data, (size_t)num_pairs, 3*sizeof(int), int_pair_bsearch_comp);
  for (int i = 0; i < num_pairs; ++i)
  {
    if (data[3*i+2] == 1) // swapped
    {
      indices[2*i]   = data[3*i+1];
      indices[2*i+1] = data[3*i];
    }
    else
    {
      indices[2*i]   = data[3*i];
      indices[2*i+1] = data[3*i+1];
    }
  }
  STOP_FUNCTION_TIMER();
}

static void mesh_distribute(mesh_t** mesh, 
                            MPI_Comm comm,
                            adj_graph_t* global_graph, 
                            int64_t* global_partition)
{
  START_FUNCTION_TIMER();
  int nprocs, rank;
  MPI_Comm_size(comm, &nprocs);
  MPI_Comm_rank(comm, &rank);

  // Make sure we're all here.
  MPI_Barrier(comm);

  // We use the int_array serializer, so we need to make sure it's 
  // registered on all processes. See serializer.h for details.
  {
    serializer_t* s = int_array_serializer();
    s = NULL;
  }

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

  // Extract the boundary faces and the original (global) ghost cell indices
  // associated with them.
  ASSERT(mesh_has_tag(local_mesh->face_tags, "parallel_boundary_faces"));
  int num_pbfaces;
  int* pbfaces = mesh_tag(local_mesh->face_tags, "parallel_boundary_faces", &num_pbfaces);
  int_array_t* pbgcells = mesh_tag_property(local_mesh->face_tags, "parallel_boundary_faces", 
                                            "ghost_cell_indices");
  ASSERT(pbgcells != NULL);
  int_array_t* global_cell_indices = mesh_property(local_mesh, "global_cell_indices");
  ASSERT(global_cell_indices != NULL);

  // Here's an inverse map for global cell indices to the local ones we have.
  int_int_unordered_map_t* inverse_cell_map = int_int_unordered_map_new();
  for (int i = 0; i < local_mesh->num_cells; ++i)
    int_int_unordered_map_insert(inverse_cell_map, global_cell_indices->data[i], i);

  // Now we create pairwise global cell indices for the exchanger.
  int_ptr_unordered_map_t* ghost_cell_indices = int_ptr_unordered_map_new();
  int next_ghost_index = local_mesh->num_cells;
  for (int i = 0; i < num_pbfaces; ++i)
  {
    int f = pbfaces[i];
    ASSERT(local_mesh->face_cells[2*f+1] < -1);
    // Get the destination process for this ghost cell.
    int proc = -local_mesh->face_cells[2*f+1] - 2;
    ASSERT(proc >= 0);
    ASSERT(proc < nprocs);

    // Generate a mapping from a global index to its ghost cell.
    int global_ghost_index = pbgcells->data[i];
    int local_ghost_index;
    int* ghost_index_p = int_int_unordered_map_get(inverse_cell_map, global_ghost_index);
    if (ghost_index_p == NULL)
    {
      local_ghost_index = next_ghost_index++;
      int_int_unordered_map_insert(inverse_cell_map, global_ghost_index, local_ghost_index);
    }
    else
      local_ghost_index = *ghost_index_p;
    local_mesh->face_cells[2*f+1] = local_ghost_index;

    if (!int_ptr_unordered_map_contains(ghost_cell_indices, proc))
      int_ptr_unordered_map_insert_with_v_dtor(ghost_cell_indices, proc, int_array_new(), DTOR(int_array_free));
    int_array_t* indices = *int_ptr_unordered_map_get(ghost_cell_indices, proc);
    int local_cell = local_mesh->face_cells[2*f];
    int global_cell = global_cell_indices->data[local_cell];
    int_array_append(indices, global_cell);
    int_array_append(indices, global_ghost_index);
  }

  // Make sure everything lines up.
  ASSERT(next_ghost_index - local_mesh->num_cells == local_mesh->num_ghost_cells);

  exchanger_t* ex = mesh_exchanger(local_mesh);
  int pos = 0, proc;
  int_array_t* indices;
  while (int_ptr_unordered_map_next(ghost_cell_indices, &pos, &proc, (void**)&indices))
  {
    if (proc != rank)
    {
      ASSERT((indices->size > 0) && ((indices->size % 2) == 0));
      int num_pairs = indices->size/2;

      // Sort the indices array lexicographically by pairs so that all of the 
      // exchanger send/receive transactions have the same order across 
      // processes. This requires a specialized sort, since we have to 
      // arrange the integers within the pairs in ascending order, sort 
      // them, and then switch them back.
      sort_global_cell_pairs(indices->data, num_pairs);
      
      int send_indices[num_pairs], recv_indices[num_pairs];
      for (int i = 0; i < num_pairs; ++i)
      {
        send_indices[i] = *int_int_unordered_map_get(inverse_cell_map, indices->data[2*i]);
        ASSERT(send_indices[i] < local_mesh->num_cells);
        recv_indices[i] = *int_int_unordered_map_get(inverse_cell_map, indices->data[2*i+1]);
        ASSERT(recv_indices[i] >= local_mesh->num_cells);
      }
      exchanger_set_send(ex, proc, send_indices, num_pairs, true);
      exchanger_set_receive(ex, proc, recv_indices, num_pairs, true);
    }
  }

  // Clean up again.
  int_ptr_unordered_map_free(ghost_cell_indices);
  int_int_unordered_map_free(inverse_cell_map);

  // Destroy the tag and global cell index properties. We don't want to have 
  // to rely on them further.
  mesh_delete_tag(local_mesh->face_tags, "parallel_boundary_faces");
  mesh_delete_property(local_mesh, "global_cell_indices");
  STOP_FUNCTION_TIMER();
}

// This helper creates an array containing an index map that removes "holes" 
// left by mapped duplicates as mapped in dup_map. The range [0, num_indices) 
// is mapped. This is used to remap faces and nodes in the fuse_submeshes() 
// helper below. A newly-allocated array of length num_indices is returned.
static int* create_index_map_with_dups_removed(int num_indices, 
                                               int_int_unordered_map_t* dup_map)
{
  START_FUNCTION_TIMER();
  int* map = polymec_malloc(sizeof(int) * num_indices);
  int j = 0; // Mapped index.
  for (int i = 0; i < num_indices; ++i)
  {
    int* k_ptr = int_int_unordered_map_get(dup_map, i);
    if (k_ptr == NULL)
      map[i] = j++;
    else
      map[i] = *k_ptr;
  }
  STOP_FUNCTION_TIMER();
  return map;
}

// This helper takes an array of submeshes and stitches them all together into 
// one single mesh (contiguous or not) on the current domain. The submeshes are 
// consumed in the process.
static mesh_t* fuse_submeshes(mesh_t** submeshes, 
                              int num_submeshes)
{
  START_FUNCTION_TIMER();
  ASSERT(num_submeshes > 0);
  int rank, nprocs;
  MPI_Comm_rank(submeshes[0]->comm, &rank);
  MPI_Comm_size(submeshes[0]->comm, &nprocs);

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

        // These nodes need to be merged with their neighbors.
        int pos = 0, node;
        while (mesh_face_next_node(submesh, f, &pos, &node))
          int_unordered_set_insert(seam_nodes, submesh_node_offsets[m] + node);
      }
      else if (submesh->face_cells[2*f+1] <= -1) // belongs to another domain
      {
        // These nodes may need to be merged with their neighbors.
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
  if (seam_faces->size > 0)
  {
    // Translate the seam faces to an array.
    int seam_face_array[seam_faces->size];
    {
      int pos = 0, i = 0, j;
      while (int_unordered_set_next(seam_faces, &pos, &j))
        seam_face_array[i++] = j;
    }

    // Make a kd-tree of face centers.
    point_t* xf = polymec_malloc(sizeof(point_t) * seam_faces->size);
    for (int i = 0; i < seam_faces->size; ++i)
    {
      // Find the submesh that contains this face.
      int m = 0, j = seam_face_array[i];
      while (j >= submesh_face_offsets[m+1]) ++m;

      // Put the face center into this array.
      int f = j - submesh_face_offsets[m];
      ASSERT((f >= 0) && (f < submeshes[m]->num_faces));
      xf[i] = submeshes[m]->face_centers[f];
    }
    kd_tree_t* face_tree = kd_tree_new(xf, seam_faces->size);
    polymec_free(xf);

    // Now examine each seam face and find the 2 nearest faces to it.
    // One will be the face itself, and the other will be a duplicate. 
    // Map the one with the higher index to the lower index. Note that 
    // the indices of the points in the kd-tree should correspond to those
    // in the seam face array.
    for (int i = 0; i < seam_faces->size; ++i)
    {
      int nearest[2], j = seam_face_array[i];

      // Get the face center we're considering (in the same way we did above).
      int m = 0;
      while (j >= submesh_face_offsets[m+1]) ++m;
      int f = j - submesh_face_offsets[m];
      point_t* x = &submeshes[m]->face_centers[f];

      // Find the 2 nearest faces.
      kd_tree_nearest_n(face_tree, x, 2, nearest);

      // Merge the faces by mapping the one with the higher index to the 
      // lower index. 
      int index0 = seam_face_array[nearest[0]];
      int index1 = seam_face_array[nearest[1]];
      ASSERT(index0 != index1);
      ASSERT((index0 == j) || (index1 == j));
      int_int_unordered_map_insert(dup_face_map, MAX(index0, index1), MIN(index0, index1));

      // ALSO: Alter the submesh so that its interior ghost cells actually 
      // refer to the correct (flattened) interior cells. This may seem strange,
      // but it helps us to get the face->cell connectivity right later on.
      int j1 = (j == index0) ? index1 : index0;
      int m1 = 0;
      while (j1 >= submesh_face_offsets[m1+1]) ++m1;
      int f1 = j1 - submesh_face_offsets[m1];
      ASSERT((f1 >= 0) && (f1 < submeshes[m1]->num_faces));
      int flattened_cell1 = submesh_cell_offsets[m1] + submeshes[m1]->face_cells[2*f1];
      submeshes[m]->face_cells[2*f+1] = flattened_cell1;
//printf("%d: Munging submesh[%d] face %d to attach to cell %d\n", rank, m, f, flattened_cell1);
    }

    // Clean up.
    kd_tree_free(face_tree);
  }
  int* face_map = create_index_map_with_dups_removed(num_faces, dup_face_map); 

  // Do the same for duplicate nodes.
  int_int_unordered_map_t* dup_node_map = int_int_unordered_map_new();
  if (seam_nodes->size > 0)
  {
    // Translate the seam nodes to an array.
    int seam_node_array[seam_nodes->size];
    {
      int pos = 0, i = 0, j;
      while (int_unordered_set_next(seam_nodes, &pos, &j))
        seam_node_array[i++] = j;
    }

    // Make a kd-tree of node positions.
    point_t* xn = polymec_malloc(sizeof(point_t) * seam_nodes->size);
    for (int i = 0; i < seam_nodes->size; ++i)
    {
      // Find the submesh that contains this face.
      int m = 0, j = seam_node_array[i];
      while (j >= submesh_node_offsets[m+1]) ++m;

      // Put the node position into this array.
      int n = j - submesh_node_offsets[m];
      xn[i] = submeshes[m]->nodes[n];
    }
    kd_tree_t* node_tree = kd_tree_new(xn, seam_nodes->size);
    polymec_free(xn);

    // Now examine each seam node and find ALL NODES that appear to be 
    // the same. NOTE: for now, we find all the nodes within a distance 
    // of epsilon, which is KNOWN NOT TO BE A ROBUST METHOD OF FUSING NODES.
    // It might be good enough for now, though...
    for (int i = 0; i < seam_nodes->size; ++i)
    {
      int j = seam_node_array[i];

      // Get the node position we're considering (in the same way we did above).
      int m = 0;
      while (j >= submesh_node_offsets[m+1]) ++m;
      int n = j - submesh_node_offsets[m];
      point_t* x = &submeshes[m]->nodes[n];

      // Find all the nodes within epsilon of this one.
      // Note that not all "seam nodes" will be merged with others, since we 
      // identify seam node candidates from a broader class of nodes including
      // nodes on the problem boundary that don't necessary belong to a seam 
      // face.
      static const real_t epsilon = 1e-12;
      int_array_t* same_nodes = kd_tree_within_radius(node_tree, x, epsilon);

      if (same_nodes->size > 1)
      {
        // Merge the nodes by mapping all those with higher indices to the 
        // lowest index.
        int min_index = INT_MAX;
        for (int k = 0; k < same_nodes->size; ++k)
        {
          int index = same_nodes->data[k];
          if (index < min_index)
            min_index = index;
        }
        for (int k = 0; k < same_nodes->size; ++k)
        {
          int index = same_nodes->data[k];
          if (index != min_index)
            int_int_unordered_map_insert(dup_node_map, index, min_index);
        }
      }
      int_array_free(same_nodes);
    }

    // Clean up.
    kd_tree_free(node_tree);
  }
  int* node_map = create_index_map_with_dups_removed(num_nodes, dup_node_map); 

  // Reduce the number of faces and nodes in the fused mesh by the ones that 
  // have been merged to others.
  num_faces -= dup_face_map->size;
  num_nodes -= dup_node_map->size;

  // We're through with the seam faces/nodes and the duplicate node maps.
  // We still need the duplicate face map below.
  int_unordered_set_free(seam_faces);
  int_unordered_set_free(seam_nodes);
  int_int_unordered_map_free(dup_node_map);

  // Now we create the fused mesh and fill it with the contents of the submeshes.
  mesh_t* fused_mesh = mesh_new(submeshes[0]->comm, num_cells, num_ghost_cells,
                                num_faces, num_nodes);

  // Allocate storage for cell faces and face nodes.
  fused_mesh->cell_face_offsets[0] = 0;
  int cell = 0;
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
      int face = face_map[flattened_face];
      int num_face_nodes = submesh->face_node_offsets[f+1] - submesh->face_node_offsets[f];
      fused_mesh->face_node_offsets[face+1] = fused_mesh->face_node_offsets[face] + num_face_nodes;
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
        int subface_index = submesh->cell_faces[submesh->cell_face_offsets[c] + f];
        bool flipped = false;
        if (subface_index < 0)
        {
          flipped = true;
          subface_index = ~subface_index;
        }
        int flattened_face = submesh_face_offsets[m] + subface_index;
        int face_index = face_map[flattened_face];
//        int* face_p = int_int_unordered_map_get(dup_face_map, flattened_face);
//        int face_index = (face_p != NULL) ? face_map[*face_p] : face_map[flattened_face];
        ASSERT(face_index < fused_mesh->num_faces);
//printf("%d: Cell %d has face %d\n", rank, cell, face_index);
        if (flipped)
          face_index = ~face_index;
        fused_mesh->cell_faces[fused_mesh->cell_face_offsets[cell]+f] = face_index;
      }
    }
  }

  // Copy face cells and face nodes.
  for (int m = 0; m < num_submeshes; ++m)
  {
    mesh_t* submesh = submeshes[m];
    for (int f = 0; f < submesh->num_faces; ++f)
    {
      int flattened_face = submesh_face_offsets[m] + f;
      int face = face_map[flattened_face];

      // Set up the cells of the face. 
      int subcell_index = submesh->face_cells[2*f];
      int flattened_cell = submesh_cell_offsets[m] + subcell_index;

      // To get the second cell attached to this face, we ask whether this 
      // face is on a seam. If it is, we can read off the second flattened 
      // cell index because we stored it directly in the submesh.
      int flattened_cell1;
      if (int_int_unordered_map_contains(dup_face_map, flattened_face))
        flattened_cell1 = submesh->face_cells[2*f+1];
      else
      {
        if (submesh->face_cells[2*f+1] < -1)
        {
          // Otherwise, if it's a ghost cell (< -1), it has its owning 
          // process encoded. If so, we simply copy it into place and 
          // deal with it later.
          flattened_cell1 = submesh->face_cells[2*f+1];
        }
        else if (submesh->face_cells[2*f+1] == -1)
        {
          // If it's on a boundary, the other cell is -1.
          flattened_cell1 = -1;
        }
        else 
        {
          // Otherwise, we construct its flattened index from this 
          // same submesh.
          ASSERT(submesh->face_cells[2*f+1] >= 0);
          flattened_cell1 = submesh_cell_offsets[m] + submesh->face_cells[2*f+1];
        }
      }
      fused_mesh->face_cells[2*face] = flattened_cell;
      fused_mesh->face_cells[2*face+1] = flattened_cell1;
//printf("%d: Face %d (from [%d, %d]) has cells %d, %d\n", rank, face, m, f, flattened_cell, flattened_cell1);

      // Set up the nodes of the face.
      int num_face_nodes = submesh->face_node_offsets[f+1] - submesh->face_node_offsets[f];
      for (int n = 0; n < num_face_nodes; ++n)
      {
        int subnode_index = submesh->face_nodes[submesh->face_node_offsets[f] + n];
        int flattened_node = submesh_node_offsets[m] + subnode_index;
        int node_index = node_map[flattened_node];
//        int* node_p = int_int_unordered_map_get(dup_node_map, flattened_node);
//        int node_index = (node_p != NULL) ? node_map[*node_p] : node_map[flattened_node];
printf("node index: %d of %d\n", node_index, fused_mesh->num_nodes);
        ASSERT(node_index < fused_mesh->num_nodes);
        fused_mesh->face_nodes[fused_mesh->face_node_offsets[face]+n] = node_index;
      }
    }
  }
  polymec_free(face_map);
  int_int_unordered_map_free(dup_face_map);

  // Copy node positions.
  for (int m = 0; m < num_submeshes; ++m)
  {
    mesh_t* submesh = submeshes[m];
    for (int n = 0; n < submesh->num_nodes; ++n)
    {
      int flattened_node = submesh_node_offsets[m] + n;
      int node = node_map[flattened_node];
      fused_mesh->nodes[node] = submesh->nodes[n];
    }
  }
  polymec_free(node_map);

  // Construct edges.
  mesh_construct_edges(fused_mesh);

  // At this point, we can verify the topological correctness of the fused mesh.
  ASSERT(mesh_verify_topology(fused_mesh, polymec_error));

  // Now compute geometry.
  mesh_compute_geometry(fused_mesh);

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
      int proc = -fused_mesh->face_cells[2*f+1] - 2;
      ASSERT(proc != rank);
      ASSERT(proc >= 0);
      ASSERT(proc < nprocs);

      // Set up the sends and receives associated with this process.
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
  STOP_FUNCTION_TIMER();
  return fused_mesh;
}

static void mesh_migrate(mesh_t** mesh, 
                         adj_graph_t* local_graph, 
                         int64_t* local_partition,
                         exchanger_t* migrator)
{
  START_FUNCTION_TIMER();
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
  int_unordered_set_free(sent_cells);
  mesh_free(m);
  *mesh = fuse_submeshes(submeshes, 1+num_receives);
  STOP_FUNCTION_TIMER();
}

#endif

exchanger_t* partition_mesh(mesh_t** mesh, MPI_Comm comm, int* weights, real_t imbalance_tol)
{
#if POLYMEC_HAVE_MPI
  START_FUNCTION_TIMER();
  ASSERT((weights == NULL) || (imbalance_tol > 0.0));
  ASSERT((weights == NULL) || (imbalance_tol <= 1.0));
  ASSERT((*mesh == NULL) || ((*mesh)->comm == MPI_COMM_SELF));
  _Static_assert(sizeof(SCOTCH_Num) == sizeof(int64_t), "SCOTCH_Num must be 64-bit.");

  int nprocs, rank;
  MPI_Comm_size(comm, &nprocs);
  MPI_Comm_rank(comm, &rank);
  ASSERT((rank != 0) || (*mesh != NULL));

  // On a single process, partitioning has no meaning, but we do replace the communicator
  // if needed. NOTE: the exchanger will still have its original communicator, but this 
  // shouldn't matter in any practical sense.
  if (nprocs == 1)
  {
    if (comm != (*mesh)->comm)
      (*mesh)->comm = comm;
    STOP_FUNCTION_TIMER();
    return exchanger_new(comm);
  }

  log_debug("partition_mesh: Partitioning mesh into %d subdomains.", nprocs);

  // If meshes on rank != 0 are not NULL, we delete them.
  mesh_t* m = *mesh;
  if ((rank != 0) && (m != NULL))
  {
    mesh_free(m);
    *mesh = m = NULL; 
  }

  // Generate a global adjacency graph for the mesh.
  adj_graph_t* global_graph = (m != NULL) ? graph_from_mesh_cells(m) : NULL;

#ifndef NDEBUG
  // Make sure there are enough cells to go around for the processes we're given.
  if (rank == 0)
  {
    ASSERT((*mesh)->num_cells > nprocs);
  }
#endif

  // Map the graph to the different domains, producing a local partition vector.
  int64_t* global_partition = (rank == 0) ? partition_graph(global_graph, comm, weights, imbalance_tol): NULL;

  // Distribute the mesh.
  log_debug("partition_mesh: Distributing mesh to %d processes.", nprocs);
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
  STOP_FUNCTION_TIMER();
  return distributor;
#else
  // Replace the communicator if needed.
  if (comm != (*mesh)->comm)
    (*mesh)->comm = comm;
  return exchanger_new(comm);
#endif
}

int64_t* partition_vector_from_mesh(mesh_t* global_mesh, MPI_Comm comm, int* weights, real_t imbalance_tol)
{
#if POLYMEC_HAVE_MPI
  START_FUNCTION_TIMER();
  ASSERT((weights == NULL) || (imbalance_tol > 0.0));
  ASSERT((weights == NULL) || (imbalance_tol <= 1.0));
  _Static_assert(sizeof(SCOTCH_Num) == sizeof(int64_t), "SCOTCH_Num must be 64-bit.");

  int nprocs, rank;
  MPI_Comm_size(comm, &nprocs);
  MPI_Comm_rank(comm, &rank);
  ASSERT((rank != 0) || (global_mesh != NULL));

  // On a single process, partitioning has no meaning.
  if (nprocs == 1)
  {
    // Dumb, but correct.
    int64_t* global_partition = polymec_malloc(sizeof(int64_t) * global_mesh->num_cells);
    memset(global_partition, 0, sizeof(int64_t) * global_mesh->num_cells);
    return global_partition;
  }

  // Generate a global adjacency graph for the mesh.
  adj_graph_t* global_graph = (rank == 0) ? graph_from_mesh_cells(global_mesh) : NULL;

#ifndef NDEBUG
  // Make sure there are enough cells to go around for the processes we're given.
  if (rank == 0)
  {
    ASSERT(global_mesh->num_cells > nprocs);
  }
#endif

  // Map the graph to the different domains, producing a local partition vector.
  int64_t* global_partition = (rank == 0) ? partition_graph(global_graph, comm, weights, imbalance_tol): NULL;

  // Get rid of the graph.
  if (global_graph != NULL)
    adj_graph_free(global_graph);

  STOP_FUNCTION_TIMER();
  return global_partition;

#else
  // This is dumb, but we were asked for it.
  int64_t* global_partition = polymec_malloc(sizeof(int64_t) * global_mesh->num_cells);
  memset(global_partition, 0, sizeof(int64_t) * global_mesh->num_cells);
  return global_partition;
#endif
}

exchanger_t* distribute_mesh(mesh_t** mesh, MPI_Comm comm, int64_t* global_partition)
{
#if POLYMEC_HAVE_MPI
  ASSERT((*mesh == NULL) || ((*mesh)->comm == MPI_COMM_SELF));

  int nprocs, rank;
  MPI_Comm_size(comm, &nprocs);
  MPI_Comm_rank(comm, &rank);
  ASSERT((rank != 0) || (global_partition != NULL));
  ASSERT((rank != 0) || (*mesh != NULL));

  // On a single process, partitioning has no meaning.
  if (nprocs == 1)
    return exchanger_new(comm);

  START_FUNCTION_TIMER();

  // If meshes on rank != 0 are not NULL, we delete them.
  mesh_t* m = *mesh;
  if ((rank != 0) && (m != NULL))
  {
    mesh_free(m);
    *mesh = m = NULL; 
  }

  // Generate a global adjacency graph for the mesh.
  adj_graph_t* global_graph = (m != NULL) ? graph_from_mesh_cells(m) : NULL;

  // Distribute the mesh.
  mesh_distribute(mesh, comm, global_graph, global_partition);

  // Set up an exchanger to distribute field data.
  int num_vertices = (m != NULL) ? adj_graph_num_vertices(global_graph) : 0;
  exchanger_t* distributor = create_distributor(comm, global_partition, num_vertices);

  // Get rid of the graph.
  if (global_graph != NULL)
    adj_graph_free(global_graph);

  STOP_FUNCTION_TIMER();
  return distributor;
#else
  return exchanger_new(comm);
#endif
}

exchanger_t* repartition_mesh(mesh_t** mesh, int* weights, real_t imbalance_tol)
{
  POLYMEC_NOT_IMPLEMENTED;

  ASSERT(imbalance_tol > 0.0);
  ASSERT(imbalance_tol <= 1.0);
  mesh_t* m = *mesh;

#if POLYMEC_HAVE_MPI
  START_FUNCTION_TIMER();
  _Static_assert(sizeof(SCOTCH_Num) == sizeof(int64_t), "SCOTCH_Num must be 64-bit.");

  int nprocs, rank;
  MPI_Comm_size(m->comm, &nprocs);
  MPI_Comm_rank(m->comm, &rank);

  // On a single process, repartitioning has no meaning.
  if (nprocs == 1)
    return exchanger_new(m->comm);

  // Generate a local adjacency graph for the mesh.
  adj_graph_t* local_graph = graph_from_mesh_cells(m);

#ifndef NDEBUG
  // Make sure there are enough cells to go around for the processes we're given.
  index_t* vtx_dist = adj_graph_vertex_dist(local_graph);
  index_t total_num_cells = vtx_dist[nprocs];
  ASSERT(total_num_cells >= nprocs);
#endif

  // Get the exchanger for the mesh.
  exchanger_t* mesh_ex = mesh_exchanger(m);

  // Map the graph to the different domains, producing a local partition vector.
  int64_t* local_partition = repartition_graph(local_graph, m->num_ghost_cells, weights, imbalance_tol, mesh_ex);

  // Set up an exchanger to migrate field data.
  int num_vertices = adj_graph_num_vertices(local_graph);
  exchanger_t* migrator = create_migrator(m->comm, local_partition, num_vertices);

  // Migrate the mesh.
  mesh_migrate(mesh, local_graph, local_partition, migrator);

  // Clean up.
  adj_graph_free(local_graph);
  polymec_free(local_partition);

  // Return the migrator.
  STOP_FUNCTION_TIMER();
  return migrator;
#else
  return exchanger_new(m->comm);
#endif
}

