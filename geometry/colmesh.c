// Copyright (c) 2012-2018, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/timer.h"
#include "core/array.h"
#include "core/array_utils.h"
#include "core/exchanger.h"
#include "core/unordered_map.h"
#include "core/partitioning.h"
#include "geometry/colmesh.h"
#include "geometry/colmesh_field.h"
#include "geometry/polymesh.h"

// Chunk xy data. Shared across all "stacked" chunks.
typedef struct 
{
  size_t num_columns;
  size_t num_ghost_columns;
  int* column_xy_face_offsets;
  int* column_xy_faces;
  int* xy_face_columns;
  size_t num_xy_faces;
  size_t num_xy_edges;
  int* xy_edge_nodes;
  size_t num_xy_nodes;
  point2_t* xy_nodes;

  // Exchanger process maps--used to construct exchangers.
  exchanger_proc_map_t* send_map;
  exchanger_proc_map_t* receive_map;
} chunk_xy_data_t;

DEFINE_ARRAY(chunk_xy_data_array, chunk_xy_data_t*)

// Creates shared xy chunk data for the local process, generating information
// for the locally-owned cells in the planar polymesh.
static chunk_xy_data_t* chunk_xy_data_new(MPI_Comm comm,
                                          planar_polymesh_t* mesh,
                                          int64_t* partition_vector,
                                          int xy_chunk_index)
{
  chunk_xy_data_t* xy_data = polymec_malloc(sizeof(chunk_xy_data_t));
  xy_data->num_columns = 0;
  xy_data->num_ghost_columns = 0;
  xy_data->num_xy_faces = 0;
  xy_data->num_xy_edges = 0;
  xy_data->num_xy_nodes = 0;
  xy_data->send_map = exchanger_proc_map_new();
  xy_data->receive_map = exchanger_proc_map_new();

  int rank;
  MPI_Comm_rank(comm, &rank);

  // Make a list of polygonal cells (which become columns in our chunks) that 
  // correspond to the given xy chunk index.
  int_array_t* local_cells = int_array_new(); // maps columns to planar mesh cells
  int_int_unordered_map_t* cell_to_col_map = 
    int_int_unordered_map_new(); // maps planar mesh cells to columns in our chunk
  for (int cell = 0; cell < mesh->num_cells; ++cell)
  {
    if ((int)partition_vector[cell] == xy_chunk_index)
    {
      int_array_append(local_cells, cell);
      int_int_unordered_map_insert(cell_to_col_map, cell, (int)xy_data->num_columns);
      ++xy_data->num_columns;
    }
  }

  // Allocate storage for xy faces attached to columns.
  xy_data->column_xy_face_offsets = polymec_malloc(sizeof(int) * (xy_data->num_columns+1));
  xy_data->column_xy_face_offsets[0] = 0;
  for (int c = 0; c < xy_data->num_columns; ++c)
  {
    int cell = local_cells->data[c];
    xy_data->column_xy_face_offsets[c+1] = xy_data->column_xy_face_offsets[c] + 
                                           planar_polymesh_cell_num_edges(mesh, cell);
  }
  xy_data->column_xy_faces = polymec_malloc(sizeof(int) * xy_data->column_xy_face_offsets[xy_data->num_columns]);

  // Here we map planar mesh nodes to chunk nodes.
  int_int_unordered_map_t* node_map = int_int_unordered_map_new(); 

  // Now determine the connectivity in a single pass.
  int_array_t* face_cols = int_array_new();
  int_array_t* edge_nodes = int_array_new();
  for (int edge = 0; edge < mesh->num_edges; ++edge)
  {
    int cell1 = mesh->edge_cells[2*edge];
    int cell2 = mesh->edge_cells[2*edge+1];

    int* col1_p = int_int_unordered_map_get(cell_to_col_map, cell1);
    int* col2_p = int_int_unordered_map_get(cell_to_col_map, cell2);
    if ((col1_p != NULL) || (col2_p != NULL)) // at least one column is present locally
    {
      // The locally present column is column 1, by definition.
      int col1;
      bool cols_reversed = false;
      if (col1_p != NULL) 
        col1 = *col1_p; 
      else
      { 
        col1 = *col2_p;
        cols_reversed = true;
      }
      int col2 = ((col1_p != NULL) && (col2_p != NULL)) ? *col2_p : -1;
      if (col2 == -1) // ghost column
        col2 = (int)(xy_data->num_columns + xy_data->receive_map->size); 

      // For the interior columns we identify which face this edge corresponds to.
      int face = (int)(xy_data->num_xy_faces);
      int f1 = 0;
      int f2 = -1;
      {
        int cell = (cols_reversed) ? cell2 : cell1;
        while ((mesh->cell_edges[mesh->cell_edge_offsets[cell]+f1] != edge) && 
               (mesh->cell_edges[mesh->cell_edge_offsets[cell]+f1] != ~edge))
          ++f1;
        ASSERT(mesh->cell_edge_offsets[cell+1] > (mesh->cell_edge_offsets[cell] + f1));
      }
      if (col2 < xy_data->num_columns)
      {
        f2 = 0;
        int cell = (cols_reversed) ? cell1 : cell2;
        while ((mesh->cell_edges[mesh->cell_edge_offsets[cell]+f2] != edge) &&
               (mesh->cell_edges[mesh->cell_edge_offsets[cell]+f2] != ~edge)) 
          ++f2;
        ASSERT(mesh->cell_edge_offsets[cell+1] > (mesh->cell_edge_offsets[cell] + f2));
      }
      else if (col2 != -1) // ghost column
      {
        // Set up a parallel exchange between the columns.
        int cell = (cols_reversed) ? cell2 : cell1;
        int neighbor_xy_index = (int)partition_vector[cell];
        exchanger_proc_map_add_index(xy_data->send_map, neighbor_xy_index, col1);
        exchanger_proc_map_add_index(xy_data->receive_map, neighbor_xy_index, col2);
        ++xy_data->num_ghost_columns;
      }

      // Create a new xy face for this edge, and hook it up to columns 
      // corresponding to the edge's adjacent planar cells.
      xy_data->column_xy_faces[xy_data->column_xy_face_offsets[col1]+f1] = face;
      int_array_append(face_cols, col1);
      if (f2 != -1)
        xy_data->column_xy_faces[xy_data->column_xy_face_offsets[col2]+f2] = ~face;
      int_array_append(face_cols, col2);
      ++xy_data->num_xy_faces;

      // Read off the 2 nodes connecting the edge for the planar cell and 
      // see if we've already added them.
      int n1 = mesh->edge_nodes[2*edge];
      int n2 = mesh->edge_nodes[2*edge+1];
      int node1, node2;
      int* node1_p = int_int_unordered_map_get(node_map, n1);
      if (node1_p == NULL)
      {
        node1 = node_map->size;
        int_int_unordered_map_insert(node_map, n1, node1);
      }
      else
        node1 = *node1_p;
      int* node2_p = int_int_unordered_map_get(node_map, n2);
      if (node2_p == NULL)
      {
        node2 = node_map->size;
        int_int_unordered_map_insert(node_map, n2, node2);
      }
      else
        node2 = *node2_p;

      // Now create an xy edge connecting these two nodes.
      int_array_append(edge_nodes, node1);
      int_array_append(edge_nodes, node2);
    }
  }

  // Surrender the data from the various arrays.
  ASSERT(xy_data->num_xy_faces == (face_cols->size / 2));
  xy_data->xy_face_columns = face_cols->data;
  int_array_release_data_and_free(face_cols);

  xy_data->xy_edge_nodes = edge_nodes->data;
  xy_data->num_xy_edges = edge_nodes->size/2;
  int_array_release_data_and_free(edge_nodes);

  // Set node positions.
  xy_data->num_xy_nodes = node_map->size;
  xy_data->xy_nodes = polymec_malloc(sizeof(point2_t) * node_map->size);
  int npos = 0, n, node;
  while (int_int_unordered_map_next(node_map, &npos, &n, &node))
    xy_data->xy_nodes[node] = mesh->nodes[n];
  int_int_unordered_map_free(node_map);

  // Clean up.
  int_int_unordered_map_free(cell_to_col_map);
  int_array_free(local_cells);

  return xy_data;
}

static void chunk_xy_data_free(chunk_xy_data_t* xy_data)
{
  polymec_free(xy_data->xy_nodes);
  polymec_free(xy_data->xy_edge_nodes);
  polymec_free(xy_data->xy_face_columns);
  polymec_free(xy_data->column_xy_faces);
  polymec_free(xy_data->column_xy_face_offsets);
  exchanger_proc_map_free(xy_data->send_map);
  exchanger_proc_map_free(xy_data->receive_map);
  polymec_free(xy_data);
}

static void free_chunk(colmesh_chunk_t* chunk)
{
  // The chunk's xy data is managed by the colmesh.
  polymec_free(chunk);
}

DEFINE_UNORDERED_MAP(chunk_map, int, colmesh_chunk_t*, int_hash, int_equals)

struct colmesh_t 
{
  MPI_Comm comm;
  int nproc, rank;

  size_t num_xy_chunks, num_z_chunks, nz_per_chunk;
  real_t z1, z2;

  // Chunk data.
  adj_graph_t* chunk_graph;
  chunk_xy_data_array_t* chunk_xy_data;
  chunk_map_t* chunks;
  int* chunk_indices;

  // Exchangers.
  exchanger_t* cell_ex;
  exchanger_t* xy_face_ex;
  exchanger_t* z_face_ex;
  exchanger_t* xy_edge_ex;
  exchanger_t* z_edge_ex;
  exchanger_t* node_ex;

  // True if the mesh is periodic along the z axis, false if not.
  bool periodic_in_z;

  // This flag is set by colmesh_finalize() after a mesh has been assembled.
  bool finalized;
};

static inline int chunk_index(colmesh_t* mesh, int xy_index, int z_index)
{
  return (int)(mesh->num_z_chunks * xy_index + z_index);
}

colmesh_t* create_empty_colmesh(MPI_Comm comm, 
                                planar_polymesh_t* columns,
                                real_t z1, real_t z2,
                                size_t num_xy_chunks, size_t num_z_chunks,
                                size_t nz_per_chunk, bool periodic_in_z)
{
  ASSERT(columns != NULL);
  ASSERT(z1 < z2);
  ASSERT(num_xy_chunks > 0);
  ASSERT(num_z_chunks > 0);
  ASSERT(nz_per_chunk > 0);

  colmesh_t* mesh = polymec_malloc(sizeof(colmesh_t));
  mesh->comm = comm;
  mesh->chunks = chunk_map_new();
  mesh->chunk_indices = NULL;
  mesh->num_xy_chunks = num_xy_chunks;
  mesh->num_z_chunks = num_z_chunks;
  mesh->nz_per_chunk = nz_per_chunk;
  mesh->z1 = z1;
  mesh->z2 = z2;
  MPI_Comm_size(comm, &mesh->nproc);
  MPI_Comm_rank(comm, &mesh->rank);
  mesh->periodic_in_z = periodic_in_z;
  mesh->finalized = false;

  // Partition the planar polymesh.
#if POLYMEC_HAVE_MPI
  adj_graph_t* planar_graph = graph_from_planar_polymesh_cells(columns);
  int64_t* P = partition_graph_n_ways(planar_graph, num_xy_chunks, NULL, 0.05);
  adj_graph_free(planar_graph);
#else
  // Use a space-filling curve, since we haven't built graph cutting in.
  point_t* centroids = polymec_malloc(sizeof(point_t) * columns->num_cells);
  polygon_t* poly = planar_polymesh_cell_polygon(columns, 0);
  for (int i = 0; i < columns->num_cells; ++i)
  {
    planar_polymesh_cell_get_polygon(columns, i, poly);
    point2_t xc;
    polygon_compute_centroid(poly, &xc);
    centroids[i].x = xc.x;
    centroids[i].y = xc.y;
    centroids[i].z = 0.0;
  }
  int64_t* P = partition_points_n_ways(centroids, (size_t)columns->num_cells, num_xy_chunks, NULL, 0.0);
  polymec_free(centroids);
  polymec_release(poly);
#endif

  // Create xy data for chunks, and create a global graph whose vertices 
  // are the mesh's chunks.
  mesh->chunk_xy_data = chunk_xy_data_array_new();
  size_t num_chunks = num_xy_chunks * num_z_chunks;
  mesh->chunk_graph = adj_graph_new(MPI_COMM_SELF, num_chunks);
  for (int xy_index = 0; xy_index < (int)num_xy_chunks; ++xy_index)
  {
    chunk_xy_data_t* xy_data = chunk_xy_data_new(comm, columns, P, xy_index);
    chunk_xy_data_array_append_with_dtor(mesh->chunk_xy_data, xy_data, 
                                         chunk_xy_data_free);

    // Extract the neighboring chunks from this xy data and add edges to our graph.
    for (int z_index = 0; z_index < (int)num_z_chunks; ++z_index)
    {
      int ch_index = chunk_index(mesh, xy_index, z_index);

      // Count up edges.
      int num_xy_neighbors = (int)xy_data->send_map->size;
      int num_z_neighbors;
      if (num_z_chunks == 1)
        num_z_neighbors = 0;
      else if ((z_index == 0) || (z_index == (int)(num_z_chunks-1)))
        num_z_neighbors = 1;
      else 
        num_z_neighbors = 2;
      adj_graph_set_num_edges(mesh->chunk_graph, ch_index, num_xy_neighbors + num_z_neighbors);

      // Add xy edges.
      int* edges = adj_graph_edges(mesh->chunk_graph, ch_index);
      int pos = 0, neighbor_xy_index, i = 0;
      int_array_t* indices;
      while (exchanger_proc_map_next(xy_data->send_map, &pos, &neighbor_xy_index, &indices))
        edges[i++] = chunk_index(mesh, neighbor_xy_index, z_index);

      // Add z edges.
      if (z_index > 0)
        edges[i++] = chunk_index(mesh, xy_index, z_index-1);
      if (z_index < (int)(num_z_chunks-1))
        edges[i++] = chunk_index(mesh, xy_index, z_index+1);

      ASSERT(i == num_xy_neighbors + num_z_neighbors);
    }
  }
  polymec_free(P);

  // No exchangers yet.
  mesh->cell_ex = NULL;
  mesh->xy_face_ex = NULL;
  mesh->z_face_ex = NULL;
  mesh->xy_edge_ex = NULL;
  mesh->z_edge_ex = NULL;
  mesh->node_ex = NULL;

  return mesh;
}

void colmesh_insert_chunk(colmesh_t* mesh, int xy_index, int z_index)
{
  ASSERT(!mesh->finalized);
  ASSERT(xy_index >= 0);
  ASSERT(xy_index < mesh->num_xy_chunks);
  ASSERT(z_index >= 0);
  ASSERT(z_index < mesh->num_z_chunks);
  ASSERT(!colmesh_has_chunk(mesh, xy_index, z_index));

  colmesh_chunk_t* chunk = polymec_malloc(sizeof(colmesh_chunk_t));

  // Chunk xy data (points to data owned by mesh->chunk_xy_data).
  chunk_xy_data_t* xy_data = mesh->chunk_xy_data->data[xy_index];
  chunk->num_columns = xy_data->num_columns;
  chunk->num_ghost_columns = xy_data->num_ghost_columns;
  chunk->column_xy_face_offsets = xy_data->column_xy_face_offsets;
  chunk->column_xy_faces = xy_data->column_xy_faces;
  chunk->num_xy_faces = xy_data->num_xy_faces;
  chunk->xy_face_columns = xy_data->xy_face_columns;
  chunk->num_xy_edges = xy_data->num_xy_edges;
  chunk->num_xy_nodes = xy_data->num_xy_nodes;
  chunk->xy_nodes = xy_data->xy_nodes;

  // Chunk z data.
  size_t nz = mesh->nz_per_chunk * mesh->num_z_chunks;
  real_t dz = (mesh->z2 - mesh->z1) / nz;
  real_t z1 = mesh->z1 + z_index * dz;
  real_t z2 = mesh->z1 + (z_index+1) * dz;
  chunk->num_z_cells = mesh->nz_per_chunk;
  chunk->z1 = z1;
  chunk->z2 = z2;

  // Add this chunk to our list.
  int index = chunk_index(mesh, xy_index, z_index);
  chunk_map_insert_with_v_dtor(mesh->chunks, index, chunk, free_chunk);
}

#if POLYMEC_HAVE_MPI
static int64_t* source_vector(colmesh_t* mesh)
{
  // Catalog all the chunks on this process.
  int_array_t* my_chunks = int_array_new();
  for (int xy_index = 0; xy_index < (int)mesh->num_xy_chunks; ++xy_index)
  {
    for (int z_index = 0; z_index < (int)mesh->num_z_chunks; ++z_index)
    {
      if (colmesh_has_chunk(mesh, xy_index, z_index))
        int_array_append(my_chunks, chunk_index(mesh, xy_index, z_index));
    }
  }

  // Gather the numbers of chunks owned by each process.
  int num_my_chunks = (int)my_chunks->size;
  int num_chunks_for_proc[mesh->nproc];
  MPI_Allgather(&num_my_chunks, 1, MPI_INT, 
                num_chunks_for_proc, 1, MPI_INT, mesh->comm);

  // Arrange for the storage of the chunk indices for the patches stored 
  // on each process.
  int proc_offsets[mesh->nproc+1];
  proc_offsets[0] = 0;
  for (int p = 0; p < mesh->nproc; ++p)
    proc_offsets[p+1] = proc_offsets[p] + num_chunks_for_proc[p];

  // Gather the indices of the chunks owned by all processes into a huge list.
  int num_all_chunks = (int)(mesh->num_xy_chunks * mesh->num_z_chunks);
  ASSERT(num_all_chunks == proc_offsets[mesh->nproc]);
  int* all_chunks = polymec_malloc(sizeof(int) * num_all_chunks);
  MPI_Allgatherv(my_chunks->data, num_my_chunks, MPI_INT, 
                 all_chunks, num_chunks_for_proc, proc_offsets,
                 MPI_INT, mesh->comm);

  // Clean up a bit.
  int_array_free(my_chunks);

  // Convert the huge list into a source vector.
  int64_t* sources = polymec_malloc(sizeof(int64_t) * num_all_chunks);
  for (int p = 0; p < mesh->nproc; ++p)
  {
    for (int offset = proc_offsets[p]; offset < proc_offsets[p+1]; ++offset)
      sources[all_chunks[offset]] = (int64_t)p;
  }

  polymec_free(all_chunks);
  return sources;
}
#endif

void colmesh_finalize(colmesh_t* mesh)
{
  ASSERT(!mesh->finalized);

  // Create a sorted list of chunk indices and compute chunk index offsets.
  mesh->chunk_indices = polymec_malloc(sizeof(int) * 2 * mesh->chunks->size);
  int chunk_offsets[mesh->chunks->size+1];
  chunk_offsets[0] = 0;
  {
    size_t k = 0;
    for (size_t i = 0; i < mesh->num_xy_chunks; ++i)
    {
      size_t num_z_chunks = 0;
      for (size_t j = 0; j < mesh->num_z_chunks; ++j)
      {
        int index = chunk_index(mesh, (int)i, (int)j);
        colmesh_chunk_t** chunk_p = chunk_map_get(mesh->chunks, index);
        if (chunk_p != NULL)
        {
          colmesh_chunk_t* chunk    = *chunk_p;
          chunk_offsets[k+1]         = chunk_offsets[k] + 
                                       (int)(chunk->num_columns * chunk->num_z_cells);
          mesh->chunk_indices[2*k]   = (int)i;
          mesh->chunk_indices[2*k+1] = (int)j;
          ++k;
          ++num_z_chunks;
        }
      }
    }
    ASSERT(k == mesh->chunks->size);
  }

#if POLYMEC_HAVE_MPI
  // Figure out which processes own what chunks.
  int64_t* owners = source_vector(mesh);
#else
  int num_all_chunks = (int)(mesh->num_xy_chunks * mesh->num_z_chunks);
  int64_t* owners = polymec_calloc(sizeof(int) * num_all_chunks);
#endif

  // Assemble a cell exchanger by traversing the locally stored chunks.
  exchanger_proc_map_t* send_map = exchanger_proc_map_new();
  exchanger_proc_map_t* receive_map = exchanger_proc_map_new();
  for (size_t i = 0; i < mesh->chunks->size; ++i)
  {
    int xy = mesh->chunk_indices[2*i];
    int z = mesh->chunk_indices[2*i+1];
    int ch_index = chunk_index(mesh, xy, z);
    colmesh_chunk_t* chunk = *chunk_map_get(mesh->chunks, ch_index);
    int chunk_offset = chunk_offsets[i];

    chunk_xy_data_t* xy_data = mesh->chunk_xy_data->data[xy];
    exchanger_proc_map_t* xy_sends = xy_data->send_map;
    exchanger_proc_map_t* xy_receives = xy_data->receive_map;

    // Hook up the xy sends/receives.
    int pos = 0, neighbor_xy_index;
    int_array_t* indices;
    while (exchanger_proc_map_next(xy_sends, &pos, &neighbor_xy_index, &indices))
    {
      for (size_t k = 0; k < indices->size; ++k)
      {
        int proc = (int)(owners[neighbor_xy_index]);
        int index = (int)(chunk_offset + chunk->num_z_cells * indices->data[k] + z);
        exchanger_proc_map_add_index(send_map, proc, index);
      }
    }
    pos = 0;
    int cell_offset = 0;
    while (exchanger_proc_map_next(xy_receives, &pos, &neighbor_xy_index, &indices))
    {
      for (size_t k = 0; k < indices->size; ++k)
      {
        int proc = (int)(owners[neighbor_xy_index]);
        int index = (int)(chunk_offset + chunk->num_z_cells * cell_offset + z);
        exchanger_proc_map_add_index(receive_map, proc, index);
        ++cell_offset;
      }
    }

    // Now hook up the z sends/receives.
    if (z > 0)
    {
      for (size_t xy1 = 0; xy1 < xy_data->num_columns; ++xy1)
      {
        int ch1_index = chunk_index(mesh, (int)xy1, z-1);
        int proc = (int)(owners[ch1_index]);
        int z1 = 1;
        int send_index = (int)(chunk_offset + chunk->num_z_cells * xy1 + z1);
        int receive_index = cell_offset; 
        exchanger_proc_map_add_index(send_map, proc, send_index);
        exchanger_proc_map_add_index(receive_map, proc, receive_index);
        ++cell_offset;
      }
    }
    if (z < (mesh->num_z_chunks-1))
    {
      for (size_t xy1 = 0; xy1 < xy_data->num_columns; ++xy1)
      {
        int ch1_index = chunk_index(mesh, (int)xy1, z+1);
        int proc = (int)(owners[ch1_index]);
        int z1 = chunk->num_z_cells;
        int send_index = (int)(chunk_offset + chunk->num_z_cells * xy1 + z1);
        int receive_index = cell_offset; 
        exchanger_proc_map_add_index(send_map, proc, send_index);
        exchanger_proc_map_add_index(receive_map, proc, receive_index);
        ++cell_offset;
      }
    }
  }
  mesh->cell_ex = exchanger_new(mesh->comm);
  exchanger_set_sends(mesh->cell_ex, send_map);
  exchanger_set_receives(mesh->cell_ex, receive_map);
  polymec_free(owners);

  // Prune unused xy data.
  for (size_t i = 0; i < mesh->num_xy_chunks; ++i)
  {
    size_t num_z_chunks = 0;
    for (size_t j = 0; j < mesh->num_z_chunks; ++j)
    {
      int index = chunk_index(mesh, (int)i, (int)j);
      if (chunk_map_contains(mesh->chunks, index))
        ++num_z_chunks;
    }

    // Prune xy data for this xy_index if we don't have any chunks here.
    if (num_z_chunks == 0)
    {
      // Punch out this data on this process, and clear the corresponding 
      // destructor so it doesn't get double-freed.
      chunk_xy_data_free(mesh->chunk_xy_data->data[i]);
      mesh->chunk_xy_data->data[i] = NULL;
      mesh->chunk_xy_data->dtors[i] = NULL;
    }
  }

  // We're finished here.
  log_debug("colmesh_finalize: finalized with %d local chunks.", mesh->chunks->size);
  mesh->finalized = true;
}

colmesh_t* colmesh_new(MPI_Comm comm,
                       planar_polymesh_t* columns,
                       real_t z1, real_t z2,
                       size_t nz, bool periodic_in_z)
{
  size_t num_xy_chunks = 1;
  size_t num_z_chunks = 1; 
  int nproc, rank;
  MPI_Comm_size(comm, &nproc);
  MPI_Comm_rank(comm, &rank);
  if (nproc > 1)
  {
    // We minimize (nz*nz - nxy) subject to the constraint
    // (Nz/nz)*(Nxy/nxy) = num_z_chunks * num_xy_chunks = nproc, where
    //  nz is the number of z cells in a chunk, 
    //  nxy is the number of xy cells in a chunk,
    //  Nz is the number of z cells in the mesh, and 
    //  Nxy is the number of xy cells in the entire xy plane. 
    int Nz = (int)nz;
    int Nxy = (int)columns->num_cells;

    // Some simple algebra tell us num_z_chunks should be the integer closest 
    // to the value pow(nproc*Nz*Nz/Nxy, 1.0/3.0).
    num_z_chunks = MAX(1, (size_t)(pow(nproc*Nz*Nz/Nxy, 1.0/3.0)));
    
    // Adjust num_z_chunks so it evenly divides nproc.
    int remainder = (int)(nproc % num_z_chunks);
    if (remainder < nproc/2)
      num_z_chunks -= remainder;
    else
      num_z_chunks += remainder;
    num_xy_chunks = nproc / num_z_chunks;
    log_info("colmesh_new: Dividing %d x %d cells up into %d x %d chunks.",
             Nxy, Nz, num_xy_chunks, num_z_chunks);
  }
  ASSERT(num_xy_chunks * num_z_chunks == nproc);

  // Now create an empty colmesh with the desired numbers of chunks, and 
  // insert all the chunks on each process. We do a "naive" placement of 
  // the chunks by allocating them sequentially to processes in a flattened 
  // index space I(xy_index, z_index) = num_z_chunks * xy_index + z_index.
  // This is definitely not ideal, but it's the easiest way to get a start.
  size_t nz_per_chunk = nz / num_z_chunks;
  colmesh_t* mesh = create_empty_colmesh(comm, columns, z1, z2, 
                                         num_xy_chunks, num_z_chunks, 
                                         nz_per_chunk, periodic_in_z);
  size_t tot_num_chunks = num_xy_chunks * num_z_chunks;
  size_t chunks_per_proc = tot_num_chunks / nproc;
  for (size_t i = 0; i < tot_num_chunks; ++i)
  {
    if ((i >= rank*chunks_per_proc) && (i < (rank+1)*chunks_per_proc))
    {
      int xy_index = (int)(i / num_z_chunks);
      int z_index = (int)(i % num_z_chunks);
      log_debug("colmesh_new: Inserting chunk (%d, %d).", xy_index, z_index);
      colmesh_insert_chunk(mesh, xy_index, z_index);
    }
  }

  // Put the lid on it and ship it.
  colmesh_finalize(mesh);
  return mesh;
}

void colmesh_free(colmesh_t* mesh)
{
  polymec_free(mesh->chunk_indices);
  chunk_map_free(mesh->chunks);
  chunk_xy_data_array_free(mesh->chunk_xy_data);
  adj_graph_free(mesh->chunk_graph);
  if (mesh->cell_ex != NULL)
    polymec_release(mesh->cell_ex);
  if (mesh->xy_face_ex != NULL)
    polymec_release(mesh->xy_face_ex);
  if (mesh->z_face_ex != NULL)
    polymec_release(mesh->z_face_ex);
  if (mesh->xy_edge_ex != NULL)
    polymec_release(mesh->xy_edge_ex);
  if (mesh->z_edge_ex != NULL)
    polymec_release(mesh->z_edge_ex);
  if (mesh->node_ex != NULL)
    polymec_release(mesh->node_ex);
  polymec_free(mesh);
}

void colmesh_get_chunk_info(colmesh_t* mesh, 
                            size_t* num_xy_chunks,
                            size_t* num_z_chunks,
                            size_t* nz_per_chunk)
{
  *num_xy_chunks = mesh->num_xy_chunks;
  *num_z_chunks = mesh->num_z_chunks;
  *nz_per_chunk = mesh->nz_per_chunk;
}

void colmesh_get_z_info(colmesh_t* mesh, 
                        real_t* z1,
                        real_t* z2,
                        bool* periodic)
{
  *z1 = mesh->z1;
  *z2 = mesh->z2;
  *periodic = mesh->periodic_in_z;
}

bool colmesh_verify_topology(colmesh_t* mesh, 
                             void (*handler)(const char* format, ...))
{
  // Verify the topology of each chunk.
  bool good = true;
  colmesh_chunk_t* chunk;
  int pos = 0, xy_index, z_index;
  while (colmesh_next_chunk(mesh, &pos, &xy_index, &z_index, &chunk))
  {
    bool result = colmesh_chunk_verify_topology(chunk, handler);
    if (!result)
    {
      good = false;
      break;
    }
  }
  return good;
}

bool colmesh_chunk_verify_topology(colmesh_chunk_t* chunk,
                                   void (*handler)(const char* format, ...))
{
  // All cells are prisms and must have at least 5 faces (3 xy faces + 2 z faces).
  for (size_t c = 0; c < chunk->num_columns; ++c)
  {
    if (colmesh_chunk_column_num_xy_faces(chunk, c) < 3)
    {
      handler("polymesh_verify_topology: column %d has only %d faces per cell.", 
              (int)c, colmesh_chunk_column_num_xy_faces(chunk, c) + 2);
      return false;
    }
  }

  // All z faces must have at least 3 nodes/edges.
  int num_z_faces = (int)chunk->num_columns;
  for (int f = 0; f < num_z_faces; ++f)
  {
    int nn = colmesh_chunk_z_face_num_nodes(chunk, f);
    if (nn == 0)
    {
      handler("colmesh_verify_topology: column %d has a polygonal z face with no edges!", f);
      return false;
    }
    if (nn < 3)
    {
      handler("colmesh_verify_topology: column %d has a polygonal face %d with only %d nodes.", f, nn);
      return false;
    }
  }

  // Make sure that all the xy faces attached to this cell have it in their list.
  for (size_t c = 0; c < chunk->num_columns; ++c)
  {
    int num_xy_faces = colmesh_chunk_column_num_xy_faces(chunk, c);
    int xy_faces[num_xy_faces];
    colmesh_chunk_column_get_xy_faces(chunk, c, xy_faces);
    for (int f = 0; f < num_xy_faces; ++f)
    {
      int face = xy_faces[f];
      if ((chunk->xy_face_columns[2*face] != c) && 
          (chunk->xy_face_columns[2*face+1] != c))
      {
        handler("colmesh_verify_topology: column %d has xy face %d in its list "
                "of faces, but that face does not have that column in its list.", c, xy_faces[f]);
        return false;
      }
    }
  }

  // Now go over all xy faces and make sure that their columns can see them, too.
  for (size_t f = 0; f < chunk->num_xy_faces; ++f)
  {
    int column = chunk->xy_face_columns[2*f];
    int num_xy_faces = colmesh_chunk_column_num_xy_faces(chunk, column);
    int xy_faces[num_xy_faces];
    colmesh_chunk_column_get_xy_faces(chunk, column, xy_faces);
    bool found_face = false;
    for (int ff = 0; ff < num_xy_faces; ++ff)
    {
      if (xy_faces[ff] == f) 
      {
        found_face = true;
        break;
      }
    }
    if (!found_face)
    {
      handler("colmesh_verify_topology: xy face %d has column %d in its list of columns, but "
              "that column does not have that face in its list.", f, chunk->xy_face_columns[2*f]);
      return false;
    }

    // Check the column on the other side, too (but only if its a non-ghost column).
    if ((chunk->xy_face_columns[2*f+1] != -1) && (chunk->xy_face_columns[2*f+1] < chunk->num_columns))
    {
      found_face = false;
      int other_col = chunk->xy_face_columns[2*f+1];
      num_xy_faces = colmesh_chunk_column_num_xy_faces(chunk, other_col);
      colmesh_chunk_column_get_xy_faces(chunk, other_col, xy_faces);
      for (int ff = 0; ff < num_xy_faces; ++ff)
      {
        if (xy_faces[ff] == f) 
        {
          found_face = true;
          break;
        }
      }
      if (!found_face)
      {
        handler("colmesh_verify_topology: xy face %d has column %d in its list of columns, "
                "but that column does not have that xy face in its list of faces.", f, other_col);
        return false;
      }
    }
  }
  return true;
}

MPI_Comm colmesh_comm(colmesh_t* mesh)
{
  return mesh->comm;
}

size_t colmesh_num_chunks(colmesh_t* mesh)
{
  return mesh->chunks->size;
}

colmesh_chunk_t* colmesh_chunk(colmesh_t* mesh, int xy_index, int z_index)
{
  int index = chunk_index(mesh, xy_index, z_index);
  colmesh_chunk_t** chunk_p = chunk_map_get(mesh->chunks, index);
  return (chunk_p != NULL) ? *chunk_p : NULL;
}

bool colmesh_has_chunk(colmesh_t* mesh, int xy_index, int z_index)
{
  return (colmesh_chunk(mesh, xy_index, z_index) != NULL);
}

polygon_t* colmesh_chunk_polygon(colmesh_chunk_t* chunk, int column)
{
  ASSERT(column < (int)chunk->num_columns);
  int z_face = column;
  int num_nodes = colmesh_chunk_z_face_num_nodes(chunk, z_face);
  int nodes[num_nodes];
  colmesh_chunk_z_face_get_nodes(chunk, z_face, nodes);
  point2_t vertices[num_nodes];
  for (int n = 0; n < num_nodes; ++n)
    vertices[n] = chunk->xy_nodes[nodes[n]];
  return polygon_new(vertices, num_nodes);
}

bool colmesh_next_chunk(colmesh_t* mesh, int* pos, 
                        int* xy_index, int* z_index,
                        colmesh_chunk_t** chunk)
{
  if (*pos >= (int)mesh->chunks->size)
    return false;
  else
  {
    int i= *pos;
    *xy_index = mesh->chunk_indices[2*i];
    *z_index = mesh->chunk_indices[2*i+1];
    int index = chunk_index(mesh, *xy_index, *z_index);
    if (chunk != NULL)
      *chunk = *chunk_map_get(mesh->chunks, index);
    ++(*pos);
    return true;
  }
}

void colmesh_chunk_get_node(colmesh_chunk_t* chunk, 
                            int xy_index, int z_index,
                            point_t* node_pos)
{
  ASSERT(xy_index >= 0);
  ASSERT(xy_index < chunk->num_xy_nodes);
  ASSERT(z_index >= 0);
  ASSERT(z_index <= chunk->num_z_cells);
  point2_t* xy = &(chunk->xy_nodes[xy_index]);
  node_pos->x = xy->x;
  node_pos->y = xy->y;
  real_t dz = (chunk->z2 - chunk->z1) / chunk->num_z_cells;
  node_pos->z = chunk->z1 + z_index * dz;
}

#if POLYMEC_HAVE_MPI
static size_t xy_data_byte_size(void* obj)
{
  chunk_xy_data_t* xy_data = obj;
  size_t size = 5 * sizeof(size_t) + 
                (xy_data->num_columns+1) * sizeof(int) + 
                xy_data->column_xy_face_offsets[xy_data->num_columns] * sizeof(int) + 
                2 * xy_data->num_xy_faces * sizeof(int) +
                2 * xy_data->num_xy_edges * sizeof(int) + 
                xy_data->num_xy_nodes * sizeof(point2_t);
  return size;
}

static void* xy_data_byte_read(byte_array_t* bytes, size_t* offset)
{
  chunk_xy_data_t* xy_data = polymec_malloc(sizeof(chunk_xy_data_t));
  byte_array_read_size_ts(bytes, 1, &xy_data->num_columns, offset);
  byte_array_read_size_ts(bytes, 1, &xy_data->num_ghost_columns, offset);
  xy_data->column_xy_face_offsets = polymec_malloc((xy_data->num_columns + 1) * sizeof(int));
  byte_array_read_ints(bytes, xy_data->num_columns+1, xy_data->column_xy_face_offsets, offset);
  xy_data->column_xy_faces = polymec_malloc(xy_data->column_xy_face_offsets[xy_data->num_columns] * sizeof(int));
  byte_array_read_ints(bytes, xy_data->column_xy_faces[xy_data->column_xy_face_offsets[xy_data->num_columns]], xy_data->column_xy_faces, offset);
  byte_array_read_size_ts(bytes, 1, &xy_data->num_xy_faces, offset);
  xy_data->xy_face_columns = polymec_malloc(2*xy_data->num_xy_faces * sizeof(int));
  byte_array_read_ints(bytes, 2*xy_data->num_xy_faces, xy_data->xy_face_columns, offset);
  byte_array_read_size_ts(bytes, 1, &xy_data->num_xy_edges, offset);
  xy_data->xy_edge_nodes = polymec_malloc(2*xy_data->num_xy_edges * sizeof(int));
  byte_array_read_ints(bytes, 2*xy_data->num_xy_edges, xy_data->xy_edge_nodes, offset);
  byte_array_read_size_ts(bytes, 1, &xy_data->num_xy_nodes, offset);
  xy_data->xy_nodes = polymec_malloc(xy_data->num_xy_nodes * sizeof(point2_t));
  byte_array_read_point2s(bytes, xy_data->num_xy_nodes, xy_data->xy_nodes, offset);
  return xy_data;
}

static void xy_data_byte_write(void* obj, byte_array_t* bytes, size_t* offset)
{
  chunk_xy_data_t* xy_data = obj;
  byte_array_write_size_ts(bytes, 1, &xy_data->num_columns, offset);
  byte_array_write_size_ts(bytes, 1, &xy_data->num_ghost_columns, offset);
  byte_array_write_ints(bytes, xy_data->num_columns+1, xy_data->column_xy_face_offsets, offset);
  byte_array_write_ints(bytes, xy_data->column_xy_faces[xy_data->column_xy_face_offsets[xy_data->num_columns]], xy_data->column_xy_faces, offset);
  byte_array_write_size_ts(bytes, 1, &xy_data->num_xy_faces, offset);
  byte_array_write_ints(bytes, 2*xy_data->num_xy_faces, xy_data->xy_face_columns, offset);
  byte_array_write_size_ts(bytes, 1, &xy_data->num_xy_edges, offset);
  byte_array_write_ints(bytes, 2*xy_data->num_xy_edges, xy_data->xy_edge_nodes, offset);
  byte_array_write_size_ts(bytes, 1, &xy_data->num_xy_nodes, offset);
  byte_array_write_point2s(bytes, xy_data->num_xy_nodes, xy_data->xy_nodes, offset);
}

static serializer_t* chunk_xy_data_serializer()
{
  return serializer_new("chunk_xy_data", xy_data_byte_size, xy_data_byte_read, xy_data_byte_write, NULL);
}

// This helper inserts an integer into an array in sorted order if it isn't 
// already in the array. If it's already there, this function has no effect.
static void insert_unique_sorted(int_array_t* array, int value)
{
  size_t index = int_lower_bound(array->data, array->size, value);
  if ((index >= array->size) || (array->data[index] != value))
    int_array_insert(array, index, value);
}

static chunk_xy_data_array_t* redistribute_chunk_xy_data(colmesh_t* old_mesh,
                                                         int64_t* partition_vector,
                                                         int64_t* source_vector)
{
  START_FUNCTION_TIMER();

  // Who are we receiving from, and who are we sending to?
  int_array_t* source_procs = int_array_new();
  int_array_t* dest_procs = int_array_new();
  for (size_t i = 0; i < old_mesh->num_xy_chunks; ++i)
  {
    if (partition_vector[i] == old_mesh->rank)
      insert_unique_sorted(source_procs, (int)(source_vector[i]));
    if (source_vector[i] == old_mesh->rank)
      insert_unique_sorted(dest_procs, (int)(partition_vector[i]));
  }

  // Pack up the xy data we're sending out.
  serializer_t* ser = chunk_xy_data_serializer();
  byte_array_t* send_buffers[dest_procs->size];
  for (size_t i = 0; i < dest_procs->size; ++i)
  {
    int proc = dest_procs->data[i];
    size_t offset = 0;
    byte_array_t* buffer = byte_array_new();
    for (size_t j = 0; j < old_mesh->num_xy_chunks; ++j)
    {
      if ((source_vector[j] == old_mesh->rank) && (partition_vector[j] == proc))
        serializer_write(ser, old_mesh->chunk_xy_data->data[j], buffer, &offset);
    }
    send_buffers[i] = buffer;
  }

  // Post receives for incoming data sizes.
  MPI_Request requests[source_procs->size + dest_procs->size];
  int receive_sizes[source_procs->size];
  for (size_t i = 0; i < source_procs->size; ++i)
  {
    int proc = source_procs->data[i];
    int err = MPI_Irecv(&receive_sizes[i], 1, MPI_INT, proc, 0, old_mesh->comm, &requests[i]);
    if (err != MPI_SUCCESS)
      polymec_error("Error receiving xy data size from rank %d", proc);
  }

  // Send out data sizes.
  for (size_t i = 0; i < dest_procs->size; ++i)
  {
    int proc = dest_procs->data[i];
    int size = (int)(send_buffers[i]->size);
    int err = MPI_Isend(&size, 1, MPI_INT, proc, 0, old_mesh->comm, &requests[source_procs->size + i]);
    if (err != MPI_SUCCESS)
      polymec_error("Error sending xy data size from rank %d", proc);
  }

  // Wait for everything to finish.
  MPI_Waitall((int)(source_procs->size + dest_procs->size), requests, MPI_STATUSES_IGNORE);

  // Post receives for incoming xy data.
  byte_array_t* receive_buffers[source_procs->size];
  for (size_t i = 0; i < source_procs->size; ++i)
  {
    int proc = source_procs->data[i];
    byte_array_t* buffer = byte_array_new_with_size(receive_sizes[i]);
    int err = MPI_Irecv(buffer->data, receive_sizes[i], MPI_BYTE, proc, 0, old_mesh->comm, &requests[i]);
    if (err != MPI_SUCCESS)
      polymec_error("Error receiving xy data from rank %d", proc);
    receive_buffers[i] = buffer;
  }

  // Send out xy data.
  for (size_t i = 0; i < dest_procs->size; ++i)
  {
    int proc = source_procs->data[i];
    byte_array_t* buffer = send_buffers[i];
    int err = MPI_Isend(buffer->data, (int)(buffer->size), MPI_BYTE, proc, 0, old_mesh->comm, &requests[source_procs->size + i]);
    if (err != MPI_SUCCESS)
      polymec_error("Error sending xy data from rank %d", proc);
  }

  // Wait for everything to finish.
  MPI_Waitall((int)(source_procs->size + dest_procs->size), requests, MPI_STATUSES_IGNORE);

  // We're finished with the send buffers.
  for (size_t i = 0; i < dest_procs->size; ++i)
    byte_array_free(send_buffers[i]);

  // Unpack the xy data and place it into our own array.
  chunk_xy_data_array_t* all_xy_data = chunk_xy_data_array_new_with_size(old_mesh->num_xy_chunks);
  memset(all_xy_data->data, 0, all_xy_data->size * sizeof(chunk_xy_data_t*));
  size_t offsets[source_procs->size];
  memset(offsets, 0, sizeof(size_t) * source_procs->size);
  for (size_t i = 0; i < old_mesh->num_xy_chunks; ++i)
  {
    if (partition_vector[i] == old_mesh->rank)
    {
      // Fetch the buffer for this process.
      size_t index = int_lower_bound(source_procs->data, source_procs->size, (int)(source_vector[i]));
      byte_array_t* buffer = receive_buffers[index]; 

      // Extract the next xy data thingy from the buffer.
      chunk_xy_data_t* xy_data = serializer_read(ser, buffer, &(offsets[i]));
      chunk_xy_data_array_assign_with_dtor(all_xy_data, i, xy_data, chunk_xy_data_free);
      byte_array_free(buffer);
    }
  }

  // Clean up the rest of the mess.
  polymec_release(ser);
  int_array_free(source_procs);
  int_array_free(dest_procs);

  STOP_FUNCTION_TIMER();
  return all_xy_data;
}

static void redistribute_colmesh(colmesh_t** mesh, 
                                 int64_t* partition,
                                 int64_t* sources)
{
  START_FUNCTION_TIMER();
  colmesh_t* old_mesh = *mesh;

  // Create a new empty mesh.
  colmesh_t* new_mesh = polymec_malloc(sizeof(colmesh_t));
  new_mesh->comm = old_mesh->comm;
  new_mesh->chunks = chunk_map_new();
  new_mesh->chunk_indices = NULL;
  new_mesh->num_xy_chunks = old_mesh->num_xy_chunks;
  new_mesh->num_z_chunks = old_mesh->num_z_chunks;
  new_mesh->nz_per_chunk = old_mesh->nz_per_chunk;
  new_mesh->z1 = old_mesh->z1;
  new_mesh->z2 = old_mesh->z2;
  MPI_Comm_size(new_mesh->comm, &new_mesh->nproc);
  MPI_Comm_rank(new_mesh->comm, &new_mesh->rank);
  new_mesh->finalized = false;

  // Create xy data for chunks.
  new_mesh->chunk_xy_data = redistribute_chunk_xy_data(old_mesh, partition, sources);

  // Insert the new patches as prescribed by the partition vector.
  size_t num_chunks = new_mesh->num_xy_chunks * new_mesh->num_z_chunks;
  for (size_t i = 0; i < num_chunks; ++i)
  {
    if (partition[i] == new_mesh->rank)
    {
      int xy_index = (int)(i / new_mesh->num_z_chunks);
      int z_index = (int)(i % new_mesh->num_z_chunks);
      colmesh_insert_chunk(new_mesh, xy_index, z_index);
    }
  }

  // Replace the old mesh with the new one.
  *mesh = new_mesh;
  STOP_FUNCTION_TIMER();
}

static void redistribute_colmesh_field(colmesh_field_t** field, 
                                       int64_t* partition,
                                       int64_t* sources,
                                       colmesh_t* new_mesh)
{
  START_FUNCTION_TIMER();

  // Create a new field from the old one.
  colmesh_field_t* old_field = *field;
  colmesh_field_t* new_field = colmesh_field_new(new_mesh,
                                                 colmesh_field_centering(old_field),
                                                 colmesh_field_num_components(old_field));

  // Copy all local chunk data from one field to the other.
  colmesh_chunk_data_t* data;
  int pos = 0, xy_index, z_index;
  while (colmesh_field_next_chunk(new_field, &pos, &xy_index, &z_index, &data))
  {
    colmesh_chunk_data_t* old_data = colmesh_field_chunk_data(old_field, xy_index, z_index);
    if (old_data != NULL)
      colmesh_chunk_data_copy(old_data, data);
  }

  // Post receives for each chunk in the new field.
  size_t num_new_local_chunks = colmesh_field_num_chunks(new_field);
  MPI_Request recv_requests[num_new_local_chunks];
  pos = 0;
  size_t num_recv_reqs = 0;
  while (colmesh_field_next_chunk(new_field, &pos, &xy_index, &z_index, &data))
  {
    int ch = chunk_index(new_mesh, xy_index, z_index);
    if (partition[ch] == new_mesh->rank)
    {
      size_t data_size = colmesh_chunk_data_size(data) / sizeof(real_t);
      int err = MPI_Irecv(data->data, (int)data_size, MPI_REAL_T, (int)sources[ch],
                          0, new_mesh->comm, &(recv_requests[num_recv_reqs]));
      if (err != MPI_SUCCESS)
        polymec_error("Error receiving field data from rank %d", (int)sources[ch]);
      ++num_recv_reqs;
    }
  }
  ASSERT(num_recv_reqs <= num_new_local_chunks);

  // Post sends.
  size_t num_old_local_chunks = colmesh_field_num_chunks(old_field);
  MPI_Request send_requests[num_old_local_chunks];
  pos = 0;
  size_t num_send_reqs = 0;
  while (colmesh_field_next_chunk(old_field, &pos, &xy_index, &z_index, &data))
  {
    int ch = chunk_index(new_mesh, xy_index, z_index);
    if (sources[ch] == new_mesh->rank)
    {
      size_t data_size = colmesh_chunk_data_size(data) / sizeof(real_t);
      int err = MPI_Isend(data->data, (int)data_size, MPI_REAL_T, (int)partition[ch],
                          0, new_mesh->comm, &(send_requests[num_send_reqs]));
      if (err != MPI_SUCCESS)
        polymec_error("Error sending field data to rank %d", (int)partition[ch]);
      ++num_send_reqs;
    }
  }
  ASSERT(num_send_reqs <= num_old_local_chunks);

  // Wait for everything to finish.
  MPI_Waitall((int)num_send_reqs, send_requests, MPI_STATUSES_IGNORE);
  MPI_Waitall((int)num_recv_reqs, recv_requests, MPI_STATUSES_IGNORE);

  // Replace the old field with the new one.
  *field = new_field;
  STOP_FUNCTION_TIMER();
}
#endif // POLYMEC_HAVE_MPI

void repartition_colmesh(colmesh_t** mesh, 
                         int* weights,
                         real_t imbalance_tol,
                         colmesh_field_t** fields,
                         size_t num_fields)
{
  ASSERT((weights == NULL) || (imbalance_tol > 0.0));
  ASSERT((weights == NULL) || (imbalance_tol <= 1.0));
  ASSERT((fields != NULL) || (num_fields == 0));
#if POLYMEC_HAVE_MPI

  // On a single process, repartitioning has no meaning.
  colmesh_t* old_mesh = *mesh;
  if (old_mesh->nproc == 1) 
    return;

  START_FUNCTION_TIMER();
  // Map the mesh's graph to the new domains, producing a partition vector.
  // We need the partition vector on all processes in the communicator, so we 
  // scatter it from rank 0.
  log_debug("repartition_colmesh: Repartitioning mesh on %d subdomains.", old_mesh->nproc);
  int64_t* P = partition_graph(old_mesh->chunk_graph, old_mesh->comm, 
                               weights, imbalance_tol, true);

  // Build a sources vector whose ith component is the rank that used to own 
  // the ith patch.
  int64_t* sources = source_vector(old_mesh);

  // Redistribute the mesh. 
  log_debug("repartition_prismmesh: Redistributing mesh.");
  redistribute_colmesh(mesh, P, sources);

  // Redistribute the fields.
  if (num_fields > 0)
    log_debug("repartition_colmesh: Redistributing %d fields.", (int)num_fields);
  for (size_t f = 0; f < num_fields; ++f)
  {
    colmesh_field_t* old_field = fields[f];
    redistribute_colmesh_field(&(fields[f]), P, sources, *mesh);
    colmesh_field_free(old_field);
  }

  // Clean up.
  colmesh_free(old_mesh);
  polymec_free(sources);
  polymec_free(P);

  STOP_FUNCTION_TIMER();
#endif
}

polymesh_t* colmesh_as_polymesh(colmesh_t* mesh)
{
  return NULL; // FIXME
}

// These functions provide access to exchangers for colmesh_fields.
exchanger_t* colmesh_exchanger(colmesh_t* mesh, colmesh_centering_t centering);
exchanger_t* colmesh_exchanger(colmesh_t* mesh, colmesh_centering_t centering)
{
  exchanger_t* ex = NULL;
  switch (centering)
  {
    case COLMESH_CELL: ex = mesh->cell_ex; break;
    case COLMESH_XYFACE: ex = mesh->xy_face_ex; break;
    case COLMESH_ZFACE: ex = mesh->z_face_ex; break;
    case COLMESH_XYEDGE: ex = mesh->xy_edge_ex; break;
    case COLMESH_ZEDGE: ex = mesh->z_edge_ex; break;
    case COLMESH_NODE: ex = mesh->node_ex; break;
  }
  return ex;
}

