// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/array.h"
#include "core/array_utils.h"
#include "core/hilbert.h"
#include "core/partitioning.h"
#include "core/timer.h"
#include "core/unordered_set.h"
#include "geometry/colmesh.h"
#include "geometry/colmesh_field.h"
#include "geometry/polymesh.h"

struct colmesh_fragment_t
{
  planar_polymesh_t* mesh;
  exchanger_proc_map_t* send_map;
  exchanger_proc_map_t* receive_map;
};

static void colmesh_fragment_free(colmesh_fragment_t* fragment)
{
  planar_polymesh_free(fragment->mesh);
  if (fragment->send_map != NULL)
    exchanger_proc_map_free(fragment->send_map);
  if (fragment->receive_map != NULL)
    exchanger_proc_map_free(fragment->receive_map);
  polymec_free(fragment);
}

void colmesh_fragment_map_add(colmesh_fragment_map_t* map,
                              int xy_index,
                              planar_polymesh_t* fragment,
                              exchanger_proc_map_t* send_map,
                              exchanger_proc_map_t* receive_map)
{
  colmesh_fragment_t* frag = polymec_malloc(sizeof(colmesh_fragment_t));
  frag->mesh = fragment;
  frag->send_map = send_map;
  frag->receive_map = receive_map;
  colmesh_fragment_map_insert_with_v_dtor(map, xy_index, frag, colmesh_fragment_free);
}

// Chunk xy data. Shared across all "stacked" chunks.
typedef struct 
{
  int num_columns;
  int num_ghost_columns;
  int* column_xy_face_offsets;
  int* column_xy_faces;
  int* xy_face_columns;
  int num_xy_faces;
  int num_xy_edges;
  int* xy_edge_nodes;
  int num_xy_nodes;
  point2_t* xy_nodes;

  // Exchanger process maps--used to construct exchangers.
  exchanger_proc_map_t* send_map;       // proc -> (send cell, edge) pairs
  exchanger_proc_map_t* receive_map;    // proc -> (receive cell, edge) pairs
} chunk_xy_data_t;

DEFINE_ARRAY(chunk_xy_data_array, chunk_xy_data_t*)

// Creates shared xy chunk data for the local process, generating information
// for the locally-owned cells in the planar polymesh.
static chunk_xy_data_t* chunk_xy_data_from_partition(MPI_Comm comm,
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
      if (col2 == -1) 
      {
        // This column isn't mapped to this chunk, so either col1 is a 
        // boundary column, or col2 is a ghost column (stored on another chunk).
        if (cell2 != -1) // this column is stored on another chunk
          col2 = (int)(xy_data->num_columns + xy_data->receive_map->size); 
      }

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
      if ((col2 >= 0) && (col2 < xy_data->num_columns))
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
        // Set up an exchange mapping between the columns.
        int cell = (cols_reversed) ? cell2 : cell1;
        int neighbor_xy_index = (int)partition_vector[cell];
        int xy_face = xy_data->num_xy_faces;
        exchanger_proc_map_add_index(xy_data->send_map, neighbor_xy_index, col1);
        exchanger_proc_map_add_index(xy_data->send_map, neighbor_xy_index, xy_face);
        exchanger_proc_map_add_index(xy_data->receive_map, neighbor_xy_index, col2);
        exchanger_proc_map_add_index(xy_data->receive_map, neighbor_xy_index, xy_face);
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
  xy_data->num_xy_edges = (int)(edge_nodes->size/2);
  ASSERT(xy_data->num_xy_edges == xy_data->num_xy_faces);
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

// Creates shared xy chunk data for the local process from the given fragment.
static chunk_xy_data_t* chunk_xy_data_from_fragment(colmesh_fragment_t* fragment)
{
  planar_polymesh_t* mesh = fragment->mesh;
  // Copy the locally stored over from the planar polymesh.
  chunk_xy_data_t* xy_data = polymec_malloc(sizeof(chunk_xy_data_t));
  xy_data->num_columns = mesh->num_cells;
  xy_data->num_ghost_columns = 0;
  xy_data->num_xy_faces = mesh->num_edges;
  xy_data->num_xy_edges = mesh->num_edges;
  xy_data->num_xy_nodes = mesh->num_nodes;
  xy_data->column_xy_face_offsets = polymec_malloc(sizeof(int) * (xy_data->num_columns+1));
  memcpy(xy_data->column_xy_face_offsets, mesh->cell_edge_offsets, sizeof(int) * (xy_data->num_columns+1));
  xy_data->column_xy_faces = polymec_malloc(sizeof(int) * xy_data->column_xy_face_offsets[xy_data->num_columns]);
  memcpy(xy_data->column_xy_faces, mesh->cell_edges, sizeof(int) * xy_data->column_xy_face_offsets[xy_data->num_columns]);
  xy_data->xy_face_columns = polymec_malloc(sizeof(int) * 2 * xy_data->num_xy_faces);
  memcpy(xy_data->xy_face_columns, mesh->edge_cells, sizeof(int) * 2 * xy_data->num_xy_faces);
  xy_data->xy_edge_nodes = polymec_malloc(sizeof(int) * 2 * xy_data->num_xy_edges);
  memcpy(xy_data->xy_edge_nodes, mesh->edge_nodes, sizeof(int) * 2 * xy_data->num_xy_edges);
  xy_data->xy_nodes = polymec_malloc(sizeof(point2_t) * xy_data->num_xy_nodes);
  memcpy(xy_data->xy_nodes, mesh->nodes, sizeof(point2_t) * xy_data->num_xy_nodes);

  // Steal the send and receive maps.
  xy_data->send_map = fragment->send_map;
  fragment->send_map = NULL;
  xy_data->receive_map = fragment->receive_map;
  fragment->receive_map = NULL;

  // Count up ghost cells.
  int pos = 0, proc;
  int_array_t* indices;
  while (exchanger_proc_map_next(xy_data->receive_map, &pos, &proc, &indices))
  {
    ASSERT((indices->size % 2) == 0);
    xy_data->num_ghost_columns += indices->size/2;
  }

  // That's it!
  return xy_data;
}

#if POLYMEC_HAVE_MPI
static int_array_t* clone_int_array(int_array_t* array)
{
  return int_array_clone(array, NULL);
}

static chunk_xy_data_t* chunk_xy_data_clone(chunk_xy_data_t* xy_data)
{
  chunk_xy_data_t* clone = polymec_malloc(sizeof(chunk_xy_data_t));
  clone->num_columns = xy_data->num_columns;
  clone->num_ghost_columns = xy_data->num_ghost_columns;
  clone->column_xy_face_offsets = polymec_malloc((xy_data->num_columns+1) * sizeof(int));
  memcpy(clone->column_xy_face_offsets, xy_data->column_xy_face_offsets, (xy_data->num_columns+1) * sizeof(int));
  clone->column_xy_faces = polymec_malloc(xy_data->column_xy_face_offsets[xy_data->num_columns] * sizeof(int));
  memcpy(clone->column_xy_faces, xy_data->column_xy_faces, xy_data->column_xy_face_offsets[xy_data->num_columns] * sizeof(int));
  clone->num_xy_faces = xy_data->num_xy_faces;
  clone->xy_face_columns = polymec_malloc(2 * xy_data->num_xy_faces * sizeof(int));
  memcpy(clone->xy_face_columns, xy_data->xy_face_columns, 2 * xy_data->num_xy_faces * sizeof(int));
  clone->num_xy_edges = xy_data->num_xy_edges;
  clone->xy_edge_nodes = polymec_malloc(2 * xy_data->num_xy_edges * sizeof(int));
  memcpy(clone->xy_edge_nodes, xy_data->xy_edge_nodes, 2 * xy_data->num_xy_edges * sizeof(int));
  clone->num_xy_nodes = xy_data->num_xy_nodes;
  clone->xy_nodes = polymec_malloc(xy_data->num_xy_nodes * sizeof(point2_t));
  memcpy(clone->xy_nodes, xy_data->xy_nodes, xy_data->num_xy_nodes * sizeof(point2_t));
  clone->send_map = exchanger_proc_map_clone(xy_data->send_map, NULL, clone_int_array, NULL, int_array_free);
  clone->receive_map = exchanger_proc_map_clone(xy_data->send_map, NULL, clone_int_array, NULL, int_array_free);
  return clone;
}
#endif

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

  int num_xy_chunks, num_z_chunks, nz_per_chunk;
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

static void colmesh_create_chunk_graph(colmesh_t* mesh)
{
  ASSERT(mesh->chunk_graph == NULL);
  int num_chunks = mesh->num_xy_chunks * mesh->num_z_chunks;
  mesh->chunk_graph = adj_graph_new(MPI_COMM_SELF, num_chunks);
  for (int xy_index = 0; xy_index < mesh->num_xy_chunks; ++xy_index)
  {
    if (xy_index < mesh->chunk_xy_data->size)
    {
      chunk_xy_data_t* xy_data = mesh->chunk_xy_data->data[xy_index];
      if (xy_data == NULL) continue;

      // Extract the neighboring chunks from this xy data and add edges to our graph.
      for (int z_index = 0; z_index < mesh->num_z_chunks; ++z_index)
      {
        int ch_index = chunk_index(mesh, xy_index, z_index);

        // Count up edges.
        int num_xy_neighbors = (int)xy_data->send_map->size;
        int num_z_neighbors;
        if (mesh->num_z_chunks == 1)
          num_z_neighbors = 0;
        else if ((z_index == 0) || (z_index == (mesh->num_z_chunks-1)))
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
        if (z_index < (mesh->num_z_chunks-1))
          edges[i++] = chunk_index(mesh, xy_index, z_index+1);

        ASSERT(i == num_xy_neighbors + num_z_neighbors);
      }
    }
  }
}

colmesh_t* create_empty_colmesh(MPI_Comm comm, 
                                planar_polymesh_t* columns,
                                real_t z1, real_t z2,
                                int num_xy_chunks, int num_z_chunks,
                                int nz_per_chunk, bool periodic_in_z)
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
  mesh->chunk_graph = NULL;
  mesh->num_xy_chunks = num_xy_chunks;
  mesh->num_z_chunks = num_z_chunks;
  mesh->nz_per_chunk = nz_per_chunk;
  mesh->z1 = z1;
  mesh->z2 = z2;
  MPI_Comm_size(comm, &mesh->nproc);
  MPI_Comm_rank(comm, &mesh->rank);
  mesh->periodic_in_z = periodic_in_z;
  mesh->cell_ex = NULL;
  mesh->xy_face_ex = NULL;
  mesh->z_face_ex = NULL;
  mesh->xy_edge_ex = NULL;
  mesh->z_edge_ex = NULL;
  mesh->node_ex = NULL;
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
  release_ref(poly);
#endif

  // Create xy data for chunks. 
  mesh->chunk_xy_data = chunk_xy_data_array_new();
  for (int xy_index = 0; xy_index < (int)num_xy_chunks; ++xy_index)
  {
    chunk_xy_data_t* xy_data = chunk_xy_data_from_partition(comm, columns, P, xy_index);
    chunk_xy_data_array_append_with_dtor(mesh->chunk_xy_data, xy_data, 
                                         chunk_xy_data_free);
  }
  polymec_free(P);

  // Create a global graph whose vertices are the mesh's chunks.
  colmesh_create_chunk_graph(mesh);

  return mesh;
}

colmesh_t* create_empty_colmesh_from_fragments(MPI_Comm comm,
                                               colmesh_fragment_map_t* local_fragments,
                                               real_t z1, real_t z2,
                                               int num_xy_chunks, 
                                               int num_z_chunks, 
                                               int nz_per_chunk, 
                                               bool periodic_in_z)
{
  ASSERT(local_fragments != NULL);
  ASSERT(local_fragments->size > 0);
  ASSERT(z1 < z2);
  ASSERT(num_xy_chunks > 0);
  ASSERT(num_z_chunks > 0);
  ASSERT(nz_per_chunk > 0);

  colmesh_t* mesh = polymec_malloc(sizeof(colmesh_t));
  mesh->comm = comm;
  mesh->chunks = chunk_map_new();
  mesh->chunk_indices = NULL;
  mesh->chunk_graph = NULL;
  mesh->num_xy_chunks = num_xy_chunks;
  mesh->num_z_chunks = num_z_chunks;
  mesh->nz_per_chunk = nz_per_chunk;
  mesh->z1 = z1;
  mesh->z2 = z2;
  MPI_Comm_size(comm, &mesh->nproc);
  MPI_Comm_rank(comm, &mesh->rank);
  mesh->periodic_in_z = periodic_in_z;
  mesh->cell_ex = NULL;
  mesh->xy_face_ex = NULL;
  mesh->z_face_ex = NULL;
  mesh->xy_edge_ex = NULL;
  mesh->z_edge_ex = NULL;
  mesh->node_ex = NULL;
  mesh->finalized = false;

  // Create xy data from the distributed fragments.
  mesh->chunk_xy_data = chunk_xy_data_array_new();
  int pos = 0, xy;
  colmesh_fragment_t* fragment;
  while (colmesh_fragment_map_next(local_fragments, &pos, &xy, &fragment))
  {
    chunk_xy_data_t* xy_data = chunk_xy_data_from_fragment(fragment);
    chunk_xy_data_array_resize(mesh->chunk_xy_data, MAX(mesh->chunk_xy_data->size, xy+1));
    chunk_xy_data_array_assign_with_dtor(mesh->chunk_xy_data, xy, 
                                         xy_data, chunk_xy_data_free);
  }

  // Destroy the fragment map.
  colmesh_fragment_map_free(local_fragments);

  // Build the chunk graph.
  colmesh_create_chunk_graph(mesh);

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
  chunk->xy_edge_nodes = xy_data->xy_edge_nodes;
  chunk->num_xy_nodes = xy_data->num_xy_nodes;
  chunk->xy_nodes = xy_data->xy_nodes;

  // Chunk z data.
  int nz = mesh->nz_per_chunk * mesh->num_z_chunks;
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

  // If you trip this assertion, you have some missing chunks in your mesh.
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

DEFINE_UNORDERED_MAP(proc_point_map, int, point_array_t*, int_hash, int_equals)
static void proc_point_map_add(proc_point_map_t* map, 
                               int process, point_t* x)
{
  point_array_t** points_p = proc_point_map_get(map, process);
  point_array_t* points = NULL;
  if (points_p != NULL)
    points = *points_p;
  else
  {
    points = point_array_new();
    proc_point_map_insert_with_v_dtor(map, process, points, point_array_free);
  }
  point_array_append(points, *x);
}

static void find_chunk_offsets(colmesh_t* mesh, 
                               colmesh_centering_t centering,
                               int* chunk_offsets)
{
  chunk_offsets[0] = 0;
  int k = 0;
  for (size_t i = 0; i < mesh->chunks->size; ++i)
  {
    int xy = mesh->chunk_indices[2*i];
    int z = mesh->chunk_indices[2*i+1];
    int index = chunk_index(mesh, xy, z);
    colmesh_chunk_t* chunk = *chunk_map_get(mesh->chunks, index);
    int nxy, nz;
    switch (centering)
    {
      case COLMESH_CELL:   
        nxy = chunk->num_columns + chunk->num_ghost_columns;
        nz  = chunk->num_z_cells + 2; 
        break;
      case COLMESH_XYFACE: 
        nxy = chunk->num_xy_faces;
        nz = chunk->num_z_cells; 
        break;
      case COLMESH_ZFACE:  
        nxy = chunk->num_columns;
        nz = chunk->num_z_cells + 1; 
        break;
      case COLMESH_XYEDGE: 
        nxy = chunk->num_xy_faces;
        nz = chunk->num_z_cells + 1; 
        break;
      case COLMESH_ZEDGE:  
        nxy = chunk->num_xy_nodes;
        nz = chunk->num_z_cells; 
        break;
      case COLMESH_NODE:   
        nxy = chunk->num_xy_nodes;
        nz = chunk->num_z_cells + 1; 
    }
    chunk_offsets[k+1] = chunk_offsets[k] + nxy * nz;
    ++k;
  }
}

// This function sorts the indices in the given exchanger_proc_map according to 
// their spatial positions as reflected in the point map.
static void sort_indices(proc_point_map_t* point_map,
                         exchanger_proc_map_t* index_map)
{
  int pos = 0, proc;
  int_array_t* ind;
  while (exchanger_proc_map_next(index_map, &pos, &proc, &ind))
  {
    point_array_t* pts = *proc_point_map_get(point_map, proc);

    // Make sure we have the same number of points and indices!
    ASSERT(pts->size == ind->size);

    bbox_t bbox = {.x1 = 0.0, .x2 = 0.0, .y1 = 0.0, .y2 = 0.0, .z1 = 0.0, .z2 = 0.0};
    for (size_t k = 0; k < ind->size; ++k)
      bbox_grow(&bbox, &pts->data[k]);
    hilbert_t* curve = hilbert_new(&bbox);
    hilbert_sort_points(curve, pts->data, ind->data, ind->size);
    release_ref(curve);
  }
}

static void create_cell_ex(colmesh_t* mesh)
{
  START_FUNCTION_TIMER();
  ASSERT(mesh->cell_ex == NULL);

  // Figure out which processes own what chunks.
#if POLYMEC_HAVE_MPI
  int64_t* owners = source_vector(mesh);
#else
  int num_all_chunks = (int)(mesh->num_xy_chunks * mesh->num_z_chunks);
  int64_t* owners = polymec_calloc(num_all_chunks, sizeof(int64_t));
#endif

  // Determine offsets for chunks in the list of chunks for this mesh.
  int chunk_offsets[mesh->chunks->size+1];
  find_chunk_offsets(mesh, COLMESH_CELL, chunk_offsets);

  // Assemble exchangers by traversing the locally stored chunks, assembling 
  // send and receive indices, and ordering those indices by the spatial locations
  // of their underlying elements. For cell exchangers, we sort the indices by the 
  // location of the faces separating send and receive cells. 

  // Set up send cells.
  exchanger_proc_map_t* send_map = exchanger_proc_map_new();
  proc_point_map_t* point_map = proc_point_map_new();
  for (size_t i = 0; i < mesh->chunks->size; ++i)
  {
    int xy = mesh->chunk_indices[2*i];
    int z = mesh->chunk_indices[2*i+1];
    int ch_index = chunk_index(mesh, xy, z);
    colmesh_chunk_t* chunk = *chunk_map_get(mesh->chunks, ch_index);
    int chunk_offset = chunk_offsets[i];
    real_t dz = (chunk->z2 - chunk->z1) / chunk->num_z_cells;

    chunk_xy_data_t* xy_data = mesh->chunk_xy_data->data[xy];
    exchanger_proc_map_t* xy_sends = xy_data->send_map;

    // Figure out the xy send cells.
    int pos = 0, neighbor_xy_index;
    int_array_t* indices;
    while (exchanger_proc_map_next(xy_sends, &pos, &neighbor_xy_index, &indices))
    {
      // Loop over the xy indices in the map.
      size_t num_indices = indices->size/2;
      for (size_t j = 0; j < num_indices; ++j)
      {
        int proc = (int)(owners[neighbor_xy_index]);
        int xy1 = indices->data[2*j];
        int edge = indices->data[2*j+1];

        // Get the x and y coordinates for the send cell's face.
        int node_indices[2] = {xy_data->xy_edge_nodes[2*edge], xy_data->xy_edge_nodes[2*edge+1]};
        point2_t nodes[2] = {xy_data->xy_nodes[node_indices[0]], xy_data->xy_nodes[node_indices[1]]};
        point_t x = {.x = 0.5 * (nodes[0].x + nodes[1].x), 
                     .y = 0.5 * (nodes[0].y + nodes[1].y), 
                     .z = 0.0}; 

        // Traverse the column and add send indices/points.
        for (int zz = 1; zz <= chunk->num_z_cells; ++zz)
        {
          int index = (int)(chunk_offset + (chunk->num_z_cells+2)*xy1 + zz);
          exchanger_proc_map_add_index(send_map, proc, index);
          x.z = chunk->z1 + ((zz-1)+0.5) * dz;
          proc_point_map_add(point_map, proc, &x);
        }
      }
    }
  }

  // Sort the xy send indices/points.
  sort_indices(point_map, send_map);

  // Now set up receive cells.
  exchanger_proc_map_t* receive_map = exchanger_proc_map_new();
  proc_point_map_clear(point_map);
  for (size_t i = 0; i < mesh->chunks->size; ++i)
  {
    int xy = mesh->chunk_indices[2*i];
    int z = mesh->chunk_indices[2*i+1];
    int ch_index = chunk_index(mesh, xy, z);
    colmesh_chunk_t* chunk = *chunk_map_get(mesh->chunks, ch_index);
    int chunk_offset = chunk_offsets[i];
    real_t dz = (chunk->z2 - chunk->z1) / chunk->num_z_cells;

    chunk_xy_data_t* xy_data = mesh->chunk_xy_data->data[xy];
    exchanger_proc_map_t* xy_receives = xy_data->receive_map;

    int pos = 0, neighbor_xy_index;
    int_array_t* indices;
    int cell_offset = (int)((chunk->num_z_cells+2) * chunk->num_columns);
    while (exchanger_proc_map_next(xy_receives, &pos, &neighbor_xy_index, &indices))
    {
      size_t num_indices = indices->size/2;
      for (size_t j = 0; j < num_indices; ++j)
      {
        int proc = (int)(owners[neighbor_xy_index]);
        int edge = indices->data[2*j+1];

        // Get the x and y coordinates for the receive cell's face.
        int node_indices[2] = {xy_data->xy_edge_nodes[2*edge], xy_data->xy_edge_nodes[2*edge+1]};
        point2_t nodes[2] = {xy_data->xy_nodes[node_indices[0]], xy_data->xy_nodes[node_indices[1]]};
        point_t x = {.x = 0.5 * (nodes[0].x + nodes[1].x), 
                     .y = 0.5 * (nodes[0].y + nodes[1].y), 
                     .z = 0.0}; 

        // Traverse the column and add receive indices/points.
        for (int zz = 0; zz < chunk->num_z_cells; ++zz)
        {
          int index = (int)chunk_offset + cell_offset;
          exchanger_proc_map_add_index(receive_map, proc, index);
          x.z = chunk->z1 + ((zz-1)+0.5) * dz;
          proc_point_map_add(point_map, proc, &x);
          ++cell_offset;
        }
      }
      ASSERT(cell_offset < chunk_offset + (chunk->num_z_cells+2) * (chunk->num_columns + chunk->num_ghost_columns));
    }
  }

  // Sort the xy receive indices/points.
  sort_indices(point_map, receive_map);
  proc_point_map_free(point_map);

  // Now hook up the z sends/receives.
  for (size_t i = 0; i < mesh->chunks->size; ++i)
  {
    int xy = mesh->chunk_indices[2*i];
    int z = mesh->chunk_indices[2*i+1];
    int ch_index = chunk_index(mesh, xy, z);
    colmesh_chunk_t* chunk = *chunk_map_get(mesh->chunks, ch_index);
    int chunk_offset = chunk_offsets[i];

    // Upper neighbor.
    if ((z > 0) || (mesh->periodic_in_z))
    {
      for (int xy1 = 0; xy1 < chunk->num_columns; ++xy1)
      {
        int ch1_index = ((z == 0) && mesh->periodic_in_z) ? 
                        chunk_index(mesh, xy, mesh->num_z_chunks-1) : chunk_index(mesh, xy, z-1);
        int proc = (int)(owners[ch1_index]);
        int z1 = 1;
        int send_index = (int)(chunk_offset + (chunk->num_z_cells+2) * xy1 + z1);
        int receive_index = (int)(chunk_offset + (chunk->num_z_cells+2) * xy1 + z1 - 1);
        exchanger_proc_map_add_index(send_map, proc, send_index);
        exchanger_proc_map_add_index(receive_map, proc, receive_index);
      }
    }

    // Lower neighbor.
    if ((z < (mesh->num_z_chunks-1)) || (mesh->periodic_in_z))
    {
      for (int xy1 = 0; xy1 < chunk->num_columns; ++xy1)
      {
        int ch1_index = ((z == (mesh->num_z_chunks-1)) && mesh->periodic_in_z) ? 
                        chunk_index(mesh, xy, 0) : chunk_index(mesh, xy, z+1);
        int proc = (int)(owners[ch1_index]);
        int z1 = chunk->num_z_cells;
        int send_index = (int)(chunk_offset + (chunk->num_z_cells+2) * xy1 + z1);
        int receive_index = (int)(chunk_offset + (chunk->num_z_cells+2) * xy1 + z1 + 1);
        exchanger_proc_map_add_index(send_map, proc, send_index);
        exchanger_proc_map_add_index(receive_map, proc, receive_index);
      }
    }
  }
  polymec_free(owners);

  // Now construct the exchanger.
  mesh->cell_ex = exchanger_new(mesh->comm);
  exchanger_set_sends(mesh->cell_ex, send_map);
  exchanger_set_receives(mesh->cell_ex, receive_map);
  STOP_FUNCTION_TIMER();
}

static void create_xy_face_ex(colmesh_t* mesh)
{
  START_FUNCTION_TIMER();

  // Figure out which processes own what chunks.
#if POLYMEC_HAVE_MPI
  int64_t* owners = source_vector(mesh);
#else
  int num_all_chunks = (int)(mesh->num_xy_chunks * mesh->num_z_chunks);
  int64_t* owners = polymec_calloc(num_all_chunks, sizeof(int64_t));
#endif

  // Determine offsets for chunks in the list of chunks for this mesh.
  int chunk_offsets[mesh->chunks->size+1];
  find_chunk_offsets(mesh, COLMESH_XYFACE, chunk_offsets);

  // Assemble exchangers by traversing the locally stored chunks, assembling 
  // send and receive indices, and ordering those indices by the spatial locations
  // of their underlying elements. 

  // Set up send faces.
  exchanger_proc_map_t* send_map = exchanger_proc_map_new();
  proc_point_map_t* point_map = proc_point_map_new();
  int_unordered_set_t* contributed_to_self = int_unordered_set_new();
  for (size_t i = 0; i < mesh->chunks->size; ++i)
  {
    int xy = mesh->chunk_indices[2*i];
    int z = mesh->chunk_indices[2*i+1];
    int ch_index = chunk_index(mesh, xy, z);
    colmesh_chunk_t* chunk = *chunk_map_get(mesh->chunks, ch_index);
    int chunk_offset = chunk_offsets[i];
    real_t dz = (chunk->z2 - chunk->z1) / chunk->num_z_cells;

    chunk_xy_data_t* xy_data = mesh->chunk_xy_data->data[xy];
    exchanger_proc_map_t* xy_sends = xy_data->send_map;

    // Figure out the send faces.
    int pos = 0, neighbor_xy_index;
    int_array_t* indices;
    int_unordered_set_clear(contributed_to_self);
    while (exchanger_proc_map_next(xy_sends, &pos, &neighbor_xy_index, &indices))
    {
      // Loop over the xy indices in the map.
      size_t num_indices = indices->size/2;
      for (size_t j = 0; j < num_indices; ++j)
      {
        int proc = (int)(owners[neighbor_xy_index]);
        int face = indices->data[2*j+1], edge = face;

        // Get the x and y coordinates for the send face.
        int node_indices[2] = {xy_data->xy_edge_nodes[2*edge], xy_data->xy_edge_nodes[2*edge+1]};
        point2_t nodes[2] = {xy_data->xy_nodes[node_indices[0]], xy_data->xy_nodes[node_indices[1]]};
        point_t x = {.x = 0.5 * (nodes[0].x + nodes[1].x), 
                     .y = 0.5 * (nodes[0].y + nodes[1].y), 
                     .z = 0.0}; 

        // Traverse the column and add send indices/points.
        for (int zz = 0; zz < chunk->num_z_cells; ++zz)
        {
          int index = (int)(chunk_offset + chunk->num_z_cells*face + zz);
          x.z = chunk->z1 + (zz+0.5) * dz;

          // Self contribution.
          if (!int_unordered_set_contains(contributed_to_self, index))
          {
            exchanger_proc_map_add_index(send_map, mesh->rank, index);
            proc_point_map_add(point_map, mesh->rank, &x);
            int_unordered_set_insert(contributed_to_self, index);
          }

          // Contribution to neighboring process.
          if (proc != mesh->rank)
          {
            exchanger_proc_map_add_index(send_map, proc, index);
            proc_point_map_add(point_map, proc, &x);
          }
        }
      }
    }
  }

  // Sort the send faces.
  sort_indices(point_map, send_map);

  // Set up receive faces.
  exchanger_proc_map_t* receive_map = exchanger_proc_map_new();
  proc_point_map_clear(point_map);
  for (size_t i = 0; i < mesh->chunks->size; ++i)
  {
    int xy = mesh->chunk_indices[2*i];
    int z = mesh->chunk_indices[2*i+1];
    int ch_index = chunk_index(mesh, xy, z);
    colmesh_chunk_t* chunk = *chunk_map_get(mesh->chunks, ch_index);
    int chunk_offset = chunk_offsets[i];
    real_t dz = (chunk->z2 - chunk->z1) / chunk->num_z_cells;

    chunk_xy_data_t* xy_data = mesh->chunk_xy_data->data[xy];
    exchanger_proc_map_t* xy_receives = xy_data->receive_map;

    // Figure out the xy receive faces.
    int pos = 0, neighbor_xy_index;
    int_array_t* indices;
    int_unordered_set_clear(contributed_to_self);
    while (exchanger_proc_map_next(xy_receives, &pos, &neighbor_xy_index, &indices))
    {
      size_t num_indices = indices->size/2;
      for (size_t j = 0; j < num_indices; ++j)
      {
        int proc = (int)(owners[neighbor_xy_index]);
        int face = indices->data[2*j+1], edge = face;

        // Get the x and y coordinates for the receive face.
        int node_indices[2] = {xy_data->xy_edge_nodes[2*edge], xy_data->xy_edge_nodes[2*edge+1]};
        point2_t nodes[2] = {xy_data->xy_nodes[node_indices[0]], xy_data->xy_nodes[node_indices[1]]};
        point_t x = {.x = 0.5 * (nodes[0].x + nodes[1].x), 
                     .y = 0.5 * (nodes[0].y + nodes[1].y), 
                     .z = 0.0}; 

        // Traverse the column and add receive indices/points.
        for (int zz = 0; zz < chunk->num_z_cells; ++zz)
        {
          int index = (int)(chunk_offset + chunk->num_z_cells*face + zz);
          x.z = chunk->z1 + (zz+0.5) * dz;

          // Self contribution.
          if (!int_unordered_set_contains(contributed_to_self, index))
          {
            exchanger_proc_map_add_index(receive_map, mesh->rank, index);
            proc_point_map_add(point_map, mesh->rank, &x);
            int_unordered_set_insert(contributed_to_self, index);
          }

          // Contribution to neighboring process.
          if (proc != mesh->rank)
          {
            proc_point_map_add(point_map, proc, &x);
            exchanger_proc_map_add_index(receive_map, proc, index);
          }
        }
      }
    }
  }
  int_unordered_set_free(contributed_to_self);
  polymec_free(owners);

  // Sort the receive faces.
  sort_indices(point_map, receive_map);
  proc_point_map_free(point_map);

  // Now construct the exchanger.
  mesh->xy_face_ex = exchanger_new(mesh->comm);
  exchanger_set_sends(mesh->xy_face_ex, send_map);
  exchanger_set_receives(mesh->xy_face_ex, receive_map);

  // By default, this exchanger uses the "min rank" reducer.
  exchanger_set_reducer(mesh->xy_face_ex, EXCHANGER_MIN_RANK);
  STOP_FUNCTION_TIMER();
}

static void create_z_face_ex(colmesh_t* mesh)
{
  START_FUNCTION_TIMER();

  // Figure out which processes own what chunks.
#if POLYMEC_HAVE_MPI
  int64_t* owners = source_vector(mesh);
#else
  int num_all_chunks = (int)(mesh->num_xy_chunks * mesh->num_z_chunks);
  int64_t* owners = polymec_calloc(num_all_chunks, sizeof(int64_t));
#endif

  // Determine offsets for chunks in the list of chunks for this mesh.
  int chunk_offsets[mesh->chunks->size+1];
  find_chunk_offsets(mesh, COLMESH_ZFACE, chunk_offsets);

  // Hook up the z send/receive faces.
  exchanger_proc_map_t* send_map = exchanger_proc_map_new();
  exchanger_proc_map_t* receive_map = exchanger_proc_map_new();
  for (size_t i = 0; i < mesh->chunks->size; ++i)
  {
    int xy = mesh->chunk_indices[2*i];
    int z = mesh->chunk_indices[2*i+1];
    int ch_index = chunk_index(mesh, xy, z);
    colmesh_chunk_t* chunk = *chunk_map_get(mesh->chunks, ch_index);
    int chunk_offset = chunk_offsets[i];

    // Upper neighbor.
    if ((z > 0) || (mesh->periodic_in_z))
    {
      for (int xy1 = 0; xy1 < chunk->num_columns; ++xy1)
      {
        int ch1_index = ((z == 0) && mesh->periodic_in_z) ? 
                        chunk_index(mesh, xy, mesh->num_z_chunks-1) : chunk_index(mesh, xy, z-1);
        int proc = (int)(owners[ch1_index]);
        int z1 = 0;
        int index = (int)(chunk_offset + (chunk->num_z_cells+1) * xy1 + z1);

        // Self contribution.
        exchanger_proc_map_add_index(send_map, mesh->rank, index);
        exchanger_proc_map_add_index(receive_map, mesh->rank, index);

        // Contribution to neighboring process.
        if (proc != mesh->rank)
        {
          exchanger_proc_map_add_index(send_map, proc, index);
          exchanger_proc_map_add_index(receive_map, proc, index);
        }
      }
    }

    // Lower neighbor.
    if ((z < (mesh->num_z_chunks-1)) || (mesh->periodic_in_z))
    {
      for (int xy1 = 0; xy1 < chunk->num_columns; ++xy1)
      {
        int ch1_index = ((z == (mesh->num_z_chunks-1)) && mesh->periodic_in_z) ? 
                        chunk_index(mesh, xy, 0) : chunk_index(mesh, xy, z+1);
        int proc = (int)(owners[ch1_index]);
        int z1 = chunk->num_z_cells;
        int index = (int)(chunk_offset + (chunk->num_z_cells+1) * xy1 + z1);
        exchanger_proc_map_add_index(send_map, mesh->rank, index);
        exchanger_proc_map_add_index(receive_map, mesh->rank, index);
        if (proc != mesh->rank)
        {
          exchanger_proc_map_add_index(send_map, proc, index);
          exchanger_proc_map_add_index(receive_map, proc, index);
        }
      }
    }
  }
  polymec_free(owners);

  // Now construct the exchanger.
  mesh->z_face_ex = exchanger_new(mesh->comm);
  exchanger_set_sends(mesh->z_face_ex, send_map);
  exchanger_set_receives(mesh->z_face_ex, receive_map);

  // By default, this exchanger uses the "min rank" reducer.
  exchanger_set_reducer(mesh->z_face_ex, EXCHANGER_MIN_RANK);
  STOP_FUNCTION_TIMER();
}

static void create_xy_edge_ex(colmesh_t* mesh)
{
  START_FUNCTION_TIMER();

  // Figure out which processes own what chunks.
#if POLYMEC_HAVE_MPI
  int64_t* owners = source_vector(mesh);
#else
  int num_all_chunks = (int)(mesh->num_xy_chunks * mesh->num_z_chunks);
  int64_t* owners = polymec_calloc(num_all_chunks, sizeof(int64_t));
#endif

  // Determine offsets for chunks in the list of chunks for this mesh.
  int chunk_offsets[mesh->chunks->size+1];
  find_chunk_offsets(mesh, COLMESH_XYEDGE, chunk_offsets);

  // Assemble exchangers by traversing the locally stored chunks, assembling 
  // send and receive indices, and ordering those indices by the spatial locations
  // of their underlying elements. 

  // Set up send edges.
  exchanger_proc_map_t* send_map = exchanger_proc_map_new();
  proc_point_map_t* point_map = proc_point_map_new();
  int_unordered_set_t* contributed_to_self = int_unordered_set_new();
  for (size_t i = 0; i < mesh->chunks->size; ++i)
  {
    int xy = mesh->chunk_indices[2*i];
    int z = mesh->chunk_indices[2*i+1];
    int ch_index = chunk_index(mesh, xy, z);
    colmesh_chunk_t* chunk = *chunk_map_get(mesh->chunks, ch_index);
    int chunk_offset = chunk_offsets[i];
    real_t dz = (chunk->z2 - chunk->z1) / chunk->num_z_cells;

    chunk_xy_data_t* xy_data = mesh->chunk_xy_data->data[xy];
    exchanger_proc_map_t* xy_sends = xy_data->send_map;

    // Figure out the send edges.
    int pos = 0, neighbor_xy_index;
    int_array_t* indices;
    int_unordered_set_clear(contributed_to_self);
    while (exchanger_proc_map_next(xy_sends, &pos, &neighbor_xy_index, &indices))
    {
      // Loop over the xy indices in the map.
      size_t num_indices = indices->size/2;
      for (size_t j = 0; j < num_indices; ++j)
      {
        int proc = (int)(owners[neighbor_xy_index]);
        int edge = indices->data[2*j+1];

        // Get the x and y coordinates for the send edge.
        int node_indices[2] = {xy_data->xy_edge_nodes[2*edge], xy_data->xy_edge_nodes[2*edge+1]};
        point2_t nodes[2] = {xy_data->xy_nodes[node_indices[0]], xy_data->xy_nodes[node_indices[1]]};
        point_t x = {.x = 0.5 * (nodes[0].x + nodes[1].x), 
                     .y = 0.5 * (nodes[0].y + nodes[1].y), 
                     .z = 0.0}; 

        // Traverse the column and add send indices/points.
        for (int zz = 0; zz <= chunk->num_z_cells; ++zz)
        {
          int index = (int)(chunk_offset + (chunk->num_z_cells+1)*edge + zz);
          x.z = chunk->z1 + zz * dz;

          // Self contribution.
          if (!int_unordered_set_contains(contributed_to_self, index))
          {
            exchanger_proc_map_add_index(send_map, mesh->rank, index);
            proc_point_map_add(point_map, mesh->rank, &x);
          }

          // Contribution to neighboring process.
          if (proc != mesh->rank)
          {
            exchanger_proc_map_add_index(send_map, proc, index);
            proc_point_map_add(point_map, proc, &x);
          }
        }
      }
    }
  }

  // Sort the send edges.
  sort_indices(point_map, send_map);

  // Set up receive edges.
  exchanger_proc_map_t* receive_map = exchanger_proc_map_new();
  proc_point_map_clear(point_map);
  for (size_t i = 0; i < mesh->chunks->size; ++i)
  {
    int xy = mesh->chunk_indices[2*i];
    int z = mesh->chunk_indices[2*i+1];
    int ch_index = chunk_index(mesh, xy, z);
    colmesh_chunk_t* chunk = *chunk_map_get(mesh->chunks, ch_index);
    int chunk_offset = chunk_offsets[i];
    real_t dz = (chunk->z2 - chunk->z1) / chunk->num_z_cells;

    chunk_xy_data_t* xy_data = mesh->chunk_xy_data->data[xy];
    exchanger_proc_map_t* xy_receives = xy_data->receive_map;

    // Figure out the xy receive edges.
    int pos = 0, neighbor_xy_index;
    int_array_t* indices;
    int_unordered_set_clear(contributed_to_self);
    while (exchanger_proc_map_next(xy_receives, &pos, &neighbor_xy_index, &indices))
    {
      size_t num_indices = indices->size/2;
      for (size_t j = 0; j < num_indices; ++j)
      {
        int proc = (int)(owners[neighbor_xy_index]);
        int edge = indices->data[2*j+1];

        // Get the x and y coordinates for the receive cell's face.
        int node_indices[2] = {xy_data->xy_edge_nodes[2*edge], xy_data->xy_edge_nodes[2*edge+1]};
        point2_t nodes[2] = {xy_data->xy_nodes[node_indices[0]], xy_data->xy_nodes[node_indices[1]]};
        point_t x = {.x = 0.5 * (nodes[0].x + nodes[1].x), 
                     .y = 0.5 * (nodes[0].y + nodes[1].y), 
                     .z = 0.0}; 

        // Traverse the column and add receive indices/points.
        for (int zz = 0; zz <= chunk->num_z_cells; ++zz)
        {
          int index = (int)(chunk_offset + (chunk->num_z_cells+1)*edge + zz);
          x.z = chunk->z1 + zz * dz;

          // Self contribution.
          if (!int_unordered_set_contains(contributed_to_self, index))
          {
            exchanger_proc_map_add_index(receive_map, mesh->rank, index);
            proc_point_map_add(point_map, mesh->rank, &x);
          }

          // Contribution to neighboring process.
          if (proc != mesh->rank)
          {
            exchanger_proc_map_add_index(receive_map, proc, index);
            proc_point_map_add(point_map, proc, &x);
          }
        }
      }
    }
  }
  int_unordered_set_free(contributed_to_self);
  polymec_free(owners);

  // Sort the receive edges.
  sort_indices(point_map, receive_map);
  proc_point_map_free(point_map);

  // Now construct the exchanger.
  mesh->xy_edge_ex = exchanger_new(mesh->comm);
  exchanger_set_sends(mesh->xy_edge_ex, send_map);
  exchanger_set_receives(mesh->xy_edge_ex, receive_map);

  // By default, this exchanger uses the "min rank" reducer.
  exchanger_set_reducer(mesh->xy_edge_ex, EXCHANGER_MIN_RANK);
  STOP_FUNCTION_TIMER();
}

static void create_z_edge_ex(colmesh_t* mesh)
{
  START_FUNCTION_TIMER();

  // Figure out which processes own what chunks.
#if POLYMEC_HAVE_MPI
  int64_t* owners = source_vector(mesh);
#else
  int num_all_chunks = (int)(mesh->num_xy_chunks * mesh->num_z_chunks);
  int64_t* owners = polymec_calloc(num_all_chunks, sizeof(int64_t));
#endif

  // Determine offsets for chunks in the list of chunks for this mesh.
  int chunk_offsets[mesh->chunks->size+1];
  find_chunk_offsets(mesh, COLMESH_ZEDGE, chunk_offsets);

  // Set up send edges.
  exchanger_proc_map_t* send_map = exchanger_proc_map_new();
  proc_point_map_t* point_map = proc_point_map_new();
  int_unordered_set_t* contributed_to_self = int_unordered_set_new();
  int_unordered_set_t* processed_edge = int_unordered_set_new();
  for (size_t i = 0; i < mesh->chunks->size; ++i)
  {
    int xy = mesh->chunk_indices[2*i];
    int z = mesh->chunk_indices[2*i+1];
    int ch_index = chunk_index(mesh, xy, z);
    colmesh_chunk_t* chunk = *chunk_map_get(mesh->chunks, ch_index);
    int chunk_offset = chunk_offsets[i];
    real_t dz = (chunk->z2 - chunk->z1) / chunk->num_z_cells;

    chunk_xy_data_t* xy_data = mesh->chunk_xy_data->data[xy];
    exchanger_proc_map_t* xy_sends = xy_data->send_map;

    // Figure out the send edges.
    int pos = 0, neighbor_xy_index;
    int_array_t* indices;
    int_unordered_set_clear(contributed_to_self);
    while (exchanger_proc_map_next(xy_sends, &pos, &neighbor_xy_index, &indices))
    {
      int_unordered_set_clear(processed_edge);

      // Loop over the xy indices in the map.
      size_t num_indices = indices->size/2;
      for (size_t j = 0; j < num_indices; ++j)
      {
        int proc = (int)(owners[neighbor_xy_index]);
        int xy_edge = indices->data[2*j+1];

        // Process each one of the nodes attached to this edge, if we haven't
        // already.
        int xy_node_indices[2] = {xy_data->xy_edge_nodes[2*xy_edge], 
                                  xy_data->xy_edge_nodes[2*xy_edge+1]};
        for (int n = 0; n < 2; ++n)
        {
          int z_edge = xy_node_indices[n];
          if (!int_unordered_set_contains(processed_edge, z_edge))
          {
            // Get the x and y coordinates for the send edge.
            point2_t xe = xy_data->xy_nodes[z_edge];
            point_t x = {.x = xe.x, .y = xe.y, .z = 0.0};

            // Traverse the column and add send indices/points.
            for (int zz = 0; zz < chunk->num_z_cells; ++zz)
            {
              int index = (int)(chunk_offset + chunk->num_z_cells*z_edge + zz);
              x.z = chunk->z1 + (zz+0.5) * dz;

              // Self contribution.
              if (!int_unordered_set_contains(contributed_to_self, index))
              {
                exchanger_proc_map_add_index(send_map, mesh->rank, index);
                proc_point_map_add(point_map, mesh->rank, &x);
              }

              // Contribution to neighboring process.
              if (proc != mesh->rank)
              {
                exchanger_proc_map_add_index(send_map, proc, index);
                proc_point_map_add(point_map, proc, &x);
              }
            }
            int_unordered_set_insert(processed_edge, z_edge);
          }
        }
      }
    }
  }

  // Sort the send nodes.
  sort_indices(point_map, send_map);

  // Set up receive nodes.
  exchanger_proc_map_t* receive_map = exchanger_proc_map_new();
  proc_point_map_clear(point_map);
  for (size_t i = 0; i < mesh->chunks->size; ++i)
  {
    int xy = mesh->chunk_indices[2*i];
    int z = mesh->chunk_indices[2*i+1];
    int ch_index = chunk_index(mesh, xy, z);
    colmesh_chunk_t* chunk = *chunk_map_get(mesh->chunks, ch_index);
    int chunk_offset = chunk_offsets[i];
    real_t dz = (chunk->z2 - chunk->z1) / chunk->num_z_cells;

    chunk_xy_data_t* xy_data = mesh->chunk_xy_data->data[xy];
    exchanger_proc_map_t* xy_receives = xy_data->receive_map;

    // Figure out the xy receive edges.
    int pos = 0, neighbor_xy_index;
    int_array_t* indices;
    int_unordered_set_clear(contributed_to_self);
    while (exchanger_proc_map_next(xy_receives, &pos, &neighbor_xy_index, &indices))
    {
      int_unordered_set_clear(processed_edge);

      size_t num_indices = indices->size/2;
      for (size_t j = 0; j < num_indices; ++j)
      {
        int proc = (int)(owners[neighbor_xy_index]);
        int xy_edge = indices->data[2*j+1];

        // Process each one of the nodes attached to this edge, if we haven't
        // already.
        int xy_node_indices[2] = {xy_data->xy_edge_nodes[2*xy_edge], 
                                  xy_data->xy_edge_nodes[2*xy_edge+1]};
        for (int n = 0; n < 2; ++n)
        {
          int z_edge = xy_node_indices[n];
          if (!int_unordered_set_contains(processed_edge, z_edge))
          {
            // Get the x and y coordinates for the send edge.
            point2_t xe = xy_data->xy_nodes[z_edge];
            point_t x = {.x = xe.x, .y = xe.y, .z = 0.0};

            // Traverse the column and add receive indices/points.
            for (int zz = 0; zz < chunk->num_z_cells; ++zz)
            {
              int index = (int)(chunk_offset + chunk->num_z_cells*z_edge + zz);
              x.z = chunk->z1 + (zz+0.5) * dz;

              // Self contribution.
              if (!int_unordered_set_contains(contributed_to_self, index))
              {
                exchanger_proc_map_add_index(receive_map, mesh->rank, index);
                proc_point_map_add(point_map, mesh->rank, &x);
              }

              // Contribution to neighboring process.
              if (proc != mesh->rank)
              {
                exchanger_proc_map_add_index(receive_map, proc, index);
                proc_point_map_add(point_map, proc, &x);
              }
            }
            int_unordered_set_insert(processed_edge, z_edge);
          }
        }
      }
    }
  }
  int_unordered_set_free(processed_edge);
  int_unordered_set_free(contributed_to_self);
  polymec_free(owners);

  // Sort the receive nodes.
  sort_indices(point_map, receive_map);
  proc_point_map_free(point_map);

  // Now construct the exchanger.
  mesh->z_edge_ex = exchanger_new(mesh->comm);
  exchanger_set_sends(mesh->z_edge_ex, send_map);
  exchanger_set_receives(mesh->z_edge_ex, receive_map);

  // By default, this exchanger uses the "min rank" reducer.
  exchanger_set_reducer(mesh->z_edge_ex, EXCHANGER_MIN_RANK);
  STOP_FUNCTION_TIMER();
}

static void create_node_ex(colmesh_t* mesh)
{
  START_FUNCTION_TIMER();

  // Figure out which processes own what chunks.
#if POLYMEC_HAVE_MPI
  int64_t* owners = source_vector(mesh);
#else
  int num_all_chunks = (int)(mesh->num_xy_chunks * mesh->num_z_chunks);
  int64_t* owners = polymec_calloc(num_all_chunks, sizeof(int64_t));
#endif

  // Determine offsets for chunks in the list of chunks for this mesh.
  int chunk_offsets[mesh->chunks->size+1];
  find_chunk_offsets(mesh, COLMESH_NODE, chunk_offsets);

  // Assemble exchangers by traversing the locally stored chunks, assembling 
  // send and receive indices, and ordering those indices by the spatial locations
  // of their underlying elements. 

  // Set up send nodes.
  exchanger_proc_map_t* send_map = exchanger_proc_map_new();
  proc_point_map_t* point_map = proc_point_map_new();
  int_unordered_set_t* contributed_to_self = int_unordered_set_new();
  int_unordered_set_t* processed_node = int_unordered_set_new();
  for (size_t i = 0; i < mesh->chunks->size; ++i)
  {
    int xy = mesh->chunk_indices[2*i];
    int z = mesh->chunk_indices[2*i+1];
    int ch_index = chunk_index(mesh, xy, z);
    colmesh_chunk_t* chunk = *chunk_map_get(mesh->chunks, ch_index);
    int chunk_offset = chunk_offsets[i];
    real_t dz = (chunk->z2 - chunk->z1) / chunk->num_z_cells;

    chunk_xy_data_t* xy_data = mesh->chunk_xy_data->data[xy];
    exchanger_proc_map_t* xy_sends = xy_data->send_map;

    // Figure out the send edges.
    int pos = 0, neighbor_xy_index;
    int_array_t* indices;
    int_unordered_set_clear(contributed_to_self);
    while (exchanger_proc_map_next(xy_sends, &pos, &neighbor_xy_index, &indices))
    {
      int_unordered_set_clear(processed_node);

      // Loop over the xy indices in the map.
      size_t num_indices = indices->size/2;
      for (size_t j = 0; j < num_indices; ++j)
      {
        int proc = (int)(owners[neighbor_xy_index]);
        int edge = indices->data[2*j+1];

        // Process each one of the nodes attached to this edge, if we haven't
        // already.
        int xy_node_indices[2] = {xy_data->xy_edge_nodes[2*edge], 
                                  xy_data->xy_edge_nodes[2*edge+1]};
        for (int n = 0; n < 2; ++n)
        {
          int xy_node = xy_node_indices[n];
          if (!int_unordered_set_contains(processed_node, xy_node))
          {
            // Get the x and y coordinates for the send edge.
            point2_t xn = xy_data->xy_nodes[xy_node];
            point_t x = {.x = xn.x, .y = xn.y, .z = 0.0};

            // Traverse the column and add send indices/points.
            for (int zz = 0; zz <= chunk->num_z_cells; ++zz)
            {
              int index = (int)(chunk_offset + (chunk->num_z_cells+1)*xy_node + zz);
              x.z = chunk->z1 + zz * dz;

              // Self contribution.
              if (!int_unordered_set_contains(contributed_to_self, index))
              {
                exchanger_proc_map_add_index(send_map, mesh->rank, index);
                proc_point_map_add(point_map, mesh->rank, &x);
              }

              // Contribution to neighboring process.
              if (proc != mesh->rank)
              {
                exchanger_proc_map_add_index(send_map, proc, index);
                proc_point_map_add(point_map, proc, &x);
              }
            }
            int_unordered_set_insert(processed_node, xy_node);
          }
        }
      }
    }
  }

  // Sort the send nodes.
  sort_indices(point_map, send_map);

  // Set up receive nodes.
  exchanger_proc_map_t* receive_map = exchanger_proc_map_new();
  proc_point_map_clear(point_map);
  for (size_t i = 0; i < mesh->chunks->size; ++i)
  {
    int xy = mesh->chunk_indices[2*i];
    int z = mesh->chunk_indices[2*i+1];
    int ch_index = chunk_index(mesh, xy, z);
    colmesh_chunk_t* chunk = *chunk_map_get(mesh->chunks, ch_index);
    int chunk_offset = chunk_offsets[i];
    real_t dz = (chunk->z2 - chunk->z1) / chunk->num_z_cells;

    chunk_xy_data_t* xy_data = mesh->chunk_xy_data->data[xy];
    exchanger_proc_map_t* xy_receives = xy_data->receive_map;

    // Figure out the xy receive edges.
    int pos = 0, neighbor_xy_index;
    int_array_t* indices;
    int_unordered_set_clear(contributed_to_self);
    while (exchanger_proc_map_next(xy_receives, &pos, &neighbor_xy_index, &indices))
    {
      int_unordered_set_clear(processed_node);

      size_t num_indices = indices->size/2;
      for (size_t j = 0; j < num_indices; ++j)
      {
        int proc = (int)(owners[neighbor_xy_index]);
        int edge = indices->data[2*j+1];

        // Process each one of the nodes attached to this edge, if we haven't
        // already.
        int xy_node_indices[2] = {xy_data->xy_edge_nodes[2*edge], 
                                  xy_data->xy_edge_nodes[2*edge+1]};
        for (int n = 0; n < 2; ++n)
        {
          int xy_node = xy_node_indices[n];
          if (!int_unordered_set_contains(processed_node, xy_node))
          {
            // Get the x and y coordinates for the send edge.
            point2_t xn = xy_data->xy_nodes[xy_node];
            point_t x = {.x = xn.x, .y = xn.y, .z = 0.0};

            // Traverse the column and add receive indices/points.
            for (int zz = 0; zz <= chunk->num_z_cells; ++zz)
            {
              int index = (int)(chunk_offset + (chunk->num_z_cells+1)*xy_node + zz);
              x.z = chunk->z1 + zz * dz;

              // Self contribution.
              if (!int_unordered_set_contains(contributed_to_self, index))
              {
                exchanger_proc_map_add_index(receive_map, mesh->rank, index);
                proc_point_map_add(point_map, mesh->rank, &x);
              }

              // Contribution to neighboring process.
              if (proc != mesh->rank)
              {
                exchanger_proc_map_add_index(receive_map, proc, index);
                proc_point_map_add(point_map, proc, &x);
              }
            }
            int_unordered_set_insert(processed_node, xy_node);
          }
        }
      }
    }
  }
  int_unordered_set_free(processed_node);
  int_unordered_set_free(contributed_to_self);
  polymec_free(owners);

  // Sort the receive nodes.
  sort_indices(point_map, receive_map);
  proc_point_map_free(point_map);

  // Now construct the exchanger.
  mesh->node_ex = exchanger_new(mesh->comm);
  exchanger_set_sends(mesh->node_ex, send_map);
  exchanger_set_receives(mesh->node_ex, receive_map);

  // By default, this exchanger uses the "min rank" reducer.
  exchanger_set_reducer(mesh->node_ex, EXCHANGER_MIN_RANK);
  STOP_FUNCTION_TIMER();
}

void colmesh_finalize(colmesh_t* mesh)
{
  ASSERT(!mesh->finalized);

  // Create a sorted list of chunk indices and compute chunk index offsets.
  mesh->chunk_indices = polymec_malloc(sizeof(int) * 2 * mesh->chunks->size);
  {
    int k = 0;
    for (int i = 0; i < mesh->num_xy_chunks; ++i)
    {
      int num_z_chunks = 0;
      for (int j = 0; j < mesh->num_z_chunks; ++j)
      {
        int index = chunk_index(mesh, (int)i, (int)j);
        colmesh_chunk_t** chunk_p = chunk_map_get(mesh->chunks, index);
        if (chunk_p != NULL)
        {
          mesh->chunk_indices[2*k]   = (int)i;
          mesh->chunk_indices[2*k+1] = (int)j;
          ++k;
          ++num_z_chunks;
        }
      }
    }
    ASSERT(k == mesh->chunks->size);
  }

  // Prune unused xy data.
  for (int i = 0; i < mesh->num_xy_chunks; ++i)
  {
    int num_z_chunks = 0;
    for (int j = 0; j < mesh->num_z_chunks; ++j)
    {
      int index = chunk_index(mesh, (int)i, (int)j);
      if (chunk_map_contains(mesh->chunks, index))
        ++num_z_chunks;
    }

    // Prune xy data for this xy_index if we don't have any chunks here.
    if ((num_z_chunks == 0) && (mesh->chunk_xy_data->data[i] != NULL))
    {
      // Punch out this data on this process, and clear the corresponding 
      // destructor so it doesn't get double-freed.
      if (i < mesh->chunk_xy_data->size)
      {
        chunk_xy_data_t* xy_data = mesh->chunk_xy_data->data[i];
        if (xy_data != NULL)
        {
          chunk_xy_data_free(xy_data);
          mesh->chunk_xy_data->data[i] = NULL;
          mesh->chunk_xy_data->dtors[i] = NULL;
        }
      }
    }
  }

  // We're finished here.
  log_debug("colmesh_finalize: finalized with %d (%d x %d) local chunks.", mesh->chunks->size,
            mesh->num_xy_chunks, mesh->num_z_chunks);
  mesh->finalized = true;
}

colmesh_t* colmesh_new(MPI_Comm comm,
                       planar_polymesh_t* columns,
                       real_t z1, real_t z2,
                       int nz, bool periodic_in_z)
{
  int num_xy_chunks = 1;
  int num_z_chunks = 1; 
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
    num_z_chunks = MAX(1, (int)(pow(nproc*Nz*Nz/Nxy, 1.0/3.0)));
    
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
  int nz_per_chunk = nz / num_z_chunks;
  colmesh_t* mesh = create_empty_colmesh(comm, columns, z1, z2, 
                                         num_xy_chunks, num_z_chunks, 
                                         nz_per_chunk, periodic_in_z);
  int tot_num_chunks = num_xy_chunks * num_z_chunks;
  int chunks_per_proc = tot_num_chunks / nproc;
  for (int i = 0; i < tot_num_chunks; ++i)
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
    release_ref(mesh->cell_ex);
  if (mesh->xy_face_ex != NULL)
    release_ref(mesh->xy_face_ex);
  if (mesh->z_face_ex != NULL)
    release_ref(mesh->z_face_ex);
  if (mesh->xy_edge_ex != NULL)
    release_ref(mesh->xy_edge_ex);
  if (mesh->z_edge_ex != NULL)
    release_ref(mesh->z_edge_ex);
  if (mesh->node_ex != NULL)
    release_ref(mesh->node_ex);
  polymec_free(mesh);
}

void colmesh_get_chunk_info(colmesh_t* mesh, 
                            int* num_xy_chunks,
                            int* num_z_chunks,
                            int* nz_per_chunk)
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
  for (int c = 0; c < chunk->num_columns; ++c)
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
  for (int c = 0; c < chunk->num_columns; ++c)
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
  for (int f = 0; f < chunk->num_xy_faces; ++f)
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

int colmesh_num_chunks(colmesh_t* mesh)
{
  return (int)(mesh->chunks->size);
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
    int i = *pos;
    *xy_index = mesh->chunk_indices[2*i];
    *z_index = mesh->chunk_indices[2*i+1];
    int index = chunk_index(mesh, *xy_index, *z_index);
    if (chunk != NULL)
    {
      *chunk = *chunk_map_get(mesh->chunks, index);
      ++(*pos);
      return true;
    }
    else
      return false;
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
  size_t size = 5 * sizeof(int) + 
                (xy_data->num_columns+1) * sizeof(int) + 
                xy_data->column_xy_face_offsets[xy_data->num_columns] * sizeof(int) + 
                2 * xy_data->num_xy_faces * sizeof(int) +
                2 * xy_data->num_xy_edges * sizeof(int) + 
                xy_data->num_xy_nodes * sizeof(point2_t);
  serializer_t* ser = exchanger_proc_map_serializer();
  size += serializer_size(ser, xy_data->send_map);
  size += serializer_size(ser, xy_data->receive_map);
  return size;
}

static void* xy_data_byte_read(byte_array_t* bytes, size_t* offset)
{
  chunk_xy_data_t* xy_data = polymec_malloc(sizeof(chunk_xy_data_t));
  byte_array_read_ints(bytes, 1, &xy_data->num_columns, offset);
  byte_array_read_ints(bytes, 1, &xy_data->num_ghost_columns, offset);
  xy_data->column_xy_face_offsets = polymec_malloc((xy_data->num_columns + 1) * sizeof(int));
  byte_array_read_ints(bytes, xy_data->num_columns+1, xy_data->column_xy_face_offsets, offset);
  xy_data->column_xy_faces = polymec_malloc(xy_data->column_xy_face_offsets[xy_data->num_columns] * sizeof(int));
  byte_array_read_ints(bytes, xy_data->column_xy_face_offsets[xy_data->num_columns], xy_data->column_xy_faces, offset);
  byte_array_read_ints(bytes, 1, &xy_data->num_xy_faces, offset);
  xy_data->xy_face_columns = polymec_malloc(2*xy_data->num_xy_faces * sizeof(int));
  byte_array_read_ints(bytes, 2*xy_data->num_xy_faces, xy_data->xy_face_columns, offset);
  byte_array_read_ints(bytes, 1, &xy_data->num_xy_edges, offset);
  xy_data->xy_edge_nodes = polymec_malloc(2*xy_data->num_xy_edges * sizeof(int));
  byte_array_read_ints(bytes, 2*xy_data->num_xy_edges, xy_data->xy_edge_nodes, offset);
  byte_array_read_ints(bytes, 1, &xy_data->num_xy_nodes, offset);
  xy_data->xy_nodes = polymec_malloc(xy_data->num_xy_nodes * sizeof(point2_t));
  byte_array_read_point2s(bytes, xy_data->num_xy_nodes, xy_data->xy_nodes, offset);
  serializer_t* ser = exchanger_proc_map_serializer();
  xy_data->send_map = serializer_read(ser, bytes, offset);
  xy_data->receive_map = serializer_read(ser, bytes, offset);
  return xy_data;
}

static void xy_data_byte_write(void* obj, byte_array_t* bytes, size_t* offset)
{
  chunk_xy_data_t* xy_data = obj;
  byte_array_write_ints(bytes, 1, &xy_data->num_columns, offset);
  byte_array_write_ints(bytes, 1, &xy_data->num_ghost_columns, offset);
  byte_array_write_ints(bytes, xy_data->num_columns+1, xy_data->column_xy_face_offsets, offset);
  byte_array_write_ints(bytes, xy_data->column_xy_face_offsets[xy_data->num_columns], xy_data->column_xy_faces, offset);
  byte_array_write_ints(bytes, 1, &xy_data->num_xy_faces, offset);
  byte_array_write_ints(bytes, 2*xy_data->num_xy_faces, xy_data->xy_face_columns, offset);
  byte_array_write_ints(bytes, 1, &xy_data->num_xy_edges, offset);
  byte_array_write_ints(bytes, 2*xy_data->num_xy_edges, xy_data->xy_edge_nodes, offset);
  byte_array_write_ints(bytes, 1, &xy_data->num_xy_nodes, offset);
  byte_array_write_point2s(bytes, xy_data->num_xy_nodes, xy_data->xy_nodes, offset);
  serializer_t* ser = exchanger_proc_map_serializer();
  serializer_write(ser, xy_data->send_map, bytes, offset);
  serializer_write(ser, xy_data->receive_map, bytes, offset);
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
  for (int i = 0; i < old_mesh->num_xy_chunks; ++i)
  {
    if ((partition_vector[i] == old_mesh->rank) && (partition_vector[i] != source_vector[i]))
      insert_unique_sorted(source_procs, (int)(source_vector[i]));
    if ((source_vector[i] == old_mesh->rank) && (partition_vector[i] != source_vector[i]))
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
    for (int j = 0; j < old_mesh->num_xy_chunks; ++j)
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
    int proc = dest_procs->data[i];
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
  for (int i = 0; i < old_mesh->num_xy_chunks; ++i)
  {
    if (partition_vector[i] == old_mesh->rank)
    {
      chunk_xy_data_t* xy_data = NULL;
      if (partition_vector[i] == source_vector[i])
      {
        // Copy the xy data over directly.
        xy_data = chunk_xy_data_clone(old_mesh->chunk_xy_data->data[i]);
      }
      else
      {
        // Fetch the buffer for this process.
        size_t index = int_lower_bound(source_procs->data, source_procs->size, (int)(source_vector[i]));
        byte_array_t* buffer = receive_buffers[index]; 

        // Extract the next xy data thingy from the buffer.
        xy_data = serializer_read(ser, buffer, &(offsets[index]));
      }

      // Stick the xy data into our list.
      chunk_xy_data_array_resize(all_xy_data, MAX(all_xy_data->size, i+1));
      chunk_xy_data_array_assign_with_dtor(all_xy_data, i, xy_data, chunk_xy_data_free);
    }
  }

  // Clean up the rest of the mess.
  for (size_t i = 0; i < source_procs->size; ++i)
    byte_array_free(receive_buffers[i]);
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
  new_mesh->periodic_in_z = old_mesh->periodic_in_z;
  MPI_Comm_size(new_mesh->comm, &new_mesh->nproc);
  MPI_Comm_rank(new_mesh->comm, &new_mesh->rank);
  new_mesh->chunk_graph = adj_graph_clone(old_mesh->chunk_graph);
  new_mesh->finalized = false;
  new_mesh->cell_ex = NULL;
  new_mesh->xy_face_ex = NULL;
  new_mesh->z_face_ex = NULL;
  new_mesh->xy_edge_ex = NULL;
  new_mesh->z_edge_ex = NULL;
  new_mesh->node_ex = NULL;

  // Create xy data for chunks.
  new_mesh->chunk_xy_data = redistribute_chunk_xy_data(old_mesh, partition, sources);

  // Insert the new patches as prescribed by the partition vector.
  int num_chunks = new_mesh->num_xy_chunks * new_mesh->num_z_chunks;
  for (int i = 0; i < num_chunks; ++i)
  {
    if (partition[i] == new_mesh->rank)
    {
      int xy_index = (int)(i / new_mesh->num_z_chunks);
      int z_index = (int)(i % new_mesh->num_z_chunks);
      colmesh_insert_chunk(new_mesh, xy_index, z_index);
    }
  }

  // Finalize the new mesh.
  colmesh_finalize(new_mesh);

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
  int num_new_local_chunks = (int)colmesh_field_num_chunks(new_field);
  MPI_Request recv_requests[num_new_local_chunks];
  pos = 0;
  int num_recv_reqs = 0;
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
  int num_old_local_chunks = (int)colmesh_field_num_chunks(old_field);
  MPI_Request send_requests[num_old_local_chunks];
  pos = 0;
  int num_send_reqs = 0;
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
  MPI_Waitall(num_send_reqs, send_requests, MPI_STATUSES_IGNORE);
  MPI_Waitall(num_recv_reqs, recv_requests, MPI_STATUSES_IGNORE);

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

// These functions provide access to exchangers for colmesh_fields.
exchanger_t* colmesh_exchanger(colmesh_t* mesh, colmesh_centering_t centering);
exchanger_t* colmesh_exchanger(colmesh_t* mesh, colmesh_centering_t centering)
{
  exchanger_t* ex = NULL;
  switch (centering)
  {
    case COLMESH_CELL: 
      if (mesh->cell_ex == NULL)
        create_cell_ex(mesh);
      ex = mesh->cell_ex; 
      break;
    case COLMESH_XYFACE: 
      if (mesh->xy_face_ex == NULL)
        create_xy_face_ex(mesh);
      ex = mesh->xy_face_ex; 
      break;
    case COLMESH_ZFACE: 
      if (mesh->z_face_ex == NULL)
        create_z_face_ex(mesh);
      ex = mesh->z_face_ex; 
      break;
    case COLMESH_XYEDGE: 
      if (mesh->xy_edge_ex == NULL)
        create_xy_edge_ex(mesh);
      ex = mesh->xy_edge_ex; 
      break;
    case COLMESH_ZEDGE: 
      if (mesh->z_edge_ex == NULL)
        create_z_edge_ex(mesh);
      ex = mesh->z_edge_ex; 
      break;
    case COLMESH_NODE: 
      if (mesh->node_ex == NULL)
        create_node_ex(mesh);
      ex = mesh->node_ex; 
      break;
  }
  return ex;
}

void colmesh_chunk_z_face_get_nodes(colmesh_chunk_t* chunk,
                                    int z_face,
                                    int* nodes)
{
  // To gather nodes, we traverse the xy edges of this z face, which is 
  // the same as traversing the xy faces of the column for the z face.
  int nfaces = colmesh_chunk_column_num_xy_faces(chunk, z_face);
  int faces[nfaces];
  colmesh_chunk_column_get_xy_faces(chunk, z_face, faces);
  for (int n = 0; n < nfaces; ++n)
  {
    int face = faces[n];
    int nn = 0;
    if (face < 0)
    {
      face = ~face;
      nn = 1;
    }
    nodes[n] = chunk->xy_edge_nodes[2*face+nn];
  }
}

// These aren't part of the public API.
exchanger_proc_map_t* colmesh_xy_data_send_map(colmesh_t* mesh, int xy_index);
exchanger_proc_map_t* colmesh_xy_data_send_map(colmesh_t* mesh, int xy_index)
{
  chunk_xy_data_t* xy_data = mesh->chunk_xy_data->data[xy_index];
  return xy_data->send_map;
}

exchanger_proc_map_t* colmesh_xy_data_receive_map(colmesh_t* mesh, int xy_index);
exchanger_proc_map_t* colmesh_xy_data_receive_map(colmesh_t* mesh, int xy_index)
{
  chunk_xy_data_t* xy_data = mesh->chunk_xy_data->data[xy_index];
  return xy_data->receive_map;
}

