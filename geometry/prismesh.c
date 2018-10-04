// Copyright (c) 2012-2018, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/timer.h"
#include "core/array.h"
#include "core/unordered_map.h"
#include "geometry/prismesh.h"
#include "geometry/prismesh_field.h"
#include "geometry/polymesh.h"

#if POLYMEC_HAVE_MPI
#include "core/partitioning.h"
#endif

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
  exchanger_t* cell_ex;
} chunk_xy_data_t;

DEFINE_ARRAY(chunk_xy_data_array, chunk_xy_data_t*)

// Creates shared xy chunk data for the local process, generating information
// for the locally-owned cells in the planar polymesh.
static chunk_xy_data_t* chunk_xy_data_new(MPI_Comm comm,
                                          planar_polymesh_t* mesh,
                                          int64_t* partition_vector,
                                          int xy_index)
{
  chunk_xy_data_t* xy_data = polymec_malloc(sizeof(chunk_xy_data_t));
  xy_data->num_columns = 0;
  xy_data->num_ghost_columns = 0;
  xy_data->num_xy_faces = 0;
  xy_data->num_xy_edges = 0;
  xy_data->num_xy_nodes = 0;

  int rank;
  MPI_Comm_rank(comm, &rank);

  // Make a list of polygonal cells (which become columns in our chunks) that 
  // correspond to the given xy chunk index.
  int_array_t* local_cells = int_array_new(); // maps columns to planar mesh cells
  int_int_unordered_map_t* cell_to_col_map = 
    int_int_unordered_map_new(); // maps planar mesh cells to columns in our chunk
  for (int cell = 0; cell < mesh->num_cells; ++cell)
  {
    if ((int)partition_vector[cell] == xy_index)
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

  // Now determine the connectivity in a single pass.
  int_array_t* face_cols = int_array_new();
  int_int_unordered_map_t* node_map = 
    int_int_unordered_map_new(); // maps planar mesh nodes to chunk nodes
  int_array_t* edge_nodes = int_array_new();
  int_array_t* ghost_xy_indices = int_array_new();
  for (int col = 0; col < xy_data->num_columns; ++col)
  {
    // Get the polygonal cell corresponding to this column.
    int cell = local_cells->data[col];

    // Loop over the neighboring cells.
    int epos = 0, edge, col_face = 0;
    while (planar_polymesh_cell_next_edge(mesh, cell, &epos, &edge))
    {
      // Create the xy face that connects these two columns (from the edge
      // connecting the two polygonal cells).
      if (mesh->edge_cells[2*edge] == cell) 
      {
        // This is our first encounter with the edge, so we create a new xy face.
        xy_data->column_xy_faces[xy_data->column_xy_face_offsets[col]+col_face] = 
          (int)(xy_data->num_xy_faces);
        ++xy_data->num_xy_faces;
        int_array_append(face_cols, col);

        // Hook up the neighboring column on the other side of our new face.
        int neighbor_cell = mesh->edge_cells[2*edge+1];
        int neighbor_col;
        bool neighbor_col_in_chunk = ((neighbor_cell != -1) && 
                                      (partition_vector[neighbor_cell] == xy_index));
        if (neighbor_col_in_chunk)
        {
          // The neighbor column is in this chunk. 
          neighbor_col = *int_int_unordered_map_get(cell_to_col_map, neighbor_cell);
        }
        else
        {
          if (neighbor_cell != -1)
          {
            // The neighbor is a ghost column, as far as this chunk is concerned.
            neighbor_col = (int)(xy_data->num_columns + ghost_xy_indices->size);
            int_array_append(ghost_xy_indices, neighbor_col);
          }
          else
          {
            // We're at an xy boundary.
            neighbor_col = -1;
          }
        }
        int_array_append(face_cols, neighbor_col);

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
          int_int_unordered_map_insert(node_map, n1, node2);
        }
        else
          node2 = *node2_p;

        // Now create an xy edge connecting these two nodes.
        int_array_append(edge_nodes, node1);
        int_array_append(edge_nodes, node2);
      }
      else
      {
        // We've already created this face, so verify that the cell corresponding 
        // to our column is on "the other side" of the planar mesh's edge.
        ASSERT(mesh->edge_cells[2*edge+1] == cell);

        // Fish out the xy face shared by this column and its neighbor.
        int neighbor_cell = mesh->edge_cells[2*edge];
        int neighbor_col = *int_int_unordered_map_get(cell_to_col_map, neighbor_cell);
        ASSERT(neighbor_col < col);

        // Find the face shared by col and neighbor_col, and add it to col's 
        // list of faces.
        bool found_col = false;
        int num_col_faces = planar_polymesh_cell_num_edges(mesh, neighbor_cell);
        for (int f = 0; f < num_col_faces; ++f)
        {
          int neighbor_col_face = xy_data->column_xy_faces[xy_data->column_xy_face_offsets[neighbor_col]+f];
          if ((neighbor_col == face_cols->data[2*neighbor_col_face]) &&
              (col == face_cols->data[2*neighbor_col_face+1]))
          {
            xy_data->column_xy_faces[xy_data->column_xy_face_offsets[col]+col_face] = neighbor_col_face;
            found_col = true;
            break;
          }
        }
        ASSERT(found_col);
      }
      ++col_face;
    }

    // Have we added a column face for every polygonal cell edge?
    ASSERT(col_face == planar_polymesh_cell_num_edges(mesh, cell));
  }

  // Surrender the data from the various arrays.
  ASSERT(xy_data->num_xy_faces == (face_cols->size / 2));
  xy_data->xy_face_columns = face_cols->data;
  int_array_release_data_and_free(face_cols);

  xy_data->xy_edge_nodes = edge_nodes->data;
  xy_data->num_xy_edges = edge_nodes->size;
  int_array_release_data_and_free(edge_nodes);

  // Set node positions.
  xy_data->num_xy_nodes = node_map->size;
  xy_data->xy_nodes = polymec_malloc(sizeof(point2_t) * node_map->size);
  int npos = 0, n, node;
  while (int_int_unordered_map_next(node_map, &npos, &n, &node))
    xy_data->xy_nodes[node] = mesh->nodes[n];
  int_int_unordered_map_free(node_map);

  // Set up a cell exchanger for this chunk.
  xy_data->num_ghost_columns = ghost_xy_indices->size;
  // FIXME
  int_array_free(ghost_xy_indices);

  // Clean up.
  int_int_unordered_map_free(cell_to_col_map);
  int_array_free(local_cells);

  return xy_data;
}

static void chunk_xy_data_free(chunk_xy_data_t* xy_data)
{
  polymec_free(xy_data->xy_nodes);
  polymec_free(xy_data->xy_face_columns);
  polymec_free(xy_data->column_xy_faces);
  polymec_free(xy_data->column_xy_face_offsets);
  polymec_free(xy_data);
}

static chunk_xy_data_array_t* redistribute_chunk_xy_data(prismesh_t* old_mesh,
                                                         int64_t* partition_vector,
                                                         int64_t* sources)
{
  chunk_xy_data_array_t* all_xy_data = chunk_xy_data_array_new();
  // FIXME
  return all_xy_data;
}

static void free_chunk(prismesh_chunk_t* chunk)
{
  // The chunk's xy data is managed by the prismesh.
  polymec_free(chunk);
}

DEFINE_UNORDERED_MAP(chunk_map, int, prismesh_chunk_t*, int_hash, int_equals)

struct prismesh_t 
{
  MPI_Comm comm;
  int nproc, rank;

  planar_polymesh_t* columns;
  size_t num_xy_chunks, num_z_chunks, nz_per_chunk;
  real_t z1, z2;

  // Chunk data.
  chunk_xy_data_array_t* chunk_xy_data;
  chunk_map_t* chunks;
  int* chunk_indices;

  // True if the mesh is periodic along the z axis, false if not.
  bool periodic_in_z;

  // This flag is set by prismesh_finalize() after a mesh has been assembled.
  bool finalized;

  // Exchangers for field data.
  int_ptr_unordered_map_t* cell_exchangers;
};

prismesh_t* create_empty_prismesh(MPI_Comm comm, 
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

  prismesh_t* mesh = polymec_malloc(sizeof(prismesh_t));
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
  adj_graph_t* graph = graph_from_planar_polymesh_cells(columns);
  int64_t* P = partition_graph_n_ways(graph, num_xy_chunks, NULL, 0.05);
  adj_graph_free(graph);
#else
  // Manually and naively partition the graph into xy chunks.
  // FIXME: If we end up experimenting with threads to do several chunks 
  // FIXME: per process, perhaps we should use a space-filling curve here?
  int64_t* P = polymec_malloc(sizeof(int64_t) * columns->num_cells);
  int num_cols_per_chunk = columns->num_cells / num_xy_chunks;
  for (int c = 0; c < columns->num_cells; ++c)
    P[c] = (int64_t)(c / num_cols_per_chunk);
#endif

  // Create xy data for chunks.
  mesh->chunk_xy_data = chunk_xy_data_array_new();
  for (int xy_index = 0; xy_index < (int)num_xy_chunks; ++xy_index)
  {
    chunk_xy_data_t* xy_data = chunk_xy_data_new(comm, columns, P, xy_index);
    chunk_xy_data_array_append_with_dtor(mesh->chunk_xy_data, xy_data, 
                                         chunk_xy_data_free);
  }

  // Initialize our exchanger maps.
  mesh->cell_exchangers = int_ptr_unordered_map_new();

  // Clean up.
  polymec_free(P);

  return mesh;
}

static inline int chunk_index(prismesh_t* mesh, int xy_index, int z_index)
{
  return (int)(mesh->num_z_chunks * xy_index + z_index);
}

void prismesh_insert_chunk(prismesh_t* mesh, int xy_index, int z_index)
{
  ASSERT(!mesh->finalized);
  ASSERT(xy_index >= 0);
  ASSERT(xy_index < mesh->num_xy_chunks);
  ASSERT(z_index >= 0);
  ASSERT(z_index < mesh->num_z_chunks);
  ASSERT(!prismesh_has_chunk(mesh, xy_index, z_index));

  prismesh_chunk_t* chunk = polymec_malloc(sizeof(prismesh_chunk_t));

  // Chunk xy data (points to data owned by mesh->chunk_xy_data).
  chunk_xy_data_t* xy_data = mesh->chunk_xy_data->data[xy_index];
  chunk->num_columns = xy_data->num_columns;
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

  // FIXME: Exchangers???
}

void prismesh_finalize(prismesh_t* mesh)
{
  ASSERT(!mesh->finalized);

  // Create a sorted list of chunk indices. And prune unused xy data.
  mesh->chunk_indices = polymec_malloc(sizeof(int) * 2 * mesh->chunks->size);
  size_t k = 0;
  for (size_t i = 0; i < mesh->num_xy_chunks; ++i)
  {
    size_t num_z_chunks = 0;
    for (size_t j = 0; j < mesh->num_z_chunks; ++j)
    {
      int index = chunk_index(mesh, (int)i, (int)j);
      if (chunk_map_contains(mesh->chunks, index))
      {
        mesh->chunk_indices[2*k]   = (int)i;
        mesh->chunk_indices[2*k+1] = (int)j;
        ++k;
        ++num_z_chunks;
      }
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
  ASSERT(k == mesh->chunks->size);

  // We're finished here.
  mesh->finalized = true;
}

prismesh_t* prismesh_new(MPI_Comm comm,
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
    // Figure out how many chunks we want in the xy and z index spaces.
    // Start by breaking the mesh up into the smallest feasible chunks.
    // Then we grow them to maximize the number of cells on each process.
    num_xy_chunks = MIN(columns->num_cells, nproc);
    num_z_chunks = MIN(nz, nproc);

#if 0
    // Our objective is one isotropic chunk per process. An isotropic chunk
    // has roughly the same number of cells in x, y, and z. In other words:
    // nx == ny == nz. Since x and y are complicated, we use the product 
    // nxny = nx*ny in our analysis. So "isotropic" means 
    // sqrt(nxny) == nz or nxny == nz*nz.
    size_t chunk_nxny = (size_t)(columns->num_cells / num_xy_chunks);
    size_t chunk_nz = nz / num_z_chunks;

    // Start growing the chunks.
    while (true)
    {
      size_t prev_num_chunks = num_xy_chunks * num_z_chunks;

      // If the aspect ratio of the chunks is skewed, make them more isotropic.
      if (chunk_nxny > (size_t)(chunk_nz*chunk_nz

      // If our number of chunks has reached nproc or our growth process has
      // stalled, terminate the iteration.
      if (((int)num_chunks == nproc) || (num_chunks == prev_num_chunks))
        break;
    }
#endif
  }

  // Now create an empty prismesh with the desired numbers of chunks, and 
  // insert all the chunks on each process. We do a "naive" placement of 
  // the chunks by allocating them sequentially to processes in a flattened 
  // index space I(xy_index, z_index) = num_z_chunks * xy_index + z_index.
  // This is definitely not ideal, but it's the easiest way to get a start.
  size_t nz_per_chunk = nz / num_z_chunks;
  prismesh_t* mesh = create_empty_prismesh(comm, columns, z1, z2, 
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
      prismesh_insert_chunk(mesh, xy_index, z_index);
    }
  }

  // Put the lid on it and ship it.
  prismesh_finalize(mesh);
  return mesh;
}

void prismesh_free(prismesh_t* mesh)
{
  int_ptr_unordered_map_free(mesh->cell_exchangers);
  polymec_free(mesh->chunk_indices);
  chunk_map_free(mesh->chunks);
  chunk_xy_data_array_free(mesh->chunk_xy_data);

  if (mesh->columns != NULL)
    planar_polymesh_free(mesh->columns);
  polymec_free(mesh);
}

void prismesh_get_chunk_info(prismesh_t* mesh, 
                             size_t* num_xy_chunks,
                             size_t* num_z_chunks,
                             size_t* nz_per_chunk)
{
  *num_xy_chunks = mesh->num_xy_chunks;
  *num_z_chunks = mesh->num_z_chunks;
  *nz_per_chunk = mesh->nz_per_chunk;
}

void prismesh_get_z_info(prismesh_t* mesh, 
                         real_t* z1,
                         real_t* z2,
                         bool* periodic)
{
  *z1 = mesh->z1;
  *z2 = mesh->z2;
  *periodic = mesh->periodic_in_z;
}

bool prismesh_verify_topology(prismesh_t* mesh, 
                              void (*handler)(const char* format, ...))
{
  // Verify the topology of each chunk.
  bool good = true;
  prismesh_chunk_t* chunk;
  int pos = 0, xy_index, z_index;
  while (prismesh_next_chunk(mesh, &pos, &xy_index, &z_index, &chunk))
  {
    bool result = prismesh_chunk_verify_topology(chunk, handler);
    if (!result)
    {
      good = false;
      break;
    }
  }
  return good;
}

bool prismesh_chunk_verify_topology(prismesh_chunk_t* chunk,
                                    void (*handler)(const char* format, ...))
{
  // All cells are prisms and must have at least 5 faces (3 xy faces + 2 z faces).
  for (size_t c = 0; c < chunk->num_columns; ++c)
  {
    if (prismesh_chunk_column_num_xy_faces(chunk, c) < 3)
    {
      handler("polymesh_verify_topology: column %d has only %d faces per cell.", 
              (int)c, prismesh_chunk_column_num_xy_faces(chunk, c) + 2);
      return false;
    }
  }

  // All z faces must have at least 3 nodes/edges.
  int num_z_faces = (int)chunk->num_columns;
  for (int f = 0; f < num_z_faces; ++f)
  {
    int nn = prismesh_chunk_z_face_num_nodes(chunk, f);
    if (nn == 0)
    {
      handler("prismesh_verify_topology: column %d has a polygonal z face with no edges!", f);
      return false;
    }
    if (nn < 3)
    {
      handler("prismesh_verify_topology: column %d has a polygonal face %d with only %d nodes.", f, nn);
      return false;
    }
  }

  // Check xy cell-face topology.
  for (size_t c = 0; c < chunk->num_columns; ++c)
  {
    int num_xy_faces = prismesh_chunk_column_num_xy_faces(chunk, c);
    int xy_faces[num_xy_faces];
    prismesh_chunk_column_get_xy_faces(chunk, c, xy_faces);
    for (int f = 0; f < num_xy_faces; ++f)
    {
      int face = xy_faces[f];
      if ((chunk->xy_face_columns[2*face] != c) && 
          (chunk->xy_face_columns[2*face+1] != c))
      {
        handler("prismesh_verify_topology: column %d has xy face %d in its list "
                "of faces, but that face does not have that column in its list.", c, xy_faces[f]);
        return false;
      }
    }
  }
  for (size_t f = 0; f < chunk->num_xy_faces; ++f)
  {
    int column = chunk->xy_face_columns[2*f];
    int num_xy_faces = prismesh_chunk_column_num_xy_faces(chunk, column);
    int xy_faces[num_xy_faces];
    prismesh_chunk_column_get_xy_faces(chunk, column, xy_faces);
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
      handler("prismesh_verify_topology: xy face %d has column %d in its list of columns, but "
              "that column does not have that face in its list.", f, chunk->xy_face_columns[2*f]);
      return false;
    }
    if (chunk->xy_face_columns[2*f+1] != -1)
    {
      found_face = false;
      int other_col = chunk->xy_face_columns[2*f+1];
      num_xy_faces = prismesh_chunk_column_num_xy_faces(chunk, other_col);
      prismesh_chunk_column_get_xy_faces(chunk, other_col, xy_faces);
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
        handler("prismesh_verify_topology: xy face %d has column %d in its list of columns, "
                "but that column does not have that xy face in its list of faces.", f, other_col);
        return false;
      }
    }
  }
  return true;
}

MPI_Comm prismesh_comm(prismesh_t* mesh)
{
  return mesh->comm;
}

size_t prismesh_num_chunks(prismesh_t* mesh)
{
  return mesh->chunks->size;
}

prismesh_chunk_t* prismesh_chunk(prismesh_t* mesh, int xy_index, int z_index)
{
  int index = chunk_index(mesh, xy_index, z_index);
  prismesh_chunk_t** chunk_p = chunk_map_get(mesh->chunks, index);
  return (chunk_p != NULL) ? *chunk_p : NULL;
}

bool prismesh_has_chunk(prismesh_t* mesh, int xy_index, int z_index)
{
  return (prismesh_chunk(mesh, xy_index, z_index) != NULL);
}

polygon_t* prismesh_chunk_polygon(prismesh_chunk_t* chunk, int column)
{
  ASSERT(column < (int)chunk->num_columns);
  int z_face = column;
  int num_nodes = prismesh_chunk_z_face_num_nodes(chunk, z_face);
  int nodes[num_nodes];
  prismesh_chunk_z_face_get_nodes(chunk, z_face, nodes);
  point2_t vertices[num_nodes];
  for (int n = 0; n < num_nodes; ++n)
    vertices[n] = chunk->xy_nodes[nodes[n]];
  return polygon_new(vertices, num_nodes);
}

bool prismesh_next_chunk(prismesh_t* mesh, int* pos, 
                         int* xy_index, int* z_index,
                         prismesh_chunk_t** chunk)
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

#if POLYMEC_HAVE_MPI
static void redistribute_prismesh(prismesh_t** mesh, 
                                  int64_t* partition,
                                  int64_t* sources)
{
  START_FUNCTION_TIMER();
  prismesh_t* old_mesh = *mesh;

  // Create a new empty mesh.
  prismesh_t* new_mesh = polymec_malloc(sizeof(prismesh_t));
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
      prismesh_insert_chunk(new_mesh, xy_index, z_index);
    }
  }

  // Replace the old mesh with the new one.
  *mesh = new_mesh;
  STOP_FUNCTION_TIMER();
}

static adj_graph_t* graph_from_prismesh_chunks(prismesh_t* mesh)
{
  // Create a graph whose vertices are the mesh's chunks. NOTE
  // that we associate this graph with the MPI_COMM_SELF communicator 
  // because it's a global graph.
  size_t num_chunks = mesh->num_xy_chunks * mesh->num_z_chunks;
  adj_graph_t* g = adj_graph_new(MPI_COMM_SELF, num_chunks);

#if 0
  // Allocate space in the graph for the edges (chunk boundaries).
  for (int i = 0; i < mesh->npx; ++i)
  {
    int num_x_edges = (i == 0) ? (i == mesh->npx-1) ? mesh->periodic_in_x ? 2 
                                                                          : 0
                                                    : 1
                               : (i == mesh->npx-1) ? 1
                                                    : 2;
    for (int j = 0; j < mesh->npy; ++j)
    {
      int num_y_edges = (j == 0) ? (j == mesh->npy-1) ? mesh->periodic_in_y ? 2 
                                                                            : 0
                                                      : 1
                                 : (j == mesh->npy-1) ? 1
                                                      : 2;
      for (int k = 0; k < mesh->npz; ++k)
      {
        int num_z_edges = (k == 0) ? (k == mesh->npz-1) ? mesh->periodic_in_z ? 2 
                                                                              : 0
                                                        : 1
                                   : (k == mesh->npz-1) ? 1
                                                        : 2;
        int num_edges = num_x_edges + num_y_edges + num_z_edges;
        int p_index = patch_index(mesh, i, j, k);
        adj_graph_set_num_edges(g, p_index, num_edges);
      }
    }
  }

  // Now fill in the edges.
  for (int i = 0; i < mesh->npx; ++i)
  {
    for (int j = 0; j < mesh->npy; ++j)
    {
      for (int k = 0; k < mesh->npz; ++k)
      {
        int p_index = patch_index(mesh, i, j, k);
        int* edges = adj_graph_edges(g, p_index);
        int offset = 0;

        if ((i == 0) && mesh->periodic_in_x)
          edges[offset++] = patch_index(mesh, mesh->npx-1, j, k);
        else if (i > 0)
          edges[offset++] = patch_index(mesh, i-1, j, k);
        if ((i == mesh->npx-1) && mesh->periodic_in_x)
          edges[offset++] = patch_index(mesh, 0, j, k);
        else if (i < mesh->npx-1)
          edges[offset++] = patch_index(mesh, i+1, j, k);

        if ((j == 0) && mesh->periodic_in_y)
          edges[offset++] = patch_index(mesh, i, mesh->npy-1, k);
        else if (j > 0)
          edges[offset++] = patch_index(mesh, i, j-1, k);
        if ((j == mesh->npy-1) && mesh->periodic_in_y)
          edges[offset++] = patch_index(mesh, i, 0, k);
        else if (j < mesh->npy-1)
          edges[offset++] = patch_index(mesh, i, j+1, k);

        if ((k == 0) && mesh->periodic_in_z)
          edges[offset++] = patch_index(mesh, i, j, mesh->npz-1);
        else if (k > 0)
          edges[offset++] = patch_index(mesh, i, j, k-1);
        if ((k == mesh->npz-1) && mesh->periodic_in_z)
          edges[offset++] = patch_index(mesh, i, j, 0);
        else if (k < mesh->npz-1)
          edges[offset++] = patch_index(mesh, i, j, k+1);
      }
    }
  }
#endif

  return g;
}

static int64_t* source_vector(prismesh_t* mesh)
{
  // Catalog all the chunks on this process.
  int_array_t* my_chunks = int_array_new();
  for (int xy_index = 0; xy_index < (int)mesh->num_xy_chunks; ++xy_index)
  {
    for (int z_index = 0; z_index < (int)mesh->num_z_chunks; ++z_index)
    {
      if (prismesh_has_chunk(mesh, xy_index, z_index))
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

static void redistribute_prismesh_field(prismesh_field_t** field, 
                                        int64_t* partition,
                                        int64_t* sources,
                                        prismesh_t* new_mesh)
{
  START_FUNCTION_TIMER();

  // Create a new field from the old one.
  prismesh_field_t* old_field = *field;
  prismesh_field_t* new_field = prismesh_field_new(new_mesh,
                                                   prismesh_field_centering(old_field),
                                                   prismesh_field_num_components(old_field));

  // Copy all local chunk data from one field to the other.
  prismesh_chunk_data_t* data;
  int pos = 0, xy_index, z_index;
  while (prismesh_field_next_chunk(new_field, &pos, &xy_index, &z_index, &data))
  {
    prismesh_chunk_data_t* old_data = prismesh_field_chunk_data(old_field, xy_index, z_index);
    if (old_data != NULL)
      prismesh_chunk_data_copy(old_data, data);
  }

  // Post receives for each chunk in the new field.
  size_t num_new_local_chunks = prismesh_field_num_chunks(new_field);
  MPI_Request recv_requests[num_new_local_chunks];
  pos = 0;
  size_t num_recv_reqs = 0;
  while (prismesh_field_next_chunk(new_field, &pos, &xy_index, &z_index, &data))
  {
    int ch = chunk_index(new_mesh, xy_index, z_index);
    if (partition[ch] == new_mesh->rank)
    {
      size_t data_size = prismesh_chunk_data_size(data) / sizeof(real_t);
      int err = MPI_Irecv(data->data, (int)data_size, MPI_REAL_T, (int)sources[ch],
                          0, new_mesh->comm, &(recv_requests[num_recv_reqs]));
      if (err != MPI_SUCCESS)
        polymec_error("Error receiving field data from rank %d", (int)sources[ch]);
      ++num_recv_reqs;
    }
  }
  ASSERT(num_recv_reqs <= num_new_local_chunks);

  // Post sends.
  size_t num_old_local_chunks = prismesh_field_num_chunks(old_field);
  MPI_Request send_requests[num_old_local_chunks];
  pos = 0;
  size_t num_send_reqs = 0;
  while (prismesh_field_next_chunk(old_field, &pos, &xy_index, &z_index, &data))
  {
    int ch = chunk_index(new_mesh, xy_index, z_index);
    if (sources[ch] == new_mesh->rank)
    {
      size_t data_size = prismesh_chunk_data_size(data) / sizeof(real_t);
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

void repartition_prismesh(prismesh_t** mesh, 
                          int* weights,
                          real_t imbalance_tol,
                          prismesh_field_t** fields,
                          size_t num_fields)
{
  ASSERT((weights == NULL) || (imbalance_tol > 0.0));
  ASSERT((weights == NULL) || (imbalance_tol <= 1.0));
  ASSERT((fields != NULL) || (num_fields == 0));
#if POLYMEC_HAVE_MPI
  START_FUNCTION_TIMER();

  // On a single process, repartitioning has no meaning.
  prismesh_t* old_mesh = *mesh;
  if (old_mesh->nproc == 1) 
  {
    STOP_FUNCTION_TIMER();
    return;
  }

  // Generate a distributed adjacency graph for the mesh's chunks.
  adj_graph_t* graph = graph_from_prismesh_chunks(old_mesh);

  // Map the graph to the different domains, producing a partition vector.
  // We need the partition vector on all processes in the communicator, so we 
  // scatter it from rank 0.
  log_debug("repartition_prismesh: Repartitioning mesh on %d subdomains.", old_mesh->nproc);
  int64_t* P = partition_graph(graph, old_mesh->comm, weights, imbalance_tol, true);

  // Build a sources vector whose ith component is the rank that used to own 
  // the ith patch.
  int64_t* sources = source_vector(old_mesh);

  // Redistribute the mesh. 
  log_debug("repartition_peximesh: Redistributing mesh.");
  redistribute_prismesh(mesh, P, sources);

  // Redistribute the fields.
  if (num_fields > 0)
    log_debug("repartition_prismesh: Redistributing %d fields.", (int)num_fields);
  for (size_t f = 0; f < num_fields; ++f)
  {
    prismesh_field_t* old_field = fields[f];
    redistribute_prismesh_field(&(fields[f]), P, sources, *mesh);
    prismesh_field_free(old_field);
  }

  // Clean up.
  prismesh_free(old_mesh);
  adj_graph_free(graph);
  polymec_free(sources);
  polymec_free(P);

  STOP_FUNCTION_TIMER();
#endif
}

polymesh_t* prismesh_as_polymesh(prismesh_t* mesh)
{
  return NULL; // FIXME
}

// These functions provide access to exchangers for prismesh_fields.
exchanger_t* prismesh_chunk_cell_exchanger(prismesh_t* mesh, int xy_index, int z_index);
exchanger_t* prismesh_chunk_cell_exchanger(prismesh_t* mesh, int xy_index, int z_index)
{
  int index = chunk_index(mesh, xy_index, z_index);
  return *int_ptr_unordered_map_get(mesh->cell_exchangers, index);
}
