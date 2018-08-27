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
  size_t num_z_cells;
  real_t z1, z2;
  int* column_xy_face_offsets;
  int* column_xy_faces;
  int* xy_face_columns;
  size_t num_xy_faces;
  size_t num_xy_edges;
  size_t num_xy_nodes;
  point2_t* xy_nodes;
} chunk_xy_data_t;

DEFINE_ARRAY(chunk_xy_data_array, chunk_xy_data_t*)

static chunk_xy_data_t* chunk_xy_data_new(planar_polymesh_t* mesh)
{
  chunk_xy_data_t* xy_data = polymec_malloc(sizeof(chunk_xy_data_t));
  xy_data->num_columns = (size_t)mesh->num_cells;

  // cell -> xy face connectivity.
  xy_data->column_xy_face_offsets = polymec_malloc(sizeof(int) * (xy_data->num_columns+1));
  memcpy(xy_data->column_xy_face_offsets, mesh->cell_edge_offsets, sizeof(int) * (xy_data->num_columns+1));

  xy_data->column_xy_faces = polymec_malloc(sizeof(int) * xy_data->column_xy_face_offsets[xy_data->num_columns]);
  memcpy(xy_data->column_xy_faces, mesh->cell_edges, sizeof(int) * xy_data->column_xy_face_offsets[xy_data->num_columns]);

  // xy face -> cell connectivity.
  xy_data->num_xy_faces = (size_t)mesh->num_edges;
  xy_data->xy_face_columns = polymec_malloc(sizeof(int) * 2 * xy_data->num_xy_faces);
  memcpy(xy_data->xy_face_columns, mesh->edge_cells, sizeof(int) * 2 * xy_data->num_xy_faces);

  // Node xy coordinates.
  xy_data->num_xy_edges = (size_t)mesh->num_edges;
  xy_data->num_xy_nodes = (size_t)mesh->num_nodes;
  xy_data->xy_nodes = polymec_malloc(sizeof(point2_t) * xy_data->num_xy_nodes);
  memcpy(xy_data->xy_nodes, mesh->nodes, sizeof(point2_t) * xy_data->num_xy_nodes);

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

static chunk_xy_data_array_t* create_chunk_xy_data(planar_polymesh_t* mesh,
                                                   int64_t* partition_vector)
{
  chunk_xy_data_array_t* all_xy_data = chunk_xy_data_array_new();
  // FIXME

  chunk_xy_data_t* xy_data = chunk_xy_data_new(mesh);
  chunk_xy_data_array_append_with_dtor(all_xy_data, xy_data, chunk_xy_data_free);
  return all_xy_data;
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
  polymec_free(chunk->xy_nodes);
  polymec_free(chunk->xy_face_columns);
  polymec_free(chunk->column_xy_faces);
  polymec_free(chunk->column_xy_face_offsets);
  polymec_free(chunk);
}

struct prismesh_t 
{
  MPI_Comm comm;
  int nproc, rank;

  planar_polymesh_t* columns;
  size_t num_xy_chunks, num_z_chunks, nz_per_chunk;
  real_t z1, z2;

  // Chunk data.
  chunk_xy_data_array_t* chunk_xy_data;
  int_ptr_unordered_map_t* chunks;
  int* chunk_indices;

  // This flag is set by prismesh_finalize() after a mesh has been assembled.
  bool finalized;

  // Exchangers for field data on each centering.
  exchanger_t *c_ex, *xyf_ex, *zf_ex, *xye_ex, *ze_ex, *n_ex;
};

prismesh_t* create_empty_prismesh(MPI_Comm comm, 
                                  planar_polymesh_t* columns,
                                  real_t z1, real_t z2,
                                  size_t num_xy_chunks, size_t num_z_chunks,
                                  size_t nz_per_chunk)
{
  ASSERT(columns != NULL);
  ASSERT(z1 < z2);
  ASSERT(num_xy_chunks > 0);
  ASSERT(num_z_chunks > 0);
  ASSERT(nz_per_chunk > 0);

  prismesh_t* mesh = polymec_malloc(sizeof(prismesh_t));
  mesh->comm = comm;
  mesh->chunks = int_ptr_unordered_map_new();
  mesh->chunk_indices = NULL;
  mesh->num_xy_chunks = num_xy_chunks;
  mesh->num_z_chunks = num_z_chunks;
  mesh->nz_per_chunk = nz_per_chunk;
  mesh->z1 = z1;
  mesh->z2 = z2;
  MPI_Comm_size(comm, &mesh->nproc);
  MPI_Comm_rank(comm, &mesh->rank);
  mesh->finalized = false;

  // Initialize exchangers. Unless we ask for anything else, we only 
  // set up the cell exchanger.
  mesh->c_ex = exchanger_new(comm);
  mesh->xyf_ex = NULL;
  mesh->zf_ex = NULL;
  mesh->xye_ex = NULL;
  mesh->ze_ex = NULL;
  mesh->n_ex = NULL;

  // Partition the planar polymesh.
#if POLYMEC_HAVE_MPI
  adj_graph_t* graph = graph_from_planar_polymesh_cells(columns);
  int64_t* P = partition_graph_n_ways(graph, num_xy_chunks, NULL, 0.05);
  adj_graph_free(graph);
#else
  int64_t* P = polymec_calloc(sizeof(int64_t) * columns->num_cells);
#endif

  // Create xy data for chunks.
  mesh->chunk_xy_data = create_chunk_xy_data(columns, P);

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
  // FIXME: Assert that we haven't already inserted this chunk!

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
  chunk->num_z_cells = nz;
  chunk->z1 = z1;
  chunk->z2 = z2;

  // Add this chunk to our list.
  int index = chunk_index(mesh, xy_index, z_index);
  int_ptr_unordered_map_insert_with_v_dtor(mesh->chunks, index, chunk, DTOR(free_chunk));
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
      if (int_ptr_unordered_map_contains(mesh->chunks, index))
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
                         size_t nz)
{
  size_t num_xy_chunks = 1;
  size_t num_z_chunks = 1; 
  int nproc, rank;
  MPI_Comm_size(comm, &nproc);
  MPI_Comm_rank(comm, &rank);
  if (nproc > 1)
  {
    // Figure out how many chunks we want in the xy and z "directions".
    // FIXME
  }

  // Now create an empty prismesh with the desired numbers of chunks, and 
  // insert all the chunks on each process. We do a "naive" placement of 
  // the chunks by allocating them sequentially to processes in a flattened 
  // index space I(xy_index, z_index) = num_z_chunks * xy_index + z_index.
  // This is definitely not ideal, but it's the easiest way to get a start.
  size_t nz_per_chunk = nz / num_z_chunks;
  prismesh_t* mesh = create_empty_prismesh(comm, columns, z1, z2, 
                                           num_xy_chunks, num_z_chunks, 
                                           nz_per_chunk);
  size_t tot_num_chunks = num_xy_chunks * num_z_chunks;
  size_t chunks_per_proc = tot_num_chunks / nproc;
  for (size_t i = 0; i < tot_num_chunks; ++i)
  {
    if ((i >= rank*chunks_per_proc) && (i < rank*(chunks_per_proc+1)))
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
  if (mesh->xyf_ex != NULL)
    polymec_release(mesh->xyf_ex);
  if (mesh->zf_ex != NULL)
    polymec_release(mesh->zf_ex);
  if (mesh->xye_ex != NULL)
    polymec_release(mesh->xye_ex);
  if (mesh->ze_ex != NULL)
    polymec_release(mesh->ze_ex);
  if (mesh->n_ex != NULL)
    polymec_release(mesh->n_ex);
  polymec_release(mesh->c_ex);

  polymec_free(mesh->chunk_indices);
  int_ptr_unordered_map_free(mesh->chunks);
  chunk_xy_data_array_free(mesh->chunk_xy_data);

  if (mesh->columns != NULL)
    planar_polymesh_free(mesh->columns);
  polymec_free(mesh);
}

MPI_Comm prismesh_comm(prismesh_t* mesh)
{
  return mesh->comm;
}

size_t prismesh_num_chunks(prismesh_t* mesh)
{
  return mesh->chunks->size;
}

size_t prismesh_num_xy_chunks(prismesh_t* mesh)
{
  return mesh->num_xy_chunks;
}

size_t prismesh_num_z_chunks(prismesh_t* mesh)
{
  return mesh->num_z_chunks;
}

bool prismesh_has_chunk(prismesh_t* mesh, int xy_index, int z_index)
{
  if ((xy_index < 0) || (xy_index >= (int)mesh->num_xy_chunks) ||
      (z_index < 0) || (z_index >= (int)mesh->num_z_chunks))
    return false;
  int index = chunk_index(mesh, xy_index, z_index);
  return int_ptr_unordered_map_contains(mesh->chunks, index);
}

real_t prismesh_z1(prismesh_t* mesh)
{
  return mesh->z1;
}

real_t prismesh_z2(prismesh_t* mesh)
{
  return mesh->z2;
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
    *chunk = *int_ptr_unordered_map_get(mesh->chunks, index);
    ++(*pos);
    return true;
  }
}

//------------------------------------------------------------------------
// Exchanger accessors. Not public, but externally available.
//------------------------------------------------------------------------
exchanger_t* prismesh_cell_exchanger(prismesh_t* mesh);
exchanger_t* prismesh_cell_exchanger(prismesh_t* mesh)
{
  return mesh->c_ex;
}

exchanger_t* prismesh_xyface_exchanger(prismesh_t* mesh);
exchanger_t* prismesh_xyface_exchanger(prismesh_t* mesh)
{
  if (mesh->xyf_ex == NULL)
  {
    // FIXME
  }
  return mesh->xyf_ex;
}

exchanger_t* prismesh_zface_exchanger(prismesh_t* mesh);
exchanger_t* prismesh_zface_exchanger(prismesh_t* mesh)
{
  if (mesh->zf_ex == NULL)
  {
    // FIXME
  }
  return mesh->zf_ex;
}

exchanger_t* prismesh_xyedge_exchanger(prismesh_t* mesh);
exchanger_t* prismesh_xyedge_exchanger(prismesh_t* mesh)
{
  if (mesh->xye_ex == NULL)
  {
    // FIXME
  }
  return mesh->xye_ex;
}

exchanger_t* prismesh_zedge_exchanger(prismesh_t* mesh);
exchanger_t* prismesh_zedge_exchanger(prismesh_t* mesh)
{
  if (mesh->ze_ex == NULL)
  {
    // FIXME
  }
  return mesh->ze_ex;
}

exchanger_t* prismesh_node_exchanger(prismesh_t* mesh);
exchanger_t* prismesh_node_exchanger(prismesh_t* mesh)
{
  if (mesh->n_ex == NULL)
  {
    // FIXME
  }
  return mesh->n_ex;
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
  new_mesh->chunks = int_ptr_unordered_map_new();
  new_mesh->chunk_indices = NULL;
  new_mesh->num_xy_chunks = old_mesh->num_xy_chunks;
  new_mesh->num_z_chunks = old_mesh->num_z_chunks;
  new_mesh->nz_per_chunk = old_mesh->nz_per_chunk;
  new_mesh->z1 = old_mesh->z1;
  new_mesh->z2 = old_mesh->z2;
  MPI_Comm_size(new_mesh->comm, &new_mesh->nproc);
  MPI_Comm_rank(new_mesh->comm, &new_mesh->rank);
  new_mesh->finalized = false;

  // Initialize exchangers. Unless we ask for anything else, we only 
  // set up the cell exchanger.
  new_mesh->c_ex = exchanger_new(new_mesh->comm);
  new_mesh->xyf_ex = NULL;
  new_mesh->zf_ex = NULL;
  new_mesh->xye_ex = NULL;
  new_mesh->ze_ex = NULL;
  new_mesh->n_ex = NULL;

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
  // Create a graph whose vertices are the mesh's patches. NOTE
  // that we associate this graph with the MPI_COMM_SELF communicator 
  // because it's a global graph.
  size_t num_chunks = mesh->num_xy_chunks * mesh->num_z_chunks;
  adj_graph_t* g = adj_graph_new(MPI_COMM_SELF, num_chunks);

#if 0
  // Allocate space in the graph for the edges (patch boundaries).
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
#if 0
  // Catalog all the chunks on this process.
  int_array_t* my_chunks = int_array_new();
  for (int i = 0; i < mesh->npx; ++i)
  {
    for (int j = 0; j < mesh->npy; ++j)
    {
      for (int k = 0; k < mesh->npz; ++k)
      {
        if (unimesh_has_patch(mesh, i, j, k))
          int_array_append(my_patches, patch_index(mesh, i, j, k));
      }
    }
  }

  // Gather the numbers of patches owned by each process.
  int num_my_patches = (int)my_patches->size;
  int num_patches_for_proc[mesh->nproc];
  MPI_Allgather(&num_my_patches, 1, MPI_INT, 
                num_patches_for_proc, 1, MPI_INT, mesh->comm);

  // Arrange for the storage of the patch indices for the patches stored 
  // on each process.
  int proc_offsets[mesh->nproc+1];
  proc_offsets[0] = 0;
  for (int p = 0; p < mesh->nproc; ++p)
    proc_offsets[p+1] = proc_offsets[p] + num_patches_for_proc[p];

  // Gæther the indices of the patches owned by all processes into a huge list.
  int num_all_patches = mesh->npx * mesh->npy * mesh->npz;
  ASSERT(num_all_patches == proc_offsets[mesh->nproc]);
  int* all_patches = polymec_malloc(sizeof(int) * num_all_patches);
  MPI_Allgatherv(my_patches->data, num_my_patches, MPI_INT, 
                 all_patches, num_patches_for_proc, proc_offsets,
                 MPI_INT, mesh->comm);

  // Clean up a bit.
  int_array_free(my_patches);

  // Convert the huge list into a source vector.
  int64_t* sources = polymec_malloc(sizeof(int64_t) * num_all_patches);
  for (int p = 0; p < mesh->nproc; ++p)
  {
    for (int offset = proc_offsets[p]; offset < proc_offsets[p+1]; ++offset)
      sources[all_patches[offset]] = (int64_t)p;
  }

  polymec_free(all_patches);
  return sources;
#endif
  return NULL;
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
#if 0

  // Copy all local chunks from one field to the other.
  unimesh_patch_t* patch;
  int pos = 0, i, j, k;
  while (unimesh_field_next_patch(new_field, &pos, &i, &j, &k, &patch, NULL))
  {
    unimesh_patch_t* old_patch = unimesh_field_patch(old_field, i, j, k);
    if (old_patch != NULL)
      unimesh_patch_copy(old_patch, patch);
  }

  // Post receives for each patch in the new field.
  int num_new_local_patches = unimesh_field_num_patches(new_field);
  MPI_Request recv_requests[num_new_local_patches];
  pos = 0;
  int num_recv_reqs = 0;
  while (unimesh_field_next_patch(new_field, &pos, &i, &j, &k, &patch, NULL))
  {
    int p = patch_index(new_mesh, i, j, k);
    if (partition[p] == new_mesh->rank)
    {
      size_t data_size = unimesh_patch_data_size(patch->centering, 
                                                 patch->nx, patch->ny, patch->nz,
                                                 patch->nc) / sizeof(real_t);
      int err = MPI_Irecv(patch->data, (int)data_size, MPI_REAL_T, (int)sources[p],
                          0, new_mesh->comm, &(recv_requests[num_recv_reqs]));
      if (err != MPI_SUCCESS)
        polymec_error("Error receiving field data from rank %d", (int)sources[p]);
      ++num_recv_reqs;
    }
  }
  ASSERT(num_recv_reqs <= num_new_local_patches);

  // Post sends.
  int num_old_local_patches = unimesh_field_num_patches(old_field);
  MPI_Request send_requests[num_old_local_patches];
  pos = 0;
  int num_send_reqs = 0;
  while (unimesh_field_next_patch(old_field, &pos, &i, &j, &k, &patch, NULL))
  {
    int p = patch_index(new_mesh, i, j, k);
    if (sources[p] == new_mesh->rank)
    {
      size_t data_size = unimesh_patch_data_size(patch->centering, 
                                                 patch->nx, patch->ny, patch->nz,
                                                 patch->nc) / sizeof(real_t);
      int err = MPI_Isend(patch->data, (int)data_size, MPI_REAL_T, (int)partition[p],
                          0, new_mesh->comm, &(send_requests[num_send_reqs]));
      if (err != MPI_SUCCESS)
        polymec_error("Error sending field data to rank %d", (int)partition[p]);
      ++num_send_reqs;
    }
  }
  ASSERT(num_send_reqs <= num_old_local_patches);

  // Wait for everything to finish.
  MPI_Waitall(num_send_reqs, send_requests, MPI_STATUSES_IGNORE);
  MPI_Waitall(num_recv_reqs, recv_requests, MPI_STATUSES_IGNORE);
#endif

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

