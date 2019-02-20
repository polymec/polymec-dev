// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/array.h"
#include "core/array_utils.h"
#include "core/timer.h"
#include "core/unordered_map.h"
#include "geometry/blockmesh.h"
#include "geometry/blockmesh_interblock_bc.h"
#include "geometry/blockmesh_field.h"
#include "geometry/blockmesh_pair.h"
#include "geometry/unimesh.h"
#include "geometry/unimesh_patch_bc.h"

#if POLYMEC_HAVE_MPI
#include "core/partitioning.h"
#endif

#if POLYMEC_HAVE_OPENMP
#include <omp.h>
#endif

DEFINE_ARRAY(unimesh_array, unimesh_t*)
DEFINE_ARRAY(bbox_array, bbox_t)

struct blockmesh_t
{
  // Parallel stuff.
  MPI_Comm comm;
  int nproc, rank;

  // Patch dimensions.
  int patch_nx, patch_ny, patch_nz;

  // Blocks and coordinate mappings.
  unimesh_array_t* blocks;
  bbox_array_t* bboxes;
  ptr_array_t* coords;

  // Inter-block boundary condition.
  blockmesh_interblock_bc_t* interblock_bc;

  // This flag is set by blockmesh_finalize() after a mesh has been assembled.
  bool finalized;
};

blockmesh_t* blockmesh_new(MPI_Comm comm, 
                           int patch_nx, 
                           int patch_ny, 
                           int patch_nz)
{
  ASSERT(patch_nx > 0);
  ASSERT(patch_ny > 0);
  ASSERT(patch_nz > 0);

  blockmesh_t* mesh = polymec_malloc(sizeof(blockmesh_t));
  mesh->comm = comm;
  MPI_Comm_size(comm, &mesh->nproc);
  MPI_Comm_rank(comm, &mesh->rank);
  mesh->patch_nx = patch_nx;
  mesh->patch_ny = patch_ny;
  mesh->patch_nz = patch_nz;
  mesh->blocks = unimesh_array_new();
  mesh->bboxes = bbox_array_new();
  mesh->coords = ptr_array_new();
  mesh->interblock_bc = blockmesh_interblock_bc_new(mesh);
  mesh->finalized = false;

  return mesh;
}

int blockmesh_add_block(blockmesh_t* mesh, 
                        bbox_t* block_domain,
                        coord_mapping_t* block_coords,
                        int num_x_patches, 
                        int num_y_patches, 
                        int num_z_patches)
{
  ASSERT(!mesh->finalized);
  ASSERT(block_domain != NULL);
  ASSERT(block_coords != NULL);
  ASSERT(coord_mapping_inverse(block_coords) != NULL);
  ASSERT(num_x_patches > 0);
  ASSERT(num_y_patches > 0);
  ASSERT(num_z_patches > 0);

  bbox_t bbox = {.x1 = 0.0, .x2 = 1.0, .y1 = 0.0, .y2 = 1.0, .z1 = 0.0, .z2 = 1.0};
  unimesh_t* block = create_empty_unimesh(mesh->comm, &bbox,
                                          num_x_patches, num_y_patches, num_z_patches,
                                          mesh->patch_nx, mesh->patch_ny, mesh->patch_nz,
                                          false, false, false);
  int index = (int)mesh->blocks->size;
  unimesh_array_append_with_dtor(mesh->blocks, block, unimesh_free);
  bbox_array_append(mesh->bboxes, *block_domain);
  ptr_array_append_with_dtor(mesh->coords, block_coords, release_ref);
  return index;
}

bool blockmesh_can_connect_blocks(blockmesh_t* mesh, 
                                  int block1_index, 
                                  int block1_nodes[4],
                                  int block2_index,
                                  int block2_nodes[4],
                                  char** reason)
{
  return blockmesh_pair_validate(mesh, block1_index, block1_nodes, 
                                 block2_index, block2_nodes, reason);
}

void blockmesh_connect_blocks(blockmesh_t* mesh, 
                              int block1_index, int block1_nodes[4],
                              int block2_index, int block2_nodes[4])
{
  ASSERT(blockmesh_can_connect_blocks(mesh, block1_index, block1_nodes,
                                            block2_index, block2_nodes, NULL));

  // Construct a block pair.
  blockmesh_pair_t* pair = blockmesh_pair_new(mesh, block1_index, block1_nodes,
                                                    block2_index, block2_nodes);

  // Traverse all the locally-stored patches in block1 and connect them to 
  // corresponding patches in block2. This is a bit grisly, since we have to 
  // account for rotations and different sets of boundary pairs.
  unimesh_t* block1 = blockmesh_block(mesh, block1_index);
  unimesh_t* block2 = blockmesh_block(mesh, block2_index);
  int pos = 0;
  int i1, j1, k1;
  while (unimesh_next_patch(block1, &pos, &i1, &j1, &k1, NULL))
  {
    // Figure out the coordinates of the corresponding patch in block2.
    int i2, j2, k2;
    blockmesh_pair_find_patch(pair, i1, j1, k1, &i2, &j2, &k2);
    if ((i2 != -1) && (j2 != -1) && (k2 != -1))
    {
      // Connect block1's local patch to block2's patch.
      blockmesh_interblock_bc_connect(mesh->interblock_bc, 
                                      block1, i1, j1, k1, 
                                      block2, i2, j2, k2, 
                                      blockmesh_pair_diffeomorphism(pair));

      // If block2's patch is locally stored, connect it to block1's.
      if (unimesh_has_patch(block2, i2, j2, k2))
      {
        blockmesh_interblock_bc_connect(mesh->interblock_bc, 
                                        block1, i1, j1, k1, 
                                        block2, i2, j2, k2,
                                        blockmesh_pair_diffeomorphism(pair));
      }
    }
  }
}

void blockmesh_assign_patches(blockmesh_t* mesh)
{
  // Count up all the patches in the global mesh.
  int num_patches = 0;
  for (size_t b = 0; b < mesh->blocks->size; ++b)
  {
    unimesh_t* block = mesh->blocks->data[b];
    int npx, npy, npz;
    unimesh_get_extents(block, &npx, &npy, &npz);
    num_patches += mesh->patch_nx * mesh->patch_ny * mesh->patch_nz * npx * npy * npz;
  }

  int_array_t* patch_list = int_array_new();
  int start_patch = 0, num_local_patches = num_patches;
#if POLYMEC_HAVE_MPI
  if (mesh->nproc > 1)
  {
    // Divide the total number of patches up amongs our processes.
    num_local_patches = num_patches / mesh->nproc;
    start_patch = mesh->rank * num_local_patches;
  }
#endif
  int patch = 0;
  for (size_t b = 0; b < mesh->blocks->size; ++b)
  {
    unimesh_t* block = mesh->blocks->data[b];
    int npx, npy, npz;
    unimesh_get_extents(block, &npx, &npy, &npz);
    for (int i = 0; i < npx; ++i)
    {
      for (int j = 0; j < npy; ++j)
      {
        for (int k = 0; k < npz; ++k)
        {
          if (patch >= start_patch)
          {
            int_array_append(patch_list, (int)b);
            int_array_append(patch_list, i);
            int_array_append(patch_list, j);
            int_array_append(patch_list, k);
            if ((int)(patch_list->size) == num_local_patches)
              goto done_selecting_patches;
          }
          ++patch;
        }
      }
    }
  }
  
done_selecting_patches:
  for (size_t p = 0; p < patch_list->size/4; ++p)
  {
    int block_index = patch_list->data[4*p];
    unimesh_t* block = mesh->blocks->data[block_index];
    int i = patch_list->data[4*p+1];
    int j = patch_list->data[4*p+2];
    int k = patch_list->data[4*p+3];
    unimesh_insert_patch(block, i, j, k);
  }
}

void blockmesh_finalize(blockmesh_t* mesh)
{
  START_FUNCTION_TIMER();
  ASSERT(!mesh->finalized);

  // Finalize the inter-block boundaries.
  blockmesh_interblock_bc_finalize(mesh->interblock_bc);

  // The boundary conditions for a given field at a block boundary depends on 
  // the structure of that field, since there's likely a change in coordinate
  // systems. These boundary conditions should be set up already. So I think we 
  // have nothing to do here except finalize the blocks within the mesh.
  for (size_t i = 0; i < mesh->blocks->size; ++i)
    unimesh_finalize(mesh->blocks->data[i]);

  mesh->finalized = true;
  STOP_FUNCTION_TIMER();
}

bool blockmesh_is_finalized(blockmesh_t* mesh)
{
  return mesh->finalized;
}

void blockmesh_free(blockmesh_t* mesh)
{
  blockmesh_interblock_bc_free(mesh->interblock_bc);
  ptr_array_free(mesh->coords);
  bbox_array_free(mesh->bboxes);
  unimesh_array_free(mesh->blocks);
  polymec_free(mesh);
}

MPI_Comm blockmesh_comm(blockmesh_t* mesh)
{
  return mesh->comm;
}

int blockmesh_num_blocks(blockmesh_t* mesh)
{
  return (int)(mesh->blocks->size);
}

unimesh_t* blockmesh_block(blockmesh_t* mesh, int index)
{
  ASSERT(index >= 0);
  ASSERT((size_t)index < mesh->blocks->size);
  return mesh->blocks->data[index];
}

bbox_t* blockmesh_block_domain(blockmesh_t* mesh, int index)
{
  return &(mesh->bboxes->data[index]);
}

coord_mapping_t* blockmesh_block_coords(blockmesh_t* mesh, int index)
{
  return mesh->coords->data[index];
}

extern bool blockmesh_interblock_bc_get_block_neighbors(blockmesh_interblock_bc_t* bc, 
                                                        int block_index, 
                                                        int block_neighbor_indices[6]);
bool blockmesh_block_is_connected(blockmesh_t* mesh,
                                  int index,
                                  unimesh_boundary_t boundary)
{
  int nblocks[6];
  blockmesh_interblock_bc_get_block_neighbors(mesh->interblock_bc, index, nblocks);
  int b = (int)boundary;
  return (nblocks[b] != -1);
}

bool blockmesh_next_block(blockmesh_t* mesh, 
                          int* pos, 
                          unimesh_t** block,
                          bbox_t* block_domain,
                          coord_mapping_t** block_coords)
{
  if (*pos < (int)mesh->blocks->size)
  {
    *block = mesh->blocks->data[*pos];
    if (block_domain != NULL)
      *block_domain = mesh->bboxes->data[*pos];
    if (block_coords != NULL)
      *block_coords = (coord_mapping_t*)(mesh->coords->data[*pos]);
    ++(*pos);
    return true;
  }
  else
    return false;
}

#if POLYMEC_HAVE_MPI
static adj_graph_t* graph_from_blocks(blockmesh_t* mesh)
{
  START_FUNCTION_TIMER();

  // Create a graph whose vertices are all of the patches within the mesh's 
  // blocks. NOTE that we associate this graph with the MPI_COMM_SELF 
  // communicator because it's a global graph.
  int num_blocks = (int)mesh->blocks->size, num_patches = 0;
  int patch_offsets[num_blocks+1]; // patch offsets by block
  patch_offsets[0] = 0;
  int pos = 0;
  unimesh_t* block;
  while (blockmesh_next_block(mesh, &pos, &block, NULL, NULL))
  {
    int npx, npy, npz;
    unimesh_get_extents(block, &npx, &npy, &npz);
    num_patches += npx * npy * npz;
    patch_offsets[pos] = num_patches;
  }
  adj_graph_t* g = adj_graph_new(MPI_COMM_SELF, num_patches);

  // Allocate storage for graph edges (patch boundaries) in the graph.
  pos = 0;
  while (blockmesh_next_block(mesh, &pos, &block, NULL, NULL))
  {
    int b = pos - 1;
    int npx, npy, npz;
    unimesh_get_extents(block, &npx, &npy, &npz);

    int nblocks[6];
    blockmesh_interblock_bc_get_block_neighbors(mesh->interblock_bc, b, nblocks);

    for (int i = 0; i < npx; ++i)
    {
      int num_x_edges = 0;
      if ((i > 0) || (nblocks[0] != -1)) ++num_x_edges;
      if ((i < npz-1) || (nblocks[1] != -1)) ++num_x_edges;
      for (int j = 0; j < npy; ++j)
      {
        int num_y_edges = 0;
        if ((j > 0) || (nblocks[2] != -1)) ++num_y_edges;
        if ((j < npy-1) || (nblocks[3] != -1)) ++num_y_edges;
        for (int k = 0; k < npz; ++k)
        {
          int num_z_edges = 0;
          if ((k > 0) || (nblocks[4] != -1)) ++num_z_edges;
          if ((k < npz-1) || (nblocks[5] != -1)) ++num_z_edges;
          int num_edges = num_x_edges + num_y_edges + num_z_edges;
          int index = patch_offsets[b] + npy*npz*i + npz*j + k;
          adj_graph_set_num_edges(g, index, num_edges);
        }
      }
    }
  }

  // Now fill in the edges.
  pos = 0;
  while (blockmesh_next_block(mesh, &pos, &block, NULL, NULL))
  {
    int b = pos - 1;
    int npx, npy, npz;
    unimesh_get_extents(block, &npx, &npy, &npz);

    int nblocks[6];
    blockmesh_interblock_bc_get_block_neighbors(mesh->interblock_bc, b, nblocks);
    int npxs[6], npys[6], npzs[6];
    for (int f = 0; f < 6; ++f)
      unimesh_get_extents(mesh->blocks->data[nblocks[f]], &npxs[f], &npys[f], &npzs[f]);

    for (int i = 0; i < npx; ++i)
    {
      for (int j = 0; j < npy; ++j)
      {
        for (int k = 0; k < npz; ++k)
        {
          int index = patch_offsets[b] + npy*npz*i + npz*j + k;
          int* edges = adj_graph_edges(g, index);
          int offset = 0;

          if ((i == 0) && (nblocks[0] != -1))
            edges[offset++] = patch_offsets[nblocks[0]] + npy*npz*(npxs[0]-1) + npz*j + k;
          else if (i > 0)
            edges[offset++] = patch_offsets[b] + npy*npz*(i-1) + npz*j + k;
          if ((i == npx-1) && (nblocks[1] != -1))
            edges[offset++] = patch_offsets[nblocks[1]] + npz*j + k;
          else if (i < npx-1)
            edges[offset++] = patch_offsets[b] + npy*npz*(i+1) + npz*j + k;

          if ((j == 0) && (nblocks[2] != -1))
            edges[offset++] = patch_offsets[nblocks[2]] + npy*npz*i + npz*(npys[1]-1) + k;
          else if (j > 0)
            edges[offset++] = patch_offsets[b] + npy*npz*i + npz*(j-1) + k;
          if ((j == npy-1) && (nblocks[3] != -1))
            edges[offset++] = patch_offsets[nblocks[3]] + npy*npz*i + k;
          else if (j < npy-1)
            edges[offset++] = patch_offsets[b] + npy*npz*i + npz*(j+1) + k;

          if ((k == 0) && (nblocks[4] != -1))
            edges[offset++] = patch_offsets[nblocks[4]] + npy*npz*i + npz*j + npzs[4]-1;
          else if (k > 0)
            edges[offset++] = patch_offsets[b] + npy*npz*i + npz*j + k-1;
          if ((k == npz-1) && (nblocks[5] != -1))
            edges[offset++] = patch_offsets[nblocks[5]] + npy*npz*i + npz*j;
          else if (k < npz-1)
            edges[offset++] = patch_offsets[b] + npy*npz*i + npz*j + k+1;
        }
      }
    }
  }

  STOP_FUNCTION_TIMER();
  return g;
}

static int64_t* source_vector(blockmesh_t* mesh)
{
  START_FUNCTION_TIMER();
  STOP_FUNCTION_TIMER();
  return NULL;
}

static void redistribute_blockmesh(blockmesh_t** mesh, 
                                   int64_t* partition)
{
  START_FUNCTION_TIMER();
  STOP_FUNCTION_TIMER();
}

// Redistributes the given block mesh using the given partition vector, but 
// does not finalize the mesh.
static void redistribute_blockmesh_field(blockmesh_field_t** field, 
                                         int64_t* partition,
                                         int64_t* sources,
                                         blockmesh_t* new_mesh)
{
  START_FUNCTION_TIMER();
  STOP_FUNCTION_TIMER();
}
#endif

void repartition_blockmesh(blockmesh_t** mesh, 
                           int* weights,
                           real_t imbalance_tol,
                           blockmesh_field_t** fields,
                           size_t num_fields)
{
  ASSERT((weights == NULL) || (imbalance_tol > 0.0));
  ASSERT((weights == NULL) || (imbalance_tol <= 1.0));
  ASSERT(imbalance_tol > 0.0);
  ASSERT(imbalance_tol <= 1.0);
  ASSERT((fields != NULL) || (num_fields == 0));
#if POLYMEC_HAVE_MPI
  START_FUNCTION_TIMER();

  // On a single process, repartitioning has no meaning.
  blockmesh_t* old_mesh = *mesh;
  if (old_mesh->nproc == 1) 
  {
    STOP_FUNCTION_TIMER();
    return;
  }

  // Generate a global adjacency graph for the mesh.
  adj_graph_t* graph = graph_from_blocks(old_mesh);

  // Map the graph to the different domains, producing a partition vector.
  // We need the partition vector on all processes, so we scatter it 
  // from rank 0.
  log_debug("repartition_blockmesh: Repartitioning mesh on %d subdomains.", old_mesh->nproc);
  int64_t* partition = partition_graph(graph, old_mesh->comm, weights, imbalance_tol, true);

  // Redistribute the mesh. 
  log_debug("repartition_blockmesh: Redistributing mesh.");
  redistribute_blockmesh(mesh, partition);
  blockmesh_finalize(*mesh);

  // Build a sources vector whose ith component is the rank that used to own 
  // the ith patch.
  int64_t* sources = source_vector(old_mesh);

  // Redistribute the fields.
  if (num_fields > 0)
    log_debug("repartition_blockmesh: Redistributing %d fields.", (int)num_fields);
  for (size_t f = 0; f < num_fields; ++f)
  {
    blockmesh_field_t* old_field = fields[f];
    redistribute_blockmesh_field(&(fields[f]), partition, sources, *mesh);
    blockmesh_field_free(old_field);
  }

  // Clean up.
  blockmesh_free(old_mesh);
  adj_graph_free(graph);
  polymec_free(sources);
  polymec_free(partition);

  STOP_FUNCTION_TIMER();
#endif
}
