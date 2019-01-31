// Copyright (c) 2012-2018, Jeffrey N. Johnson
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
#include "geometry/blockmesh_field.h"
#include "geometry/unimesh.h"
#include "geometry/unimesh_patch_bc.h"

#if POLYMEC_HAVE_MPI
#include "core/partitioning.h"
#endif

#if POLYMEC_HAVE_OPENMP
#include <omp.h>
#endif

DEFINE_ARRAY(unimesh_array, unimesh_t*)
DEFINE_ARRAY(patch_bc_array, unimesh_patch_bc_t*)

struct blockmesh_t
{
  // Parallel stuff.
  MPI_Comm comm;
  int nproc, rank;

  // Patch dimensions.
  int patch_nx, patch_ny, patch_nz;

  // Blocks.
  unimesh_array_t* blocks;

  // Inter-block boundary conditions.
  patch_bc_array_t* bcs;

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
  mesh->bcs = patch_bc_array_new();
  mesh->finalized = false;

  return mesh;
}

static void free_bc(unimesh_patch_bc_t* bc)
{
  release_ref(bc);
}

// Only certain combos of block faces are acceptible.
static int _valid_block_face_nodes[6][4] = {{0, 4, 3, 7},  // -x
                                            {1, 2, 6, 5},  // +x
                                            {0, 1, 5, 4},  // -y
                                            {2, 3, 7, 6},  // +y
                                            {0, 1, 2, 3},  // -z
                                            {4, 5, 6, 7}}; // +z

int blockmesh_block_boundary_for_nodes(blockmesh_t* mesh, int block_nodes[4])
{
  int face = -1;
  for (int f = 0; f < 6; ++f)
  {
    bool face_matches = true;
    for (int n = 0; n < 4; ++n)
    {
      bool node_matches = false;
      for (int nn = 0; nn < 4; ++nn)
      {
        if (block_nodes[n] == _valid_block_face_nodes[f][nn])
        {
          node_matches = true;
          break;
        }
      }
      if (!node_matches)
        face_matches = false;
    }
    if (face_matches)
    {
      face = f;
      break;
    }
  }
  return face;
}

extern unimesh_patch_bc_t* interblock_bc_new(unimesh_t* block);
int blockmesh_add_block(blockmesh_t* mesh, 
                        int num_x_patches, 
                        int num_y_patches, 
                        int num_z_patches)
{
  ASSERT(!mesh->finalized);
  ASSERT(num_x_patches > 0);
  ASSERT(num_x_patches > 1);
  ASSERT(num_x_patches > 2);
  bbox_t bbox = {.x1 = 0.0, .x2 = 1.0, .y1 = 0.0, .y2 = 1.0, .z1 = 0.0, .z2 = 1.0};
  unimesh_t* block = create_empty_unimesh(mesh->comm, &bbox,
                                          num_x_patches, num_y_patches, num_z_patches,
                                          mesh->patch_nx, mesh->patch_ny, mesh->patch_nz,
                                          false, false, false);
  int index = (int)mesh->blocks->size;
  unimesh_array_append_with_dtor(mesh->blocks, block, unimesh_free);
  patch_bc_array_append_with_dtor(mesh->bcs, interblock_bc_new(block), free_bc);
  return index;
}

extern bool interblock_bcs_can_connect(unimesh_patch_bc_t* bc1, int bnodes1[4],
                                       unimesh_patch_bc_t* bc2, int bnodes2[4]);
bool blockmesh_blocks_can_connect(blockmesh_t* mesh, 
                                  int block1_index, 
                                  int block1_nodes[4],
                                  int block2_index,
                                  int block2_nodes[4])
{
  // Fetch the inter-block boundary conditions for the given blocks, and ask them 
  // whether they can connect.
  unimesh_patch_bc_t* bc1 = mesh->bcs->data[block1_index];
  unimesh_patch_bc_t* bc2 = mesh->bcs->data[block2_index];
  return interblock_bcs_can_connect(bc1, block1_nodes, bc2, block2_nodes);
}

extern void interblock_bcs_connect(unimesh_patch_bc_t* bc1, int bnodes1[4],
                                   unimesh_patch_bc_t* bc2, int bnodes2[4]);
void blockmesh_connect_blocks(blockmesh_t* mesh, 
                              int block1_index, int block1_nodes[4],
                              int block2_index, int block2_nodes[4])
{
  ASSERT(blockmesh_blocks_can_connect(mesh, block1_index, block1_nodes,
                                            block2_index, block2_nodes));
  unimesh_patch_bc_t* bc1 = mesh->bcs->data[block1_index];
  unimesh_patch_bc_t* bc2 = mesh->bcs->data[block2_index];
  interblock_bcs_connect(bc1, block1_nodes, bc2, block2_nodes);
}

extern void interblock_bcs_finalize(unimesh_patch_bc_t** bcs, size_t num_bcs);
void blockmesh_finalize(blockmesh_t* mesh)
{
  START_FUNCTION_TIMER();
  ASSERT(!mesh->finalized);

  // Finalize the inter-block bcs.
  interblock_bcs_finalize(mesh->bcs->data, mesh->bcs->size);

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
  patch_bc_array_free(mesh->bcs);
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

bool blockmesh_next_block(blockmesh_t* mesh, int* pos, unimesh_t** block)
{
  if (*pos < (int)mesh->blocks->size)
  {
    *block = mesh->blocks->data[*pos];
    ++(*pos);
    return true;
  }
  else
    return false;
}

#if POLYMEC_HAVE_MPI
extern void interblock_bc_get_connected_blocks(unimesh_patch_bc_t* bc, int block_index, int connected_block_indices[6]);
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
  while (blockmesh_next_block(mesh, &pos, &block))
  {
    int npx, npy, npz;
    unimesh_get_extents(block, &npx, &npy, &npz);
    num_patches += npx * npy * npz;
    patch_offsets[pos] = num_patches;
  }
  adj_graph_t* g = adj_graph_new(MPI_COMM_SELF, num_patches);

  // Allocate storage for graph edges (patch boundaries) in the graph.
  pos = 0;
  while (blockmesh_next_block(mesh, &pos, &block))
  {
    int b = pos - 1;
    int npx, npy, npz;
    unimesh_get_extents(block, &npx, &npy, &npz);

    int connected_blocks[6];
    interblock_bc_get_connected_blocks(mesh->bcs->data[b], b, connected_blocks);

    for (int i = 0; i < npx; ++i)
    {
      int num_x_edges = 0;
      if ((i > 0) || (connected_blocks[0] != -1)) ++num_x_edges;
      if ((i < npz-1) || (connected_blocks[1] != -1)) ++num_x_edges;
      for (int j = 0; j < npy; ++j)
      {
        int num_y_edges = 0;
        if ((j > 0) || (connected_blocks[2] != -1)) ++num_y_edges;
        if ((j < npy-1) || (connected_blocks[3] != -1)) ++num_y_edges;
        for (int k = 0; k < npz; ++k)
        {
          int num_z_edges = 0;
          if ((k > 0) || (connected_blocks[4] != -1)) ++num_z_edges;
          if ((k < npz-1) || (connected_blocks[5] != -1)) ++num_z_edges;
          int num_edges = num_x_edges + num_y_edges + num_z_edges;
          int index = patch_offsets[b] + npy*npz*i + npz*j + k;
          adj_graph_set_num_edges(g, index, num_edges);
        }
      }
    }
  }

  // Now fill in the edges.
  pos = 0;
  while (blockmesh_next_block(mesh, &pos, &block))
  {
    int b = pos - 1;
    int npx, npy, npz;
    unimesh_get_extents(block, &npx, &npy, &npz);

    int connected_blocks[6];
    interblock_bc_get_connected_blocks(mesh->bcs->data[b], b, connected_blocks);
    int npxs[6], npys[6], npzs[6];
    for (int f = 0; f < 6; ++f)
    {
      unimesh_t* nblock = mesh->blocks->data[connected_blocks[f]];
      unimesh_get_extents(nblock, &npxs[f], &npys[f], &npzs[f]);
    }

    for (int i = 0; i < npx; ++i)
    {
      for (int j = 0; j < npy; ++j)
      {
        for (int k = 0; k < npz; ++k)
        {
          int index = patch_offsets[b] + npy*npz*i + npz*j + k;
          int* edges = adj_graph_edges(g, index);
          int offset = 0;

          if ((i == 0) && (connected_blocks[0] != -1))
            edges[offset++] = patch_offsets[connected_blocks[0]] + npy*npz*(npxs[0]-1) + npz*j + k;
          else if (i > 0)
            edges[offset++] = patch_offsets[b] + npy*npz*(i-1) + npz*j + k;
          if ((i == npx-1) && (connected_blocks[1] != -1))
            edges[offset++] = patch_offsets[connected_blocks[1]] + npz*j + k;
          else if (i < npx-1)
            edges[offset++] = patch_offsets[b] + npy*npz*(i+1) + npz*j + k;

          if ((j == 0) && (connected_blocks[2] != -1))
            edges[offset++] = patch_offsets[connected_blocks[2]] + npy*npz*i + npz*(npys[1]-1) + k;
          else if (j > 0)
            edges[offset++] = patch_offsets[b] + npy*npz*i + npz*(j-1) + k;
          if ((j == npy-1) && (connected_blocks[3] != -1))
            edges[offset++] = patch_offsets[connected_blocks[3]] + npy*npz*i + k;
          else if (j < npy-1)
            edges[offset++] = patch_offsets[b] + npy*npz*i + npz*(j+1) + k;

          if ((k == 0) && (connected_blocks[4] != -1))
            edges[offset++] = patch_offsets[connected_blocks[4]] + npy*npz*i + npz*j + npzs[4]-1;
          else if (k > 0)
            edges[offset++] = patch_offsets[b] + npy*npz*i + npz*j + k-1;
          if ((k == npz-1) && (connected_blocks[5] != -1))
            edges[offset++] = patch_offsets[connected_blocks[5]] + npy*npz*i + npz*j;
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
