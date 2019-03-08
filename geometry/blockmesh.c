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
#include "geometry/unimesh.h"
#include "geometry/unimesh_field.h"
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

// Only certain combos of block faces are acceptible.
static int _valid_block_face_nodes[6][4] = {{0, 4, 7, 3},  // -x
                                            {2, 6, 5, 1},  // +x
                                            {0, 1, 5, 4},  // -y
                                            {7, 6, 2, 3},  // +y
                                            {0, 1, 2, 3},  // -z
                                            {7, 6, 5, 4}}; // +z

static int determine_rotation(int block1_boundary, int block1_nodes[4],
                              int block2_boundary, int block2_nodes[4])
{
  // We calculate "twists" for each pair of nodes, consisting of an integer number of
  // counter-clockwise turns from 0 to 3. If all the twists are the same, that's the
  // valid twist. Otherwise, the twist is invalid.
  int twists[4];
  for (int i = 0; i < 4; ++i)
  {
    int n1 = block1_nodes[i];
    int* valid_nodes1 = _valid_block_face_nodes[block1_boundary];
    int offset1 = (int)(int_lsearch(valid_nodes1, 4, n1) - valid_nodes1);
    int n2 = block2_nodes[i];
    int* valid_nodes2 = _valid_block_face_nodes[block2_boundary];
    int offset2 = (int)(int_lsearch(valid_nodes2, 4, n2) - valid_nodes2);
    int twist = (offset2 - offset1 + 4) % 4;
    twists[i] = twist;
    if ((i > 0) && (twists[i] != twists[0]))
      return -1;
  }
  return twists[0];
}

static int block_boundary_for_nodes(int block_nodes[4])
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

bool blockmesh_can_connect_blocks(blockmesh_t* mesh,
                                  int block1_index,
                                  int block1_nodes[4],
                                  int block2_index,
                                  int block2_nodes[4],
                                  char** reason)
{
  ASSERT(block1_index >= 0);
  ASSERT(block1_index < blockmesh_num_blocks(mesh));
  ASSERT(block2_index >= 0);
  ASSERT(block2_index < blockmesh_num_blocks(mesh));

  // This string holds the reason for failed block connections.
  static char _reason[1025];

  static const char* boundary_names[6] = {"-x", "+x", "-y", "+y", "-z", "+z"};

  int b1 = block_boundary_for_nodes(block1_nodes);
  if (b1 == -1)
  {
    if (reason != NULL)
    {
      snprintf(_reason, 1024, "First block's nodes don't correspond to a "
                              "block boundary.");
      *reason = _reason;
    }
    return false;
  }

  int b2 = block_boundary_for_nodes(block2_nodes);
  if (b2 == -1)
  {
    if (reason != NULL)
    {
      snprintf(_reason, 1024, "Second block's nodes don't correspond to a "
                              "block boundary.");
      *reason = _reason;
    }
    return false;
  }

  // A block can connect to itself, but only if the connection is between two
  // different block faces.
  if ((block1_index == block2_index) && (b1 == b2))
  {
    if (reason != NULL)
    {
      snprintf(_reason, 1024, "A block can't connect to itself via a single "
                              "boundary (%s).", boundary_names[b1]);
      *reason = _reason;
    }
    return false;
  }

  // Now make sure the two blocks are compatible on the shared boundary.
  int rotation = determine_rotation(b1, block1_nodes, b2, block2_nodes);
  if (rotation == -1)
  {
    // Try again with block2_nodes reversed.
    int rev_block2_nodes[4] = {block2_nodes[3], block2_nodes[2],
                               block2_nodes[1], block2_nodes[0]};
    rotation = determine_rotation(b1, block1_nodes, b2, rev_block2_nodes);
  }
  if (rotation == -1)
  {
    if (reason != NULL)
    {
      snprintf(_reason, 1024, "Block boundaries aren't connected in a valid "
                              "way.");
      *reason = _reason;
    }
    return false;
  }

  // Fetch the patch extents.
  unimesh_t* block1 = blockmesh_block(mesh, block1_index);
  int npx1, npy1, npz1;
  unimesh_get_extents(block1, &npx1, &npy1, &npz1);

  unimesh_t* block2 = blockmesh_block(mesh, block2_index);
  int npx2, npy2, npz2;
  unimesh_get_extents(block2, &npx2, &npy2, &npz2);

  // All patches within the mesh are the same size.
  int nx, ny, nz;
  unimesh_get_patch_size(block1, &nx, &ny, &nz);

  // Identify the relevant patch extents/sizes and make sure they're compatible.
  int N1, NP1_1, NP1_2;
  if ((b1 == 0) || (b1 == 1))
  {
    N1 = nx;
    NP1_1 = npy1;
    NP1_2 = npz1;
  }
  else if ((b1 == 2) || (b1 == 3))
  {
    N1 = ny;
    NP1_1 = npz1;
    NP1_2 = npx1;
  }
  else
  {
    N1 = nz;
    NP1_1 = npx1;
    NP1_2 = npy1;
  }

  int N2, NP2_1, NP2_2;
  if ((b2 == 0) || (b2 == 1))
  {
    N2 = nx;
    NP2_1 = npy2;
    NP2_2 = npz2;
  }
  else if ((b2 == 2) || (b2 == 3))
  {
    N2 = ny;
    NP2_1 = npz2;
    NP2_2 = npx2;
  }
  else
  {
    N2 = nz;
    NP2_1 = npx2;
    NP2_2 = npy2;
  }

  if (((rotation == 0) || (rotation == 2)) &&
      ((N1 != N2) || (NP1_1 != NP2_1) || (NP1_2 != NP2_2)))
  {
    if (reason != NULL)
    {
      snprintf(_reason, 1024, "Second block's patch extents for %s boundary "
                              "(%d and %d) don't match first block's extents "
                              "for %s boundary (%d and %d)",
               boundary_names[b2], NP2_1, NP2_2, boundary_names[b1], NP1_1,
               NP1_2);
      *reason = _reason;
    }
    return false;
  }
  else if (((rotation == 1) || (rotation == 3)) &&
           ((N1 != N2) || (NP1_1 != NP2_2) || (NP1_2 != NP2_1)))
  {
    if (reason != NULL)
    {
      snprintf(_reason, 1024, "Second block's patch extents for %s boundary "
                              "(%d and %d) don't match first block's extents "
                              "for %s boundary (%d and %d)",
               boundary_names[b2], NP2_2, NP2_1, boundary_names[b1], NP1_1,
               NP1_2);
      *reason = _reason;
    }
    return false;
  }

  // I guess everything's okay.
  return true;
}

static void find_connected_patch(blockmesh_t* mesh,
                                 int block1_index,
                                 unimesh_boundary_t block1_boundary,
                                 int i1, int j1, int k1,
                                 int rotation,
                                 int block2_index,
                                 unimesh_boundary_t block2_boundary,
                                 int* i2, int* j2, int* k2)
{
}

void blockmesh_connect_blocks(blockmesh_t* mesh,
                              int block1_index, int block1_nodes[4],
                              int block2_index, int block2_nodes[4])
{
  ASSERT(blockmesh_can_connect_blocks(mesh, block1_index, block1_nodes,
                                            block2_index, block2_nodes, NULL));

  // Figure out the block boundaries that connect to one another, and
  // the (topological) rotation going from the first to the second.
  int b1 = block_boundary_for_nodes(block1_nodes);
  ASSERT(b1 != -1);

  int b2 = block_boundary_for_nodes(block2_nodes);
  ASSERT(b2 != -1);

  int rotation = determine_rotation(b1, block1_nodes,
                                    b2, block2_nodes);
  if (rotation == -1)
  {
    // Try again with block2_nodes reversed.
    int rev_block2_nodes[4] = {block2_nodes[3], block2_nodes[2],
                               block2_nodes[1], block2_nodes[0]};
    rotation = determine_rotation(b1, block1_nodes,
                                  b2, rev_block2_nodes);
  }
  ASSERT(rotation != -1);

  // Traverse all the locally-stored patches in block1 and connect them to
  // corresponding patches in block2. This is a bit grisly, since we have to
  // account for rotations and different sets of boundary pairs.
  unimesh_t* block1 = blockmesh_block(mesh, block1_index);
  unimesh_boundary_t block1_boundary = (unimesh_boundary_t)b1;

  unimesh_t* block2 = blockmesh_block(mesh, block2_index);
  unimesh_boundary_t block2_boundary = (unimesh_boundary_t)b2;

  int pos = 0;
  int i1, j1, k1;
  while (unimesh_next_patch(block1, &pos, &i1, &j1, &k1, NULL))
  {
    // Figure out the coordinates of the corresponding patch in block2.
    int i2, j2, k2;
    find_connected_patch(mesh, block1_index, block1_boundary, i1, j1, k1,
                         rotation,
                         block2_index, block2_boundary, &i2, &j2, &k2);
    if ((i2 != -1) && (j2 != -1) && (k2 != -1))
    {
      // Connect block1's local patch to block2's patch.
      blockmesh_interblock_bc_connect(mesh->interblock_bc,
                                      block1_index, block1_boundary, i1, j1, k1,
                                      rotation,
                                      block2_index, block2_boundary, i2, j2, k2);

      // If block2's patch is locally stored, connect it to block1's.
      if (unimesh_has_patch(block2, i2, j2, k2))
      {
        int opp_rotation = (rotation + 2) % 4;
        blockmesh_interblock_bc_connect(mesh->interblock_bc,
                                        block2_index, block2_boundary, i2, j2, k2,
                                        opp_rotation,
                                        block1_index, block1_boundary, i1, j1, k1);
      }
    }
  }
}

static void assign_patches(blockmesh_t* mesh)
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
    // Divide the total number of patches up amongst our processes.
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
  int_array_free(patch_list);
}

void blockmesh_finalize(blockmesh_t* mesh)
{
  START_FUNCTION_TIMER();
  ASSERT(!mesh->finalized);

  // Assign patches to processes.
  assign_patches(mesh);

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

int blockmesh_block_index(blockmesh_t* mesh, unimesh_t* block)
{
  for (size_t b = 0; b < mesh->blocks->size; ++b)
  {
    if (block == mesh->blocks->data[b])
      return (int)b;
  }
  return -1;
}

bbox_t* blockmesh_block_domain(blockmesh_t* mesh, int index)
{
  return &(mesh->bboxes->data[index]);
}

coord_mapping_t* blockmesh_block_coords(blockmesh_t* mesh, int index)
{
  return mesh->coords->data[index];
}

void blockmesh_get_patch_size(blockmesh_t* mesh, int* nx, int* ny, int* nz)
{
  *nx = mesh->patch_nx;
  *ny = mesh->patch_ny;
  *nz = mesh->patch_nz;
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
  ASSERT(mesh->finalized);
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

  // Catalog all the patches on this process.
  int_array_t* my_patches = int_array_new();
  int num_blocks = (int)(mesh->blocks->size);
  int block_offsets[num_blocks+1];
  block_offsets[0] = 0;
  for (int b = 0; b < mesh->blocks->size; ++b)
  {
    unimesh_t* block = mesh->blocks->data[b];
    int npx, npy, npz;
    unimesh_get_extents(block, &npx, &npy, &npz);
    block_offsets[b+1] = block_offsets[b] + npx*npy*npz;

    int pos = 0, i, j, k;
    while (unimesh_next_patch(block, &pos, &i, &j, &k, NULL))
    {
      if (unimesh_has_patch(block, i, j, k))
      {
        int index = block_offsets[b] + npy*npz*i + npz*j + k;
        int_array_append(my_patches, index);
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

  // Gather the indices of the patches owned by all processes into a huge list.
  int num_all_patches = block_offsets[num_blocks];
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
  STOP_FUNCTION_TIMER();
  return sources;
}

static void redistribute_blockmesh(blockmesh_t** mesh,
                                   int64_t* partition)
{
  START_FUNCTION_TIMER();

  // Create a new mesh from the old one.
  blockmesh_t* old_mesh = *mesh;
  blockmesh_t* new_mesh = blockmesh_new(old_mesh->comm,
                                        old_mesh->patch_nx,
                                        old_mesh->patch_ny,
                                        old_mesh->patch_nz);

  // Add blocks.
  int num_blocks = (int)old_mesh->blocks->size;
  int block_offsets[num_blocks+1];
  block_offsets[0] = 0;
  for (int b = 0; b < num_blocks; ++b)
  {
    unimesh_t* block = blockmesh_block(old_mesh, b);
    int npx, npy, npz;
    unimesh_get_extents(block, &npx, &npy, &npz);
    block_offsets[b+1] = block_offsets[b] + npx*npy*npz;
    blockmesh_add_block(new_mesh, &old_mesh->bboxes->data[b],
                        old_mesh->coords->data[b], npx, npy, npz);
  }

  // Insert patches as prescribed by the partition vector.
  int num_patches = block_offsets[num_blocks];
  for (int p = 0; p < num_patches; ++p)
  {
    if (partition[p] == new_mesh->rank)
    {
      // Obtain the block index.
      int b = 0;
      while (block_offsets[b] < p) ++b;
      --b;

      unimesh_t* block = blockmesh_block(new_mesh, b);
      int npx, npy, npz;
      unimesh_get_extents(block, &npx, &npy, &npz);

      // Get the patch indices.
      int i = (p - block_offsets[b]) / (npy*npz);
      int j = (p - block_offsets[b] - npy*npz*i) / npz;
      int k = p - block_offsets[b] - npy*npz*i - npz*j;

      // Insert the patch into the block.
      unimesh_insert_patch(block, i, j, k);
    }
  }

  // Replace the old mesh with the new one.
  *mesh = new_mesh;
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

  // Create a new field from the old one.
  blockmesh_field_t* old_field = *field;
  blockmesh_field_t* new_field = blockmesh_field_new(new_mesh,
                                                     blockmesh_field_centering(old_field),
                                                     blockmesh_field_num_components(old_field));

  // Copy all local patches from one field to the other.
  int num_blocks = (int)new_mesh->blocks->size;
  int block_offsets[num_blocks+1];
  block_offsets[0] = 0;
  for (int b = 0; b < num_blocks; ++b)
  {
    unimesh_field_t* new_bfield = blockmesh_field_for_block(new_field, b);
    unimesh_field_t* old_bfield = blockmesh_field_for_block(old_field, b);
    unimesh_patch_t* patch;
    int pos = 0, i, j, k;
    while (unimesh_field_next_patch(new_bfield, &pos, &i, &j, &k, &patch, NULL))
    {
      unimesh_patch_t* old_patch = unimesh_field_patch(old_bfield, i, j, k);
      if (old_patch != NULL)
        unimesh_patch_copy(old_patch, patch);
    }

    // Compute patch offsets within blocks.
    unimesh_t* block = blockmesh_block(new_mesh, b);
    int npx, npy, npz;
    unimesh_get_extents(block, &npx, &npy, &npz);
    block_offsets[b+1] = block_offsets[b] + npx*npy*npz;
  }

  // Send and receive data one block at a time. This isn't optimal, but
  // neither is our message-per-patch strategy as a whole.
  for (int b = 0; b < num_blocks; ++b)
  {
    unimesh_t* block = blockmesh_block(new_mesh, b);
    int npx, npy, npz;
    unimesh_get_extents(block, &npx, &npy, &npz);

    // Post receives for each patch in the new field.
    unimesh_field_t* new_bfield = blockmesh_field_for_block(new_field, b);
    int num_new_local_patches = unimesh_field_num_patches(new_bfield);
    MPI_Request recv_requests[num_new_local_patches];
    int pos = 0, i, j, k;
    unimesh_patch_t* patch;
    int num_recv_reqs = 0;
    while (unimesh_field_next_patch(new_bfield, &pos, &i, &j, &k, &patch, NULL))
    {
      int p = block_offsets[b] + npy*npz*i + npz*j + k;
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
    unimesh_field_t* old_bfield = blockmesh_field_for_block(old_field, b);
    int num_old_local_patches = unimesh_field_num_patches(old_bfield);
    MPI_Request send_requests[num_old_local_patches];
    pos = 0;
    int num_send_reqs = 0;
    while (unimesh_field_next_patch(old_bfield, &pos, &i, &j, &k, &patch, NULL))
    {
      int p = block_offsets[b] + npy*npz*i + npz*j + k;
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
  }

  // Replace the old field with the new one.
  *field = new_field;
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
