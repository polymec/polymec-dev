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

#if POLYMEC_HAVE_MPI
#include "core/partitioning.h"
#endif

#if POLYMEC_HAVE_OPENMP
#include <omp.h>
#endif

DEFINE_ARRAY(unimesh_array, unimesh_t*)

// This type represents counterclockwise "twists" between connected blocks.
typedef enum
{
  NO_TURN = 0,
  QUARTER_TURN = 1,
  HALF_TURN = 2,
  THREE_QUARTERS_TURN = 3,
  INVALID_TWIST
} cxn_twist_t;

// This type represents a connection between two blocks.
typedef struct
{
  int blocks[6];
  int boundaries[6]; 
  cxn_twist_t twists[6];
} cxn_t;

static cxn_t* cxn_new()
{
  cxn_t* cxn = polymec_malloc(sizeof(cxn_t));
  for (int b = 0; b < 6; ++b)
  {
    cxn->blocks[b] = -1;
    cxn->boundaries[b] = -1;
    cxn->twists[b] = INVALID_TWIST;
  }
  return cxn;
}

static void cxn_free(cxn_t* cxn)
{
  polymec_free(cxn);
}

DEFINE_UNORDERED_MAP(cxn_map, int, cxn_t*, int_hash, int_equals)

typedef struct
{
  // Parallel stuff.
  MPI_Comm comm;
  int nproc, rank;

  // Mapping of block indices to connection structs.
  cxn_map_t* cxns;

  // This flag is set by blockmesh_finalize() after a mesh has been assembled.
  bool finalized;
} interblock_bc_t;

unimesh_patch_bc_t* intermesh_bc_new(unimesh_t* block)
{
  return NULL;
}

static cxn_twist_t determine_twist(int block1_face,
                                   int block1_nodes[4],
                                   int block2_face,
                                   int block2_nodes[4])
{
  // By the time this function gets called, we know that both sets of nodes 
  // correspond to valid block faces. So we need only identify which nodes 
  // are identified between the two block faces.

  // We calculate "twists" for each pair of nodes. If all the twists are the same,
  // that's the valid twist. Otherwise, the twist is invalid.
  cxn_twist_t twists[4];
  for (int i = 0; i < 4; ++i)
  {
    int n1 = block1_nodes[i];
    int* valid_nodes1 = _valid_block_face_nodes[block1_face];
    int offset1 = (int)(int_lsearch(valid_nodes1, 4, n1) - valid_nodes1);
    int n2 = block2_nodes[i];
    int* valid_nodes2 = _valid_block_face_nodes[block2_face];
    int offset2 = (int)(int_lsearch(valid_nodes2, 4, n2) - valid_nodes2);
    int diff = (offset2 - offset1 + 4) % 4;
    twists[i] = (cxn_twist_t)diff;
    if ((i > 0) && (twists[i] != twists[0]))
      return INVALID_TWIST;
  }
  return twists[0];
}

bool interblock_bcs_can_connect(unimesh_patch_bc_t* bc1,
                                int bnodes1[4],
                                unimesh_patch_bc_t* bc2,
                                int bnodes2[4])
{
  int b1 = blockmesh_block_boundary_for_nodes(mesh, block1_nodes);
  if (b1 == -1) return false;
  int b2 = blockmesh_block_boundary_for_nodes(mesh, block2_nodes);
  if (b2 == -1) return false;

  // A block can connect to itself, but only if the connection is between two 
  // different block faces.
  if ((bc1 == bc2) && (b1 == b2))
    return false;

  // Now make sure the patch dimensions between the two blocks are compatible on 
  // the shared boundary.
  cxn_twist_t twist = determine_twist(b1, block1_nodes, b2, block2_nodes);
  if (twist == INVALID_TWIST)
    return false;

  // Fetch the patch extents.
  int npx1, npy1, npz1;
  unimesh_get_extents(bc1->block, &npx1, &npy1, &npz1);
  int npx2, npy2, npz2;
  unimesh_get_extents(bc2->block, &npx2, &npy2, &npz2);

  // All patches within the mesh are the same size.
  int nx, ny, nz;
  unimesh_get_patch_size(bc1->block, &nx, &ny, &nz);

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

  if (((twist == NO_TURN) || (twist == HALF_TURN)) && 
      ((N1 != N2) || (NP1_1 != NP2_1) || (NP1_2 != NP2_2)))
    return false;
  else if ((N1 != N2) || (NP1_1 != NP2_2) || (NP1_2 != NP2_1))
    return false;

  // I guess everything's okay.
  return true;
}

void interblock_bcs_connect(unimesh_patch_bc_t* bc1, int bnodes1[4],
                            unimesh_patch_bc_t* bc2, int bnodes2[4])
{
  ASSERT(interblock_bcs_can_connect(bc1, bnodes1, bc2, bnodes2));
  int b1 = blockmesh_block_boundary_for_nodes(mesh, bnodes1);
  int b2 = blockmesh_block_boundary_for_nodes(mesh, bnodes2);
  cxn_twist_t twist = determine_twist(b1, block1_nodes, b2, block2_nodes);

  // Find or add a connection for each block.
  cxn_t** cxn1_p = cxn_map_get(mesh->cxns, block1_index);
  cxn_t* cxn1 = NULL;
  if (cxn1_p == NULL)
  {
    cxn1 = cxn_new();
    cxn_map_insert_with_v_dtor(mesh->cxns, block1_index, cxn1, cxn_free);
  }
  else
    cxn1 = *cxn1_p;
  cxn1->blocks[b1] = block2_index;
  cxn1->boundaries[b1] = b2;
  cxn1->twists[b1] = (int)twist;

  cxn_t** cxn2_p = cxn_map_get(mesh->cxns, block2_index);
  cxn_t* cxn2 = NULL;
  if (cxn2_p == NULL)
  {
    cxn2 = cxn_new();
    cxn_map_insert_with_v_dtor(mesh->cxns, block2_index, cxn2, cxn_free);
  }
  else
    cxn2 = *cxn2_p;
  cxn2->blocks[b2] = block1_index;
  cxn2->boundaries[b2] = b1;

  // The twist for block1 w.r.t. block2 is reversed.
  if ((twist == NO_TURN) || (twist == HALF_TURN))
    cxn2->twists[b2] = twist;
  else 
    cxn2->twists[b2] = ((((int)twist) + 2) % 4);
}

void interblock_bcs_finalize(unimesh_patch_bc_t** bcs, size_t num_bcs)
{
  START_FUNCTION_TIMER();
  STOP_FUNCTION_TIMER();
}

