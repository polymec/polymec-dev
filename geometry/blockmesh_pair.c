// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/array_utils.h"
#include "geometry/blockmesh.h"
#include "geometry/blockmesh_pair.h"
#include "geometry/unimesh_patch.h"

// This represents a rotation that must be performed for a quantity to pass 
// through the face connecting the two blocks. Expressed in units of 
// counterclockwise revolutions (turns).
typedef enum 
{
  NO_ROTATION,
  QUARTER_TURN,
  HALF_TURN,
  THREE_QUARTERS_TURN,
  INVALID_ROTATION
} block_rotation_t;

// Only certain combos of block faces are acceptible.
static int _valid_block_face_nodes[6][4] = {{0, 4, 7, 3},  // -x
                                            {2, 6, 5, 1},  // +x
                                            {0, 1, 5, 4},  // -y
                                            {7, 6, 2, 3},  // +y
                                            {0, 1, 2, 3},  // -z
                                            {7, 6, 5, 4}}; // +z

static block_rotation_t determine_rotation(bbox_t* block1_domain,
                                           coord_mapping_t* block1_coords,
                                           int block1_boundary,
                                           int block1_nodes[4],
                                           bbox_t* block2_domain,
                                           coord_mapping_t* block2_coords,
                                           int block2_boundary,
                                           int block2_nodes[4])
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
      return INVALID_ROTATION; 
  }
  block_rotation_t rotation = INVALID_ROTATION;
  switch (twists[0])
  {
    case 0: rotation = NO_ROTATION; break;
    case 1: rotation = QUARTER_TURN; break;
    case 2: rotation = HALF_TURN; break;
    case 3: rotation = THREE_QUARTERS_TURN;
  }
  return rotation;
}

struct blockmesh_pair_t 
{
  unimesh_t* block1;
  unimesh_t* block2;
  bbox_t block1_domain, block2_domain;
  coord_mapping_t* block1_coords;
  coord_mapping_t* block2_coords;
  unimesh_boundary_t block1_boundary, block2_boundary;
  int npx1, npy1, npz1;
  int npx2, npy2, npz2;
  block_rotation_t rotation;
};

static void blockmesh_pair_free(void* context)
{
  blockmesh_pair_t* pair = context;
  release_ref(pair->block1_coords);
  release_ref(pair->block2_coords);
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

bool blockmesh_pair_validate(blockmesh_t* mesh, 
                             int block1_index, int block1_nodes[4],
                             int block2_index, int block2_nodes[4],
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
      snprintf(_reason, 1024, "First block's nodes don't correspond to a block boundary.");
      *reason = _reason;
    }
    return false;
  }

  int b2 = block_boundary_for_nodes(block2_nodes);
  if (b2 == -1)
  {
    if (reason != NULL)
    {
      snprintf(_reason, 1024, "Second block's nodes don't correspond to a block boundary.");
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
      snprintf(_reason, 1024, "A block can't connect to itself via a single boundary (%s).", boundary_names[b1]);
      *reason = _reason;
    }
    return false;
  }

  // Now make sure the two blocks are compatible on the shared boundary.
  bbox_t* block1_domain = blockmesh_block_domain(mesh, block1_index);
  coord_mapping_t* block1_coords = blockmesh_block_coords(mesh, block1_index);
  bbox_t* block2_domain = blockmesh_block_domain(mesh, block2_index);
  coord_mapping_t* block2_coords = blockmesh_block_coords(mesh, block2_index);

  block_rotation_t rotation = determine_rotation(block1_domain, block1_coords, b1, block1_nodes, 
                                                 block1_domain, block2_coords, b2, block2_nodes);
  if (rotation == INVALID_ROTATION)
  {
    // Try again with block2_nodes reversed.
    int rev_block2_nodes[4] = {block2_nodes[3], block2_nodes[2], 
                               block2_nodes[1], block2_nodes[0]};
    rotation = determine_rotation(block1_domain, block1_coords, b1, block1_nodes, 
                                  block2_domain, block2_coords, b2, rev_block2_nodes);
  }
  if (rotation == INVALID_ROTATION)
  {
    if (reason != NULL)
    {
      snprintf(_reason, 1024, "Block boundaries aren't connected in a valid way.");
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

  if (((rotation == NO_ROTATION) || (rotation == HALF_TURN)) && 
      ((N1 != N2) || (NP1_1 != NP2_1) || (NP1_2 != NP2_2)))
  {
    if (reason != NULL)
    {
      snprintf(_reason, 1024, "Second block's patch extents for %s boundary (%d and %d) "
                              "don't match first block's extents for %s boundary (%d and %d)", 
               boundary_names[b2], NP2_1, NP2_2, boundary_names[b1], NP1_1, NP1_2);
      *reason = _reason;
    }
    return false;
  }
  else if (((rotation == QUARTER_TURN) || (rotation == THREE_QUARTERS_TURN)) && 
           ((N1 != N2) || (NP1_1 != NP2_2) || (NP1_2 != NP2_1)))
  {
    if (reason != NULL)
    {
      snprintf(_reason, 1024, "Second block's patch extents for %s boundary (%d and %d) "
                              "don't match first block's extents for %s boundary (%d and %d)", 
               boundary_names[b2], NP2_2, NP2_1, boundary_names[b1], NP1_1, NP1_2);
      *reason = _reason;
    }
    return false;
  }

  // I guess everything's okay.
  return true;
}

blockmesh_pair_t* blockmesh_pair_new(blockmesh_t* mesh, 
                                     int block1_index, int block1_nodes[4],
                                     int block2_index, int block2_nodes[4])
{
  bool valid = blockmesh_pair_validate(mesh, block1_index, block1_nodes, 
                                       block2_index, block2_nodes, NULL);
  if (!valid)
    return NULL;

  blockmesh_pair_t* pair = polymec_refcounted_malloc(sizeof(blockmesh_pair_t), 
                                                     blockmesh_pair_free);
  unimesh_t* block1 = blockmesh_block(mesh, block1_index);
  pair->block1 = block1;
  unimesh_get_extents(block1, &pair->npx1, &pair->npy1, &pair->npz1);
  pair->block1_domain = *blockmesh_block_domain(mesh, block1_index);
  pair->block1_coords = blockmesh_block_coords(mesh, block1_index);
  retain_ref(pair->block1_coords);

  unimesh_t* block2 = blockmesh_block(mesh, block2_index);
  pair->block2 = block2;
  pair->block2_domain = *blockmesh_block_domain(mesh, block2_index);
  pair->block2_coords = blockmesh_block_coords(mesh, block2_index);
  retain_ref(pair->block2_coords);

  int b1 = block_boundary_for_nodes(block1_nodes);
  ASSERT(b1 != -1);
  pair->block1_boundary = (unimesh_boundary_t)b1;
  int b2 = block_boundary_for_nodes(block2_nodes);
  ASSERT(b2 != -1);
  pair->block2_boundary = (unimesh_boundary_t)b2;
  pair->rotation = determine_rotation(&pair->block1_domain, pair->block1_coords, pair->block1_boundary, block1_nodes,
                                      &pair->block2_domain, pair->block2_coords, pair->block2_boundary, block2_nodes);
  if (pair->rotation == INVALID_ROTATION)
  {
    // Try again with block2_nodes reversed.
    int rev_block2_nodes[4] = {block2_nodes[3], block2_nodes[2], 
                               block2_nodes[1], block2_nodes[0]};
    pair->rotation = determine_rotation(&pair->block1_domain, pair->block1_coords, (unimesh_boundary_t)b1, block1_nodes,
                                        &pair->block2_domain, pair->block2_coords, (unimesh_boundary_t)b2, rev_block2_nodes);
  }
  ASSERT(pair->rotation != INVALID_ROTATION);
  return pair;
}

unimesh_t* blockmesh_pair_block1(blockmesh_pair_t* pair)
{
  return pair->block1;
}

unimesh_boundary_t blockmesh_pair_block1_boundary(blockmesh_pair_t* pair)
{
  return pair->block1_boundary;
}

unimesh_t* blockmesh_pair_block2(blockmesh_pair_t* pair)
{
  return pair->block2;
}

unimesh_boundary_t blockmesh_pair_block2_boundary(blockmesh_pair_t* pair)
{
  return pair->block2_boundary;
}

static void find_far_patch(blockmesh_pair_t* pair,
                           int i1, int j1, 
                           int i2_max, int j2_max, 
                           int* i2, int* j2)
{
  if (pair->rotation == NO_ROTATION)
  {
    *i2 = i1;
    *j2 = j1;
  }
  else if (pair->rotation == QUARTER_TURN)
  {
  }
  else if (pair->rotation == HALF_TURN)
  {
  }
  else // pair->rotation == THREE_QUARTERS_TURN
  {
  }
}

void blockmesh_pair_find_patch(blockmesh_pair_t* pair,
                               int i1, int j1, int k1,
                               int* i2, int* j2, int* k2)
{
  int b1 = pair->block1_boundary;
  int b2 = pair->block2_boundary;
  int npx1 = pair->npx1, npy1 = pair->npy1, npz1 = pair->npz1;
//  int npx2 = pair->npx2, npy2 = pair->npy2, npz2 = pair->npz2;
  if ((b1 == 0) && (i1 == 0))
  {
    if (b2 == 0)      // -x <-> -x connection
    {
      *i2 = 0;
      find_far_patch(pair, j1, k1, npy1, npz1, j2, k2);
    }
    else if (b2 == 1) // -x <-> +x connection
    {
      *i2 = npx1-1;
      find_far_patch(pair, j1, k1, npy1, npz1, j2, k2);
    }
    else if (b2 == 2) // -x <-> -y connection
    {
    }
    else if (b2 == 3) // -x <-> +y connection
    {
    }
    else if (b2 == 4) // -x <-> -z connection
    {
    }
    else if (b2 == 5) // -x <-> +z connection
    {
    }
  }
  else if ((b1 == 1) && (i1 == (npx1 - 1))) 
  {
    if (b2 == 0)      // +x <-> -x connection
    {
      *i2 = 0;
      find_far_patch(pair, j1, k1, npy1, npz1, j2, k2);
    }
    else if (b2 == 1) // +x <-> +x connection
    {
      *i2 = npx1-1;
      find_far_patch(pair, j1, k1, npy1, npz1, j2, k2);
    }
    else if (b2 == 2) // +x <-> -y connection
    {
    }
    else if (b2 == 3) // +x <-> +y connection
    {
    }
    else if (b2 == 4) // +x <-> -z connection
    {
    }
    else if (b2 == 5) // +x <-> +z connection
    {
    }
  }
  else if ((b1 == 2) && (j1 == 0))
  {
    if (b2 == 0)      // -y <-> -x connection
    {
    }
    else if (b2 == 1) // -y <-> +x connection
    {
    }
    else if (b2 == 2) // -y <-> -y connection
    {
      *j2 = 0;
      find_far_patch(pair, k1, i1, npz1, npx1, k2, i2);
    }
    else if (b2 == 3) // -y <-> +y connection
    {
      *j2 = npy1-1;
      find_far_patch(pair, k1, i1, npz1, npx1, k2, i2);
    }
    else if (b2 == 4) // -y <-> -z connection
    {
    }
    else if (b2 == 5) // -y <-> +z connection
    {
    }
  }
  else if ((b1 == 3) && (j1 == (npy1 - 1))) 
  {
    if (b2 == 0)      // +y <-> -x connection
    {
    }
    else if (b2 == 1) // +y <-> +x connection
    {
    }
    else if (b2 == 2) // +y <-> -y connection
    {
      *j2 = 0;
      find_far_patch(pair, k1, i1, npz1, npx1, k2, i2);
    }
    else if (b2 == 3) // +y <-> +y connection
    {
      *j2 = npy1-1;
      find_far_patch(pair, k1, i1, npz1, npx1, k2, i2);
    }
    else if (b2 == 4) // +y <-> -z connection
    {
    }
    else if (b2 == 5) // +y <-> +z connection
    {
    }
  }
  else if ((b1 == 4) && (k1 == 0))
  {
    if (b2 == 0)      // -z <-> -x connection
    {
    }
    else if (b2 == 1) // -z <-> +x connection
    {
    }
    else if (b2 == 2) // -z <-> -y connection
    {
    }
    else if (b2 == 3) // -z <-> +y connection
    {
    }
    else if (b2 == 4) // -z <-> -z connection
    {
      *k2 = 0;
      find_far_patch(pair, i1, j1, npx1, npy1, i2, j2);
    }
    else if (b2 == 5) // -z <-> +z connection
    {
      *k2 = npz1-1;
      find_far_patch(pair, i1, j1, npx1, npy1, i2, j2);
    }
  }
  else if ((b1 == 5) && (k1 == (npz1 - 1))) 
  {
    if (b2 == 0)      // +z <-> -x connection
    {
    }
    else if (b2 == 1) // +z <-> +x connection
    {
    }
    else if (b2 == 2) // +z <-> -y connection
    {
    }
    else if (b2 == 3) // +z <-> +y connection
    {
    }
    else if (b2 == 4) // +z <-> -z connection
    {
      *k2 = 0;
      find_far_patch(pair, i1, j1, npx1, npy1, i2, j2);
    }
    else if (b2 == 5) // +z <-> +z connection
    {
      *k2 = npz1-1;
      find_far_patch(pair, i1, j1, npx1, npy1, i2, j2);
    }
  }
}

size_t blockmesh_pair_data_size(blockmesh_pair_t* pair,
                                unimesh_centering_t centering,
                                int num_comp)
{
  return 0; // FIXME
}

void blockmesh_pair_copy_in(blockmesh_pair_t* pair,
                            int i, int j, int k,
                            unimesh_patch_t* source_patch,
                            void* buffer)
{
}

void blockmesh_pair_copy_out(blockmesh_pair_t* pair,
                             void* buffer,
                             int i, int j, int k,
                             unimesh_patch_t* dest_patch)
{
}

