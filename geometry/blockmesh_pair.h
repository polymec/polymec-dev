// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_BLOCKMESH_PAIR_H
#define POLYMEC_BLOCKMESH_PAIR_H

/// \addtogroup geometry geometry
///@{

/// \struct blockmesh_diffeomorphism_t
/// This struct represents a diffeomorphism that maps quantities from one 
/// block (block1) to another (block2) in a blockmesh. The diffeomorphism is 
/// defined by the coordinate systems for the two blocks and a rotation 
/// in the plane of the boundary connecting them.
typedef struct
{
  /// The coordinate mapping for the first block.
  coord_mapping_t* block1_coords;

  /// The coordinate mapping for the second block.
  coord_mapping_t* block2_coords;

  /// The rotation that must be performed for a quantity to pass through the
  /// face connecting the two blocks. Expressed in units of counterclockwise
  /// revolutions (turns).
  enum
  {
    NO_ROTATION,
    QUARTER_TURN,
    HALF_TURN,
    THREE_QUARTERS_TURN,
    INVALID_ROTATION
  } rotation;
} blockmesh_diffeomorphism_t;

/// \class blockmesh_pair
/// This class represents the pairing of two blocks within a blockmesh.
/// \refcounted
typedef struct blockmesh_pair_t blockmesh_pair_t;

/// Creates a new blockmesh pair for two blocks in the given blockmesh.
/// \param [in] mesh The blockmesh for the newly-created pair.
/// \param [in] block1 The index of the first block in the pair (within the mesh).
/// \param [in] boundary1 The boundary of the first block that connects to the second.
/// \param [in] block2 The index of the second block in the pair (within the mesh).
/// \param [in] boundary2 The boundary of the second block that connects to the first.
/// \memberof blockmesh_pair
blockmesh_pair_t* blockmesh_pair_new(blockmesh_t* mesh, 
                                     int block1, unimesh_boundary_t boundary1,
                                     int block2, unimesh_boundary_t boundary2);

/// Returns an internal pointer to a diffeomorphism that maps quantities from the first 
/// block to the second block in this pair. 
/// \memberof blockmesh_pair
blockmesh_diffeomorphism_t* blockmesh_pair_diffeomorphism(blockmesh_pair_t* pair);

/// Determines the patch (i2, j2, k2) within the second block in the pair
/// that corresponds to the patch (i1, j1, k1) in the first block, given the
/// specifics of the connectivity.
/// \param [in] i1 The first logical coordinate for a patch in the first block.
/// \param [in] j1 The second logical coordinate for a patch in the first block.
/// \param [in] k1 The third logical coordinate for a patch in the first block.
/// \param [out] i2 Stores the first logical coordinate for the patch in the second block.
/// \param [out] j2 Stores the second logical coordinate for the patch in the second block.
/// \param [out] k2 Stores the third logical coordinate for the patch in the second block.
/// \memberof blockmesh_pair
void blockmesh_pair_find_patch(blockmesh_pair_t* pair,
                               int i1, int j1, int k1,
                               int* i2, int* j2, int* k2);

///@}

#endif

