// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_BLOCKMESH_INTERBLOCK_BC_H
#define POLYMEC_BLOCKMESH_INTERBLOCK_BC_H

#include "geometry/coord_mapping.h"
#include "geometry/unimesh.h"
#include "geometry/unimesh_patch_bc.h"

/// \addtogroup geometry geometry
///@{

typedef struct blockmesh_t blockmesh_t;

/// \class blockmesh_interblock_bc
/// This patch BC connects patches in different unimeshes (blocks) within 
/// a blockmesh, mapping quantities between these patches with a 
/// diffeomorphism defined by the respective coordinate systems of the blocks.

/// \struct blockmesh_diffeomorphism_t
/// This struct represents a diffeomorphism that maps quantities from one 
/// block (block1) to another (block2) in a blockmesh. The diffeomorphism is 
/// defined by the two coordinate systems for the two blocks and a rotation 
/// in the plane of the face connecting them.
typedef struct
{
  /// The coordinate system for the first block, defined by a mapping from 
  /// the block's logical coordinates [0,1] x [0,1] x [0,1].
  coord_mapping_t* block1_coords;

  /// The coordinate system for the second block, defined by a mapping from 
  /// the block's logical coordinates [0,1] x [0,1] x [0,1].
  coord_mapping_t* block2_coords;

  /// The rotation that must be performed for a quantity to pass through the
  /// face connecting the two blocks. Expressed in units of counterclockwise
  /// revolutions (turns).
  enum
  {
    NO_ROTATION,
    QUARTER_TURN,
    HALF_TURN,
    THREE_QUARTERS_TURN
  } rotation;
} blockmesh_diffeomorphism_t;

/// Constructs a new inter-block BC for the given mesh (block).
/// \param [in] mesh The block mesh on which this BC operates.
/// \memberof blockmesh_interblock_bc
unimesh_patch_bc_t* blockmesh_interblock_bc_new(blockmesh_t* mesh);

/// Establishes a connection between a patch in the block associated with 
/// this BC and another block.
/// \param [in] block1 The first block in the pair.
/// \param [in] i1 The i index identifying a patch in block1.
/// \param [in] j1 The j index identifying a patch in block1.
/// \param [in] k1 The k index identifying a patch in block1.
/// \param [in] b1 The boundary of the patch (i1, j1, k1) within block1 to be connected to block2.
/// \param [in] block2 The second block in the pair.
/// \param [in] i2 The i index identifying a patch in block2.
/// \param [in] j2 The j index identifying a patch in block2
/// \param [in] k2 The k index identifying a patch in block2.
/// \param [in] b2 The boundary of the patch (i2, j2, k2) within block2 to be connected to block1.
/// \param [in] diff A diffeomorphism defining the mapping of quantities 
///                  from block1 to block2.
/// \memberof blockmesh_interblock_bc
void blockmesh_interblock_bc_connect(unimesh_patch_bc_t* bc,
                                     unimesh_t* block1,
                                     int i1, int j1, int k1, 
                                     unimesh_boundary_t b1,
                                     unimesh_t* block2,
                                     int i2, int j2, int k2, 
                                     unimesh_boundary_t b2,
                                     blockmesh_diffeomorphism_t diff);

/// Finalizes the BC once all connections have been established.
/// \memberof blockmesh_interblock_bc
void blockmesh_interblock_bc_finalize(unimesh_patch_bc_t* bc);

/// Traverses the connections in the blockmesh associated with this BC.
/// \param [out] block1 Stores the first block in the connected pair.
/// \param [out] i1 Stores the i index of the patch in the first block.
/// \param [out] j1 Stores the j index of the patch in the first block.
/// \param [out] k1 Stores the k index of the patch in the first block.
/// \param [out] b1 Stores the boundary of the patch (i1, j1, k1) within block1 connected to block2.
/// \param [out] block2 Stores the second block in the connected pair.
/// \param [out] i2 Stores the i index of the patch in the first block.
/// \param [out] j2 Stores the j index of the patch in the first block.
/// \param [out] k2 Stores the k index of the patch in the first block.
/// \param [out] b2 Stores the boundary of the patch (i2, j2, k2) within block2 connected to block1.
/// \param [out] diff Stores the diffeomorphism defining the mapping of 
///                   quantities from block1 to block2.
/// \memberof blockmesh_interblock_bc
bool blockmesh_interblock_bc_next_connection(unimesh_patch_bc_t* bc,
                                             int* pos,
                                             unimesh_t** block1, 
                                             int* i1, int* j1, int* k1,
                                             unimesh_boundary_t* b1,
                                             unimesh_t** block2, 
                                             int* i2, int* j2, int* k2,
                                             unimesh_boundary_t* b2,
                                             blockmesh_diffeomorphism_t* diff);

///@}

#endif

