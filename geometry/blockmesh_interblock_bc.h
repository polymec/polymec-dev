// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_BLOCKMESH_INTERBLOCK_BC_H
#define POLYMEC_BLOCKMESH_INTERBLOCK_BC_H

#include "geometry/unimesh.h"

/// \addtogroup geometry geometry
///@{

typedef struct blockmesh_t blockmesh_t;

/// \class blockmesh_interblock_bc
/// This patch BC connects patches in different unimeshes (blocks) within
/// a blockmesh, mapping quantities between these patches with a
/// diffeomorphism defined by the respective coordinate systems of the blocks.
/// \experimental
typedef struct blockmesh_interblock_bc_t blockmesh_interblock_bc_t;

/// Constructs a new inter-block BC for the given mesh (block).
/// \param [in] mesh The block mesh on which this BC operates.
/// \memberof blockmesh_interblock_bc
blockmesh_interblock_bc_t* blockmesh_interblock_bc_new(blockmesh_t* mesh);

/// Destroys the given inter-block BC.
void blockmesh_interblock_bc_free(blockmesh_interblock_bc_t* bc);

/// Establishes a connection between a patch in the block associated with
/// this BC and another block.
/// \param [in] block1_index The index of the first block within the mesh.
/// \param [in] block1_boundary The boundary of the first block along which
///                             the patch is located.
/// \param [in] i1 The i index identifying a patch in the first block.
/// \param [in] j1 The j index identifying a patch in the first block.
/// \param [in] k1 The k index identifying a patch in the first block.
/// \param [in] rotation The number of topological counterclockwise quarter
///                      turns experienced by a body that moves through the
///                      block-block interface, from the first block to the
///                      second.
/// \param [in] block2_index The index of the second block within the mesh.
/// \param [in] block2_boundary The boundary of the second block along which
///                             the patch is located.
/// \param [in] i2 The i index identifying a patch in the second block.
/// \param [in] j2 The j index identifying a patch in the second block
/// \param [in] k2 The k index identifying a patch in the second block.
/// \memberof blockmesh_interblock_bc
void blockmesh_interblock_bc_connect(blockmesh_interblock_bc_t* bc,
                                     int block1_index,
                                     unimesh_boundary_t block1_boundary,
                                     int i1, int j1, int k1,
                                     int rotation,
                                     int block2_index,
                                     unimesh_boundary_t block2_boundary,
                                     int i2, int j2, int k2);

/// Finalizes the BC once all connections have been established.
/// \memberof blockmesh_interblock_bc
void blockmesh_interblock_bc_finalize(blockmesh_interblock_bc_t* bc);

///@}

#endif

