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

/// \addtogroup geometry geometry
///@{

typedef struct blockmesh_t blockmesh_t;
typedef struct blockmesh_pair_t blockmesh_pair_t;

/// \class blockmesh_interblock_bc
/// This patch BC connects patches in different unimeshes (blocks) within 
/// a blockmesh, mapping quantities between these patches with a 
/// diffeomorphism defined by the respective coordinate systems of the blocks.
typedef struct blockmesh_interblock_bc_t blockmesh_interblock_bc_t;

/// Constructs a new inter-block BC for the given mesh (block).
/// \param [in] mesh The block mesh on which this BC operates.
/// \memberof blockmesh_interblock_bc
blockmesh_interblock_bc_t* blockmesh_interblock_bc_new(blockmesh_t* mesh);

/// Destroys the given inter-block BC.
void blockmesh_interblock_bc_free(blockmesh_interblock_bc_t* bc);

/// Establishes a connection between a patch in the block associated with 
/// this BC and another block.
/// \param [in] block_pair The object representing the paired blocks.
/// \param [in] i1 The i index identifying a patch in the first block.
/// \param [in] j1 The j index identifying a patch in the first block.
/// \param [in] k1 The k index identifying a patch in the first block.
/// \param [in] i2 The i index identifying a patch in the second block.
/// \param [in] j2 The j index identifying a patch in the second block
/// \param [in] k2 The k index identifying a patch in the second block.
/// \memberof blockmesh_interblock_bc
void blockmesh_interblock_bc_connect(blockmesh_interblock_bc_t* bc,
                                     blockmesh_pair_t* block_pair,
                                     int i1, int j1, int k1, 
                                     int i2, int j2, int k2);

/// Finalizes the BC once all connections have been established.
/// \memberof blockmesh_interblock_bc
void blockmesh_interblock_bc_finalize(blockmesh_interblock_bc_t* bc);

/// Traverses the connections in the blockmesh associated with this BC.
/// \param [out] block_pair Stores the block pair object representing the 
///                         pairing of block1 with block2.
/// \param [out] i1 Stores the i index of the patch in the first block.
/// \param [out] j1 Stores the j index of the patch in the first block.
/// \param [out] k1 Stores the k index of the patch in the first block.
/// \param [out] i2 Stores the i index of the patch in the first block.
/// \param [out] j2 Stores the j index of the patch in the first block.
/// \param [out] k2 Stores the k index of the patch in the first block.
/// \memberof blockmesh_interblock_bc
bool blockmesh_interblock_bc_next_connection(blockmesh_interblock_bc_t* bc,
                                             int* pos,
                                             blockmesh_pair_t** block_pair,
                                             int* i1, int* j1, int* k1,
                                             int* i2, int* j2, int* k2);

///@}

#endif

