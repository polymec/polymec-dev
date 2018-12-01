// Copyright (c) 2012-2018, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_BLOCKMESH_H
#define POLYMEC_BLOCKMESH_H

#include "geometry/unimesh.h"

/// \addtogroup geometry geometry
///@{

/// \class blockmesh
/// A blockmesh, or block mesh, is a collection of three-dimensional cartesian 
/// meshes (unimeshes) called blocks. These blocks can be connected to one 
/// another as long as their boundaries have compatible dimensions (numbers of
/// cells). The block mesh manages these blocks and their connectivity.
///
/// A blockmesh is distributed across processes in an MPI communicator. Every 
/// process contains all of the blocks in a blockmesh. Each of these blocks
/// is domain-decomposed across these processes.
typedef struct blockmesh_t blockmesh_t;

//------------------------------------------------------------------------
//                          Construction methods
//------------------------------------------------------------------------
// The following methods are used to construct block meshs.
// blockmesh_finalize() must be called after a mesh has been properly
// constructed.
//------------------------------------------------------------------------

/// Creates a new blockmesh with no blocks, defined on the given communicator.
/// \memberof blockmesh
blockmesh_t* blockmesh_new(MPI_Comm comm);

/// Adds a new block to the blockmesh. 
/// \param [in] block A \ref unimesh presenting a block.
/// \returns the index of the new block.
/// \memberof blockmesh
int blockmesh_add_block(blockmesh_t* mesh, unimesh_t* block);

/// \enum blockmesh_cxn_t
/// Describes the way that two blocks are connected within a block mesh.
typedef enum
{
  /// Connects a block to another block without rotating. 
  BLOCKMESH_CXN_UNROTATED = 0,
  /// Connects a block to another block, rotating counterclockwise by 1/4. 
  BLOCKMESH_CXN_QUARTER_TURN = 1,
  /// Connects a block to another block, rotating counterclockwise by 1/2. 
  BLOCKMESH_CXN_HALF_TURN = 2,
  /// Connects a block to another block, rotating counterclockwise by 3/4. 
  BLOCKMESH_CXN_THREE_QUARTER_TURN = 3
} blockmesh_cxn_t;

/// Connects two blocks with the given indices within a block mesh in a manner 
/// specified by parameters. Must be called on every process within the 
/// mesh's communicator. The blocks must be connected in such a way that all 
/// block dimensions are compatible. The two indices can refer to the same 
/// block, in which case the block connects to itself.
/// \param [in] index1 The index of the first of the two blocks to be connected
///                    within the mesh.
/// \param [in] boundary1 The boundary of the first block to connect to the 
///                       second.
/// \param [in] index2 The index of the second of the two blocks to be connected
///                    within the mesh.
/// \param [in] boundary2 The boundary of the second_block to connect to the 
///                       second.
/// \param [in] connection The type of connection between the two blocks.
/// \memberof blockmesh
void blockmesh_connect_blocks(blockmesh_t* mesh, 
                              int index1, unimesh_boundary_t boundary1,
                              int index2, unimesh_boundary_t boundary2,
                              blockmesh_cxn_t connection);

/// Finalizes the construction process for the block mesh. This must be called 
/// before any of the mesh's usage methods are invoked. Should only 
/// be called once.
/// \memberof blockmesh
void blockmesh_finalize(blockmesh_t* mesh);

//------------------------------------------------------------------------
//                          Usage methods
//------------------------------------------------------------------------
// The following methods can only be used after a blockmesh has been 
// fully constructed and finalized.
//------------------------------------------------------------------------

/// Returns true if the given mesh has been finalized, false otherwise.
/// \memberof blockmesh
bool blockmesh_is_finalized(blockmesh_t* mesh);

/// Destroys the given mesh and all of its patches.
/// \memberof blockmesh
void blockmesh_free(blockmesh_t* mesh);

/// Returns the MPI communicator on which the blockmesh is defined.
/// \memberof blockmesh
MPI_Comm blockmesh_comm(blockmesh_t* mesh);

/// Returns the number of blocks within the mesh.
/// \memberof blockmesh
size_t blockmesh_num_blocks(blockmesh_t* mesh);

/// Returns the block with the given index within the mesh.
/// \param [in] index The index of the requested block.
/// \memberof blockmesh
unimesh_t* blockmesh_block(blockmesh_t* mesh, int index);

typedef struct blockmesh_field_t blockmesh_field_t;

/// Repartitions the given blockmesh and redistributes data to each of the 
/// given fields. Here, the old meshes and fields are consumed, and new ones 
/// are created in their place. Weights can be provided for each patch, and 
/// the partitioning is performed so that the load imbalance does not exceed 
/// the given tolerance.
/// \note In addition, each repartitioned field needs to have any boundary 
/// conditions reinstated, since these boundary conditions are not 
/// transmitted between processes.
/// \relates blockmesh
/// \collective Collective on mesh's communicator.
void repartition_blockmesh(blockmesh_t** mesh, 
                           int* weights,
                           real_t imbalance_tol,
                           blockmesh_field_t** fields,
                           size_t num_fields);

///@}

#endif

