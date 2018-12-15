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
///
/// Blocks in a block mesh are stitched together by identifying nodes on the 
/// common face between two blocks. Each block has 8 nodes. Looking down 
/// (in the -z direction) on a block in the xy plane, nodes 0-3 traverse the 
/// bottom block face counterclockwise, starting with the "lower left" node.
/// Nodes 4-7 traverse the top block face in the same way. This is the standard
/// way that nodes are indexed in hexahedral elements in the finite element 
/// method.
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
/// \param [in] block A \ref unimesh presenting a block. This unimesh must not 
///                   yet be finalized, since boundary conditions may be added
///                   to it when it is connected to other blocks.
/// \returns the index of the new block.
/// \memberof blockmesh
int blockmesh_add_block(blockmesh_t* mesh, unimesh_t* block);

/// Returns the face associated with the given set of block nodes, or -1 
/// if the nodes do not match any of the block's faces. 
/// \param [in] block_nodes An array of 4 block nodes that supposedly match those
///                         in one of the block's 6 faces. The order in which 
///                         the nodes are specified doesn't matter.
int blockmesh_face_for_nodes(blockmesh_t* mesh, int block_nodes[4]);

/// Connects two blocks with the given indices within a block mesh in a manner 
/// specified by parameters. You must call this function on every process 
/// within the mesh's communicator, and the blocks must be connected in such a 
/// way that all block dimensions are compatible. The two indices can refer to 
/// the same block, in which case the block connects to itself.
/// \param [in] block1_index The index of the first of the two blocks to be 
///                          connected within the mesh.
/// \param [in] block1_nodes An array containing the 4 nodes in the first block
///                          to be identified with the corresponding nodes 
///                          in the second block (block2_nodes).
/// \param [in] block2_index The index of the second of the two blocks to be 
///                          connected within the mesh.
/// \param [in] block2_nodes An array containing the 4 nodes in the second block
///                          to be identified with the corresponding nodes 
///                          in the first block (block1_nodes).
/// \memberof blockmesh
/// \collective
void blockmesh_connect_blocks(blockmesh_t* mesh, 
                              int block1_index, 
                              int block1_nodes[4],
                              int block2_index,
                              int block2_nodes[4]);

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

