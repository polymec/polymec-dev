#ifndef POLYMEC_BLOCK_H
#define POLYMEC_BLOCK_H

#include "core/point.h"

#ifdef __cplusplus
extern "C" {
#endif

// This file contains data structures that aid in constructing block-structured
// meshes.

// A block is a region of space that is topologically a hexahedron, with 
// nodes along its faces that define an isoparametric mapping. Objects of this
// type can be created independently, or retrieved from a block assembly and 
// queried. They are garbage-collected.
//
// Block faces are indexed as follows:
//  0 - negative 'x' face
//  1 - positive 'x' face
//  2 - negative 'y' face
//  3 - positive 'y' face
//  4 - negative 'z' face
//  5 - positive 'z' face
// Note that the coordinates are in quotes because 'x', 'y', and 'z' are 
// topological, not geometric coordinates. In other words, 'x' may not be 
// aligned with x in physical space, etc.
//
// It is perhaps easiest to think of a block as being equivalent to a 
// hexahedral finite element of a given geometric order (linear, quadratic, 
// etc).
typedef struct block_t block_t;

// Creates a block of the given order.
block_t* block_new(int order);

// Returns the geometric order of the block.
int block_order(block_t* block);

// Returns the number of points that define the geometry/topology of the 
// block. This is related to the order of the block, but provides a shortcut.
int block_num_points(block_t* block);

// Retrieves the logical coordinates of the points within the block.
void block_get_points(block_t* block, point_t* points);

// Maps a point xi within the logical space of the block to its 
// equivalent point x in physical space, using the isoparametric mapping 
// for the block defined by the given points.
void block_map(block_t* block, point_t* points, point_t* xi, point_t* x);

// Adds a tag to the block's given face.
void block_add_tag(block_t* block, int face_index, const char* tag);

// Delete the given tag from the block's given face. If the tag doesn't 
// exist there, this function has no effect.
void block_delete_tag(block_t* block, int face_index, const char* tag);

// Traverses the set of tags attached to the given face of the block.
bool block_next_tag(block_t* block, int face_index, int* pos, char** tag);

// A block assembly is a collection of blocks, coupled with a set of points 
// that define their positions in space and their topology. This is akin to 
// a hexahedral block of finite elements in the sense that a block is akin to 
// a hexahedral finite element.
typedef struct block_assembly_t block_assembly_t;

// Allocates a new block assembly with blocks of the given geometric order, 
// and a given number of blocks and points.
block_assembly_t* block_assembly_new(int order, int num_blocks, int num_points);

// Frees the block assembly.
void block_assembly_free(block_assembly_t* assembly);

// Returns the geometric order of blocks within the block assembly.
int block_assembly_order(block_assembly_t* assembly);

// Returns the number of blocks within the block assembly.
int block_assembly_num_blocks(block_assembly_t* assembly);

// Returns the ith block within the block assembly.
block_t* block_assembly_block(block_assembly_t* assembly, int i);

// Returns the number of points within the block assembly.
int block_assembly_num_points(block_assembly_t* assembly);

// Returns the points within the block assembly.
point_t* block_assembly_points(block_assembly_t* assembly);

// A block grid is an object that specifies a grid resolution within a block.
// A set of block grids plus a block assembly defines a mesh. N1, N2, and N3
// are respectively the number of cells on a side for the primary, secondary, 
// and tertiary directions of the block (as defined by its points).
typedef struct 
{
  int N1, N2, N3;
} block_grid_t;

#ifdef __cplusplus
}
#endif

#endif

