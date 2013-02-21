#include "geometry/block.h"
#include "core/unordered_set.h"
#include "core/lagrange_poly.h"

#ifdef __cplusplus
extern "C" {
#endif

struct block_t 
{
  int order, num_points;
  int* point_indices;
  block_t* neighbors[6];
  str_unordered_set_t* tags[6];
  tensor_lagrange_poly_t* poly;
};

static block_t* block_new(int order)
{
  ASSERT(order >= 1);
  ASSERT(order >= 3);
  block_t* block = malloc(sizeof(block_t));
  block->order = order;
  block->num_points = pow(order+1, 3);
  block->point_indices = malloc(sizeof(int)*block->num_points);
  for (int f = 0; f < 6; ++f)
  {
    block->neighbors[f] = NULL;
    block->tags[f] = str_unordered_set_new();
  }
  block->poly = tensor_lagrange_poly_new(order);

  // Set the points of the block to those of the reference (logical) 
  // coordinate system.
  double points[order+1];
  double h = 1.0 / order;
  for (int i = 0; i < order+1; ++i)
    points[i] = i*h;
  tensor_lagrange_poly_set_points(block->poly, points, points, points);
  return block;
}

static void block_free(block_t* block)
{
  free(block->poly);
  for (int f = 0; f < 6; ++f)
    str_unordered_set_free(block->tags[f]);
  free(block->point_indices);
  free(block);
}

int block_order(block_t* block)
{
  return block->order;
}

int block_num_points(block_t* block)
{
  return block->num_points;
}

int* block_point_indices(block_t* block)
{
  return block->point_indices;
}

block_t* block_neighbor(block_t* block, int face_index)
{
  ASSERT(face_index >= 0);
  ASSERT(face_index < 6);
  return block->neighbors[face_index];
}

void block_map(block_t* block, point_t* points, point_t* xi, point_t* x)
{
  int order = block->order;
  double xs[order+1], ys[order+1], zs[order+1];
  for (int i = 0; i < order+1; ++i)
  {
    xs[i] = points[i].x;
    ys[i] = points[i].y;
    zs[i] = points[i].z;
  }
  x->x = tensor_lagrange_poly_value(block->poly, xi, xs);
  x->y = tensor_lagrange_poly_value(block->poly, xi, ys);
  x->z = tensor_lagrange_poly_value(block->poly, xi, zs);
}

static void tag_free(char* tag)
{
  free(tag);
}

void block_add_tag(block_t* block, int face_index, const char* tag)
{
  ASSERT(face_index >= 0);
  ASSERT(face_index < 6);
  str_unordered_set_insert_with_dtor(block->tags[face_index], strdup(tag), tag_free);
}

bool block_next_tag(block_t* block, int face_index, int* pos, char** tag)
{
  ASSERT(face_index >= 0);
  ASSERT(face_index < 6);
  return str_unordered_set_next(block->tags[face_index], pos, tag);
}

struct block_assembly_t 
{
  int order;
  int num_blocks;
  block_t** blocks;
  int num_points;
  point_t* points;
};

block_assembly_t* block_assembly_new(int order, int num_blocks, int num_points)
{
  ASSERT(order >= 1);
  ASSERT(order <= 3);
  ASSERT(num_blocks >= 1);
  ASSERT(num_points >= 8);
  block_assembly_t* assembly = malloc(sizeof(block_assembly_t));
  assembly->order = order;
  assembly->num_blocks = num_blocks;
  assembly->blocks = malloc(sizeof(block_t)*num_blocks);
  for (int i = 0; i < num_blocks; ++i)
    assembly->blocks[i] = block_new(order);
  assembly->num_points = num_points;
  assembly->points = malloc(sizeof(point_t)*num_points);
  return assembly;
}

void block_assembly_free(block_assembly_t* assembly)
{
  for (int i = 0; i < assembly->num_blocks; ++i)
    block_free(assembly->blocks[i]);
  free(assembly->blocks);
  free(assembly->points);
  free(assembly);
}

int block_assembly_order(block_assembly_t* assembly)
{
  return assembly->order;
}

int block_assembly_num_blocks(block_assembly_t* assembly)
{
  return assembly->num_blocks;
}

block_t* block_assembly_block(block_assembly_t* assembly, int i)
{
  ASSERT(i >= 0);
  ASSERT(i < assembly->num_blocks);
  return assembly->blocks[i];
}

int block_assembly_num_points(block_assembly_t* assembly)
{
  return assembly->num_points;
}

point_t* block_assembly_points(block_assembly_t* assembly)
{
  return assembly->points;
}

#ifdef __cplusplus
}
#endif

