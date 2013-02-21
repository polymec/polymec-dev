#include "geometry/block.h"
#include "core/unordered_set.h"

#ifdef __cplusplus
extern "C" {
#endif

struct block_t 
{
  int order, num_points;
  int* point_indices;
  block_t* neighbors[6];
  str_unordered_set_t* tags[6];
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
  return block;
}

static void block_free(block_t* block)
{
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

void block_map(block_t* block, point_t* xi, point_t* x)
{
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

