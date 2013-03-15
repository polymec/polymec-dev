#include <gc/gc.h>
#include "geometry/block.h"
#include "geometry/hexahedron.h"
#include "core/unordered_set.h"

#ifdef __cplusplus
extern "C" {
#endif

struct block_t 
{
  hexahedron_t* hex;
  str_unordered_set_t* tags[6];
};

static void block_free(void* ctx, void* dummy)
{
  block_t* block = ctx;
  block->hex = NULL;
  for (int f = 0; f < 6; ++f)
    str_unordered_set_free(block->tags[f]);
}


block_t* block_new(int order)
{
  ASSERT(order >= 1);
  ASSERT(order <= 3);
  block_t* block = GC_MALLOC(sizeof(block_t));
  block->hex = hexahedron_new(order);
  for (int f = 0; f < 6; ++f)
    block->tags[f] = str_unordered_set_new();

  GC_register_finalizer(block, block_free, block, NULL, NULL);
  return block;
}

int block_order(block_t* block)
{
  return hexahedron_order(block->hex);
}

int block_num_points(block_t* block)
{
  return hexahedron_num_points(block->hex);
}

void block_get_points(block_t* block, point_t* points)
{
  hexahedron_get_points(block->hex, points);
}

void block_map(block_t* block, point_t* points, point_t* xi, point_t* x)
{
  hexahedron_map(block->hex, points, xi, x);
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

void block_delete_tag(block_t* block, int face_index, const char* tag)
{
  ASSERT(face_index >= 0);
  ASSERT(face_index < 6);
  str_unordered_set_delete(block->tags[face_index], (char*)tag);
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
    assembly->blocks[i] = NULL;
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

