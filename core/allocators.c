// Copyright (c) 2012-2014, Jeffrey N. Johnson
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without 
// modification, are permitted provided that the following conditions are met:
// 
// 1. Redistributions of source code must retain the above copyright notice, this 
// list of conditions and the following disclaimer.
// 
// 2. Redistributions in binary form must reproduce the above copyright notice, 
// this list of conditions and the following disclaimer in the documentation 
// and/or other materials provided with the distribution.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE 
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include "core/allocators.h"
#include "core/slist.h"
#include "arena/proto.h"
#include "arena/pool.h"

struct polymec_allocator_t 
{
  char* name;
  void* context;
  polymec_allocator_vtable vtable;
};

polymec_allocator_t* polymec_allocator_new(const char* name,
                                           void* context,
                                           polymec_allocator_vtable vtable)
{
  ASSERT(vtable.malloc != NULL);
  ASSERT(vtable.realloc != NULL);
  ASSERT(vtable.free != NULL);

  polymec_allocator_t* alloc = malloc(sizeof(polymec_allocator_t)); // Oh, the irony...
  alloc->name = malloc(sizeof(char) * strlen(name) + 1); // Need to avoid using other allocators.
  memcpy(alloc->name, name, strlen(name));
  alloc->context = context;
  alloc->vtable = vtable;
  return alloc;
}

void polymec_allocator_free(polymec_allocator_t* alloc)
{
  if ((alloc->context != NULL) && (alloc->vtable.dtor != NULL))
    alloc->vtable.dtor(alloc->context);
  free(alloc);
}

// Allocator stack.
static ptr_slist_t* alloc_stack = NULL;

static void free_alloc_stack()
{
  ASSERT(alloc_stack != NULL);
  ptr_slist_free(alloc_stack);
}

void push_allocator(polymec_allocator_t* alloc)
{
  if (alloc_stack == NULL)
  {
    alloc_stack = ptr_slist_new();
    polymec_atexit(free_alloc_stack);
  }
  ptr_slist_push_with_dtor(alloc_stack, alloc, DTOR(polymec_allocator_free));
}

polymec_allocator_t* pop_allocator()
{
  ASSERT(alloc_stack != NULL);
  return ptr_slist_pop(alloc_stack, NULL);
}

void* polymec_malloc(size_t size)
{
  if ((alloc_stack == NULL) || (alloc_stack->size == 0))
    return malloc(size);
  else
  {
    polymec_allocator_t* alloc = alloc_stack->front->value;
    return alloc->vtable.malloc(alloc->context, size);
  }
}

void* polymec_aligned_alloc(size_t alignment, size_t size)
{
  if ((alloc_stack == NULL) || (alloc_stack->size == 0))
#ifdef APPLE
  {
    // Mac OS X doesn't have this yet!
    polymec_error("Mac OS standard allocator does not support aligned allocation.");
  }
#else
    return aligned_alloc(alignment, size);
#endif
  else
  {
    polymec_allocator_t* alloc = alloc_stack->front->value;
    if (alloc->vtable.aligned_alloc != NULL)
      return alloc->vtable.aligned_alloc(alloc->context, alignment, size);
    else
      polymec_error("%s allocator does not support aligned allocation.");
  }
}

void* polymec_realloc(void* memory, size_t size)
{
  if ((alloc_stack == NULL) || (alloc_stack->size == 0))
    return realloc(memory, size);
  else
  {
    polymec_allocator_t* alloc = alloc_stack->front->value;
    return alloc->vtable.realloc(alloc->context, memory, size);
  }
}

void* polymec_aligned_realloc(void* memory, size_t alignment, size_t size)
{
  if ((alloc_stack == NULL) || (alloc_stack->size == 0))
  {
    // Not possible in general, since we can't get the old memory size.
    polymec_error("No allocator is available, so aligned realloc() is not possible.");
  }
  else
  {
    polymec_allocator_t* alloc = alloc_stack->front->value;
    if (alloc->vtable.aligned_realloc != NULL)
      return alloc->vtable.aligned_realloc(alloc->context, memory, alignment, size);
    else
      polymec_error("%s allocator does not support aligned allocation.");
  }
}

void polymec_free(void* memory)
{
  if ((alloc_stack == NULL) || (alloc_stack->size == 0))
    free(memory);
  else 
  {
    polymec_allocator_t* alloc = alloc_stack->front->value;
    alloc->vtable.free(alloc->context, memory);
  }
}

static void* std_malloc(void* context, size_t size)
{
  return malloc(size);
}

static void* std_aligned_alloc(void* context, size_t alignment, size_t size)
{
#ifdef APPLE
  polymec_error("Aligned allocations are not available on Mac OS at this time.");
#else
  return aligned_alloc(alignment, size);
#endif
}

static void* std_realloc(void* context, void* memory, size_t size)
{
  return realloc(memory, size);
}

static void std_free(void* context, void* memory)
{
  free(memory);
}

polymec_allocator_t* std_allocator_new()
{
  polymec_allocator_vtable vtable = {.malloc = std_malloc,
                                     .aligned_alloc = std_aligned_alloc,
                                     .realloc = std_realloc,
                                     .free = std_free};
  return polymec_allocator_new("Standard", NULL, vtable);
}

static void* my_arena_malloc(void* context, size_t size)
{
  ARENA* arena = context;
  return arena_malloc(arena, size, 0);
}

static void* my_arena_aligned_alloc(void* context, size_t alignment, size_t size)
{
  ARENA* arena = context;
  return arena_malloc(arena, size, alignment);
}

static void* my_arena_realloc(void* context, void* memory, size_t size)
{
  ARENA* arena = context;
  return arena_realloc(arena, memory, size, 0);
}

static void* my_arena_aligned_realloc(void* context, void* memory, size_t alignment, size_t size)
{
  ARENA* arena = context;
  return arena_realloc(arena, memory, size, alignment);
}

static void my_arena_free(void* context, void* memory)
{
  ARENA* arena = context;
  arena_free(arena, memory);
}

static void my_arena_dtor(void* context)
{
  ARENA* arena = context;
  arena_close(arena);
}

polymec_allocator_t* arena_allocator_new()
{
  polymec_allocator_vtable vtable = {.malloc = my_arena_malloc,
                                     .aligned_alloc = my_arena_aligned_alloc,
                                     .realloc = my_arena_realloc,
                                     .free = my_arena_free,
                                     .dtor = my_arena_dtor};
  ARENA* arena = arena_open(&arena_defaults, NULL);
  return polymec_allocator_new("Arena", arena, vtable);
}

static void* my_pool_malloc(void* context, size_t size)
{
  POOL* pool = context;
  return pool_get(pool, size, 0);
}

static void* my_pool_aligned_alloc(void* context, size_t alignment, size_t size)
{
  POOL* pool = context;
  return pool_get(pool, size, alignment);
}

static void* my_pool_realloc(void* context, void* memory, size_t size)
{
  POOL* pool = context;
  return pool_realloc(pool, memory, size, 0);
}

static void* my_pool_aligned_realloc(void* context, void* memory, size_t alignment, size_t size)
{
  POOL* pool = context;
  return pool_realloc(pool, memory, size, alignment);
}

static void my_pool_free(void* context, void* memory)
{
  POOL* pool = context;
  pool_put(pool, memory);
}

static void my_pool_dtor(void* context)
{
  POOL* pool = context;
  pool_close(pool);
}

polymec_allocator_t* pool_allocator_new()
{
  polymec_allocator_vtable vtable = {.malloc = my_pool_malloc,
                                     .aligned_alloc = my_pool_aligned_alloc,
                                     .realloc = my_pool_realloc,
                                     .free = my_pool_free,
                                     .dtor = my_pool_dtor};
  POOL* pool = pool_open(&pool_defaults, NULL);
  return polymec_allocator_new("Pool", pool, vtable);
}

