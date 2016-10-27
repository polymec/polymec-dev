// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "gc/gc.h"
#include "core/allocators.h"
#include "core/slist.h"
#include "core/logging.h"
#include "arena/proto.h"
#include "arena/pool.h"

struct polymec_allocator_t 
{
  char* name;
  void* context;
  polymec_allocator_vtable vtable;
};

static polymec_allocator_t* polymec_allocator_new(const char* name,
                                                  void* context,
                                                  polymec_allocator_vtable vtable)
{
  ASSERT(vtable.malloc != NULL);
  ASSERT(vtable.realloc != NULL);
  ASSERT(vtable.free != NULL);

  polymec_allocator_t* alloc = malloc(sizeof(polymec_allocator_t)); // Oh, the irony...
  alloc->name = malloc(sizeof(char) * (strlen(name) + 1)); // Need to avoid using other allocators.
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

  // Before we free the allocator stack, we have to empty it. This prevents
  // ptr_slist_free() from using any allocator to free itself.
  ptr_slist_clear(alloc_stack);

  // Now do the deed.
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

//------------------------------------------------------------------------
static void* std_malloc(void* context, size_t size)
{
  // We want to trace pointers within this memory, so we need to use the 
  // garbage collector's allocator.
  return GC_malloc_uncollectable(size);
}

static void* std_aligned_alloc(void* context, size_t alignment, size_t size)
{
  // FIXME: Aligned allocations are not traceable at the moment.
#if defined(APPLE) || defined(__INTEL_COMPILER)
  polymec_error("Aligned allocations are not available on MacOS/Intel at this time.");
  return NULL;
#else
  return aligned_alloc(alignment, size);
#endif
}

static void* std_realloc(void* context, void* memory, size_t size)
{
  return GC_REALLOC(memory, size);
}

typedef void (*gc_finalizer)(void* context, void* dummy);
static void* std_gc_malloc(void* context, size_t size, void (*dtor)(void* context))
{
  void* memory = GC_MALLOC(size);
  if (dtor != NULL)
  {
    gc_finalizer finalizer = (gc_finalizer)dtor;
    GC_register_finalizer(memory, finalizer, memory, NULL, NULL);
  }
  return memory;
}

static void std_free(void* context, void* memory)
{
  GC_FREE(memory);
}
//------------------------------------------------------------------------

void* polymec_malloc(size_t size)
{
  if ((alloc_stack == NULL) || (alloc_stack->size == 0))
    return std_malloc(NULL, size);
  else
  {
    polymec_allocator_t* alloc = alloc_stack->front->value;
    return alloc->vtable.malloc(alloc->context, size);
  }
}

void* polymec_aligned_alloc(size_t alignment, size_t size)
{
  if ((alloc_stack == NULL) || (alloc_stack->size == 0))
    return std_aligned_alloc(NULL, alignment, size);
  else
  {
    polymec_allocator_t* alloc = alloc_stack->front->value;
    if (alloc->vtable.aligned_alloc != NULL)
      return alloc->vtable.aligned_alloc(alloc->context, alignment, size);
    else
    {
      polymec_error("%s allocator does not support aligned allocation.");
      return NULL;
    }
  }
}

void* polymec_realloc(void* memory, size_t size)
{
  if ((alloc_stack == NULL) || (alloc_stack->size == 0))
    return std_realloc(NULL, memory, size);
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
    return NULL;
  }
  else
  {
    polymec_allocator_t* alloc = alloc_stack->front->value;
    if (alloc->vtable.aligned_realloc != NULL)
      return alloc->vtable.aligned_realloc(alloc->context, memory, alignment, size);
    else
    {
      polymec_error("%s allocator does not support aligned allocation.");
      return NULL;
    }
  }
}

void* polymec_gc_malloc(size_t size, void (*dtor)(void* memory))
{
  if ((alloc_stack == NULL) || (alloc_stack->size == 0))
    return std_gc_malloc(NULL, size, dtor);
  polymec_allocator_t* alloc = alloc_stack->front->value;
  if (alloc->vtable.gc_malloc != NULL) 
    return alloc->vtable.gc_malloc(alloc->context, size, dtor);
  else
  {
    static bool first_time = true;
    if (first_time)
    {
      log_urgent("%s allocator lacks garbage collection!", alloc->name);
      log_urgent("Using regular malloc.");
      first_time = false;
    }
    return alloc->vtable.malloc(alloc->context, size);
  }
}

void polymec_free(void* memory)
{
  if ((alloc_stack == NULL) || (alloc_stack->size == 0))
    std_free(NULL, memory);
  else 
  {
    polymec_allocator_t* alloc = alloc_stack->front->value;
    alloc->vtable.free(alloc->context, memory);
  }
}

polymec_allocator_t* std_allocator_new()
{
  polymec_allocator_vtable vtable = {.malloc = std_malloc,
                                     .aligned_alloc = std_aligned_alloc,
                                     .realloc = std_realloc,
                                     .gc_malloc = std_gc_malloc,
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

