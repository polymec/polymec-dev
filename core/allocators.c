// Copyright (c) 2012-2017, Jeffrey N. Johnson
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
  ASSERT(name != NULL);
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
  free(alloc->name);
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
  return GC_MALLOC_UNCOLLECTABLE(size);
}

static void* std_realloc(void* context, void* memory, size_t size)
{
  return GC_REALLOC(memory, size);
}

typedef struct
{
  void (*dtor)(void* context);
} gc_adaptor_t;

static void gc_finalizer_adaptor(void* obj, void* client_data)
{
  gc_adaptor_t* adaptor = client_data;
  adaptor->dtor(obj);
  polymec_free(adaptor);
}

static void* std_gc_malloc(void* context, size_t size, void (*dtor)(void* context))
{
  void* memory = GC_MALLOC(size);
  if (dtor != NULL)
  {
    gc_adaptor_t* adaptor = std_malloc(NULL, sizeof(gc_adaptor_t));
    adaptor->dtor = dtor;
    GC_REGISTER_FINALIZER(memory, gc_finalizer_adaptor, adaptor, NULL, NULL);
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
                                     .realloc = my_pool_realloc,
                                     .free = my_pool_free,
                                     .dtor = my_pool_dtor};
  POOL* pool = pool_open(&pool_defaults, NULL);
  return polymec_allocator_new("Pool", pool, vtable);
}

