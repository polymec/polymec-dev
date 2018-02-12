// Copyright (c) 2012-2018, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/allocators.h"
#include "core/slist.h"
#include "core/logging.h"
#include "arena/proto.h"
#include "arena/pool.h"

#include "lua.h"
#include "lualib.h"
#include "lauxlib.h"

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
  alloc_stack = NULL;
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
  return malloc(size);
}

static void* std_realloc(void* context, void* memory, size_t size)
{
  return realloc(memory, size);
}

static void std_free(void* context, void* memory)
{
  free(memory);
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

void* polymec_calloc(size_t size)
{
  void* memory = polymec_malloc(size);
  memset(memory, 0, size);
  return memory;
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

//------------------------------------------------------------------------ 
//                       Garbage collection
//------------------------------------------------------------------------ 
// We use Lua's garbage collector, since this is the only way to provide 
// a single context for lifetimes of collected objects shared between 
// Lua and C.
//------------------------------------------------------------------------ 

// This unpublished function lets us access the main Lua state initiated 
// by polymec_init.
extern lua_State* polymec_lua_State(void);

// We use this to store a function pointer portably.
typedef struct
{
  void (*finalize)(void* context);
} gc_finalizer_t;

// This is the destructor/finalizer that gets called by Lua when an 
// object is collected.
static int lua_finalize_object(lua_State* L)
{
  // Fetch our userdata.
  void* storage = lua_touserdata(L, 1);
  ASSERT(storage != NULL);

  // Fetch our destructor.
  gc_finalizer_t* storage_f = (gc_finalizer_t*)storage;
  void (*finalize)(void*) = storage_f[2].finalize;

  // Fetch our object.
  void* obj = &(storage_f[3]);

  // Finalize the object if needed.
  if (finalize != NULL)
    finalize(obj);
  
  return 0;
}

// Type identifier for garbage-collected C objects.
static const char* GC_C_OBJ = "gc_C_obj";

void* polymec_gc_malloc(size_t size, void (*finalize)(void* memory))
{
  lua_State* L = polymec_lua_State();
  ASSERT(L != NULL);

  // The first time we do this, we create a metatable for all 
  // C garbage-collected objects.
  static bool first_time = true;
  if (first_time)
  {
    luaL_newmetatable(L, GC_C_OBJ);
    lua_pushcfunction(L, lua_finalize_object);
    lua_setfield(L, -2, "__gc");
    first_time = false;
  }

  // We store a reference to each garbage-collected object in Lua's registry.
  // We refer to this reference using an integer key. We also count the number 
  // of C references generated by polymec_retain. So we allocate storage 
  // for a userdata of the desired size preceded by 2 ints and a finalizer.
  void* storage = lua_newuserdata(L, 3 * sizeof(gc_finalizer_t) + size);

  // Set the metatable for this object.
  luaL_getmetatable(L, GC_C_OBJ);
  lua_setmetatable(L, -2);

  // Get a key for this object and stash it at the front of the userdata.
  int ref = luaL_ref(L, LUA_REGISTRYINDEX);
  int* storage_ref = (int*)storage;
  storage_ref[0] = ref;

  // The reference count for a newly created object is 1.
  storage_ref[1] = 1;

  // Stash the destructor right behind it.
  gc_finalizer_t* storage_f = (gc_finalizer_t*)storage;
  storage_f[2].finalize = finalize;

  // The actual object goes behind the ref and the destructor.
  storage = &(storage_f[3]);
  return storage;
}

// Use this to get to the head of a collected object's block of memory.
static inline void* collected_obj_header(void* memory)
{
  gc_finalizer_t* storage_f = memory;
  return (void*)&(storage_f[-3]);
}

void polymec_retain(void* memory)
{
  // Fetch the reference for this object and make sure it exists in the 
  // Lua registry.
  int* storage_ref = (int*)collected_obj_header(memory);
  int ref = storage_ref[0];
  lua_State* L = polymec_lua_State();
  lua_rawgeti(L, LUA_REGISTRYINDEX, ref);
  if (luaL_checkudata(L, -1, GC_C_OBJ) == NULL)
    polymec_error("polymec_retain called on non-garbage-collected object!");
  lua_pop(L, 1);

  // Increment the C reference count for this object.
  ASSERT(storage_ref[1] >= 0);
  ++(storage_ref[1]);
}

void polymec_release(void* memory)
{
  ASSERT(memory != NULL);

  // Fetch the reference for this object and make sure it exists in the 
  // Lua registry.
  int* storage_ref = (int*)collected_obj_header(memory);
  int ref = storage_ref[0];
  lua_State* L = polymec_lua_State();
  lua_rawgeti(L, LUA_REGISTRYINDEX, ref);
  if (luaL_checkudata(L, -1, GC_C_OBJ) == NULL)
    polymec_error("polymec_release called on non-garbage-collected object!");
  lua_pop(L, 1);

  // Decrement our C reference count.
  --storage_ref[1];
  int ref_count = storage_ref[1];

  // If we've hit zero, delete the reference from the registry.
  if (ref_count <= 0)
  {
    ASSERT(L != NULL);
    luaL_unref(L, LUA_REGISTRYINDEX, ref);
  }
}

