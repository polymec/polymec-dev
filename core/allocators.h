// Copyright (c) 2012-2018, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_ALLOCATORS_H
#define POLYMEC_ALLOCATORS_H

#include <stdlib.h>

/// \addtogroup core core
///@{

/// \class polymec_allocator
/// This type represents a memory allocator which services memory allocations 
/// and deallocations. There are two kinds of allocators: an "arena," which 
/// represents a contiguous block of memory in the heap which grows as needed, 
/// and a "pool," which allocates pools of memory of fixed sizes.
typedef struct polymec_allocator_t polymec_allocator_t;

/// \struct polymec_allocator_vtable
/// Any implementation of an allocator must implement this virtual table.
typedef struct 
{
  /// Allocates and returns size bytes of memory.
  void* (*malloc)(void* context, size_t size);

  /// Reallocates size bytes of the given chunk memory.
  void* (*realloc)(void* context, void* memory, size_t size);

  /// Frees the given chunk of memory.
  void (*free)(void* context, void* memory);
  
  /// Context destructor.
  void (*dtor)(void* context);
} polymec_allocator_vtable;

/// Returns an allocator that uses the standard C malloc/free functions.
/// \relates polymec_allocator
polymec_allocator_t* std_allocator_new(void);

/// Returns an arena allocator.
/// \relates polymec_allocator
polymec_allocator_t* arena_allocator_new(void);

/// Returns a pool allocator.
/// \relates polymec_allocator
polymec_allocator_t* pool_allocator_new(void);

/// Destroys the given allocator.
/// \memberof polymec_allocator
void polymec_allocator_free(polymec_allocator_t* alloc);

/// This function allocates memory in the same fashion as malloc(), using the 
/// allocator on the top of polymec's allocator stack. If the stack is empty, 
/// calls to polymec_malloc() simply use malloc().
void* polymec_malloc(size_t size);

/// This function allocates and zeros memory in the same way as calloc(), using the
/// allocator on top of polymec's allocator stack. 
void* polymec_calloc(size_t size);

/// This reallocates existing memory using the allocator on top of the allocator 
/// stack, or calls realloc() if the stack is empty.
void* polymec_realloc(void* memory, size_t size);

/// Frees memory, returning it to the allocator on top of the allocator stack 
/// (or calling free() if the stack is empty).
void polymec_free(void* memory);

/// Returns memory that will be garbage-collected. This memory should not be 
/// freed--use polymec_release instead of polymec_free 
/// to tell the collector that you're finished with a garbage-collected 
/// resource. When an object is collected, the supplied finalizer is 
/// invoked.
/// Garbage collection is appropriate for objects that are shared by many 
/// systems, and that don't consume large amounts of resources. There is no 
/// garbage-collected realloc (polymec_gc_realloc), since we don't encourage 
/// garbage collection for data-intensive objects.
void* polymec_gc_malloc(size_t size, void (*finalize)(void* memory));

/// Call polymec_retain when you need to retain a reference to a garbage-
/// collected object in a C data structure or context. Calling polymec_retain
/// on a resource not allocated by polymec_gc_malloc throws a fatal 
/// runtime error.
void polymec_retain(void* memory);

/// Call polymec_release when you are finished with a garbage-collected 
/// resource so that it can be collected. Calling polymec_release on a 
/// resource not allocated by polymec_gc_malloc throws a fatal runtime error.
void polymec_release(void* memory);

/// Pushes a new memory allocator to the allocator stack, using this allocator 
/// to allocate all memory with polymec_malloc() until it is popped.
void push_allocator(polymec_allocator_t* allocator);

/// Pops an allocator off of the stack, returning it and using the next one to 
/// allocate all memory with polymec_malloc(). 
polymec_allocator_t* pop_allocator(void);

///@}

#endif
