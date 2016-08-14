// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_ALLOCATORS_H
#define POLYMEC_ALLOCATORS_H

#include <stdlib.h>

// This type represents a memory allocator which services memory allocations 
// and deallocations. There are two kinds of allocators: an "arena," which 
// represents a contiguous block of memory in the heap which grows as needed, 
// and a "pool," which allocates pools of memory of fixed sizes.
typedef struct polymec_allocator_t polymec_allocator_t;

// Any implementation of an allocator must implement this virtual table.
typedef struct 
{
  // Allocates and returns size bytes of memory.
  void* (*malloc)(void* context, size_t size);

  // Allocates and returns size bytes of memory with the given alignment.
  void* (*aligned_alloc)(void* context, size_t alignment, size_t size);

  // Reallocates size bytes of the given chunk memory.
  void* (*realloc)(void* context, void* memory, size_t size);

  // Reallocates size bytes of the given chunk memory according to the given alignment.
  // This might not be available on every allocator.
  void* (*aligned_realloc)(void* context, void* memory, size_t alignment, size_t size);

  // Frees the given chunk of memory.
  void (*free)(void* context, void* memory);
  
  // Context destructor.
  void (*dtor)(void* context);
} polymec_allocator_vtable;

// Returns an allocator that uses the standard C malloc/free functions.
polymec_allocator_t* std_allocator_new(void);

// Returns an arena allocator.
polymec_allocator_t* arena_allocator_new(void);

// Returns a pool allocator.
polymec_allocator_t* pool_allocator_new(void);

// Destroys the given allocator.
void polymec_allocator_free(polymec_allocator_t* alloc);

// This function allocates memory in the same fashion as malloc(), using the 
// allocator on the top of polymec's allocator stack. If the stack is empty, 
// calls to polymec_malloc() simply use malloc().
void* polymec_malloc(size_t size);

// This version of polymec_malloc() returns memory with the given alignment 
// using the allocator on the top of the allocator stack, or calls 
// aligned_alloc() if the stack is empty.
void* polymec_aligned_alloc(size_t alignment, size_t size);

// This reallocates existing memory using the allocator on top of the allocator 
// stack, or calls realloc() if the stack is empty.
void* polymec_realloc(void* memory, size_t size);

// This version of polymec_realloc guarantees the given alignment of the 
// reallocated memory, and using the allocator on top of the allocator stack. 
// If the stack is empty, a fatal error is issued, since there is no 
// portable way to perform an aligned realloc with the standard C allocator.
void* polymec_aligned_realloc(void* memory, size_t alignment, size_t size);

// Frees memory, returning it to the allocator on top of the allocator stack 
// (or calling free() if the stack is empty).
void polymec_free(void* memory);

// Pushes a new memory allocator to the allocator stack, using this allocator 
// to allocate all memory with polymec_malloc() until it is popped.
void push_allocator(polymec_allocator_t* allocator);

// Pops an allocator off of the stack, returning it and using the next one to 
// allocate all memory with polymec_malloc(). 
polymec_allocator_t* pop_allocator(void);

#endif
