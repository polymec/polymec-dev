// Copyright (c) 2012-2019, Jeffrey N. Johnson
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

  /// Reallocates size bytes of the given chunk of memory.
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
/// allocator on top of polymec's allocator stack to allocate a contiguous chunk of 
/// memory.
/// \param [in] count The number of objects for which memory is allocated.
/// \param [in] size The size (in bytes) for each object allocated.
void* polymec_calloc(size_t count, size_t size);

/// This reallocates existing memory using the allocator on top of the allocator 
/// stack, or calls realloc() if the stack is empty.
void* polymec_realloc(void* memory, size_t size);

/// Frees memory, returning it to the allocator on top of the allocator stack 
/// (or calling free() if the stack is empty).
void polymec_free(void* memory);

/// Returns storage for a resource to be reference-counted. Reference 
/// counted (or "refcounted") resources keep track of the number of entities 
/// using them. This number is the resource's reference count (or "refcount"). 
/// A newly created refcounted resource has a refcount of 1. Refcounted 
/// resources can be _borrowed_, _retained_, and _released_. 
/// * Borrowing a refcounted object doesn't æffect its reference count. 
///   Thuc, the borrower cannot assume the object is alive under all 
///   circumstances, so borrowing refcounted objects is only safe when 
///   it's known that the lifetime of a refcounted object is stable. It's 
///   good manners to use the \ref borrow_ref annotation when you borrow 
///   a refcounted objbect (but that function doesn't do anything).
/// * Retaining a refcounted object increments that object's refcount. Do this
///   with \ref retain_ref when you want to keep an object alive to guarantee 
///   access to it.
/// * Releasing a refcounted object decrements that object's refcount. Do this
///   with \ref release_ref when you no longer need an object that you've 
///   retained with \ref retain_ref.
///
/// When a refcounted resource's refcount reaches 0, that resource is 
/// destroyed. 
/// \param [in] size The number of bytes to allocate for the refcounted 
///                  resource.
/// \param [in] dtor A destructor to be called on the refcounted resource 
///                  when its refcount reaches zero.
///
/// \refcounted
void* polymec_refcounted_malloc(size_t size, void (*dtor)(void* memory));

/// Articulates that a refcounted resource is being borrowed, returning the 
/// pointer to the borrowed resource.
/// \refcounted
void* borrow_ref(void* refcounted_resource); 

/// Retains a reference to a refcounted resource, incrementing its refcount 
/// and returning a pointer to it.
/// \refcounted
void* retain_ref(void* refcounted_resource);

/// Releases a reference to a refcounted resource, decrementing its refcount.
/// \refcounted
void release_ref(void* refcounted_resource);

/// Returns the reference count for a refcounted resource.
/// \refcounted
int ref_count(void* refcounted_resource);

/// Pushes a new memory allocator to the allocator stack, using this allocator 
/// to allocate all memory with polymec_malloc() until it is popped.
void push_allocator(polymec_allocator_t* allocator);

/// Pops an allocator off of the stack, returning it and using the next one to 
/// allocate all memory with polymec_malloc(). 
polymec_allocator_t* pop_allocator(void);

///@}

#endif
