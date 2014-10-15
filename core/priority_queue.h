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

#ifndef POLYMEC_PRIORITY_QUEUE_H
#define POLYMEC_PRIORITY_QUEUE_H

#include "core/heap.h"

// A priority queue is a container that maintains its elements in order 
// according to their given (integer-valued) priority. One defines a priority 
// queue using 
// DEFINE_PRIORITY_QUEUE(queue_name, element)
//
// Interface for a type x_priority_queue_t (with element type x) defined with 
// DEFINE_PRIORITY_QUEUE(x_priority_queue, x):
// 
// x_priority_queue_t* x_priority_queue_new() - Creates a new empty priority queue.
// void x_priority_queue_free(x_priority_queue_t* queue) - Destroys the priority queue.
// void x_priority_queue_push_with_dtor(x_priority_queue_t* queue, x datum, int priority, dtor) - Inserts a datum into the queue with a destructor.
// void x_priority_queue_push(x_priority_queue_t* queue, x datum, int priority) - Inserts a datum into the queue.
// int x_priority_queue_pop(x_priority_queue_t* queue) - Removes the element with the largest priority (and its priority) from the queue.
// int x_priority_queue_highest(x_priority_queue_t* queue) - Returns thelargest priority in the queue.
// x x_priority_queue_front(x_priority_queue_t* queue) - Returns the element with the largest priority in the queue.
// bool x_priority_queue_empty(x_priority_queue_t* queue) - Returns true if the queue is empty, false otherwise.
// void x_priority_queue_clear(x_priority_queue_t* queue) - Clears the contents of the queue.
#define DEFINE_PRIORITY_QUEUE(queue_name, element) \
typedef element queue_name##_element_t; \
typedef struct queue_name##_t queue_name##_t; \
typedef int (*queue_name##_comparator_t)(element left, element right); \
typedef struct \
{ \
  int priority; \
  element value; \
  void (*dtor)(element e); \
} queue_name##_datum_t; \
\
static int queue_name##_cmp(queue_name##_datum_t l, queue_name##_datum_t r) \
{ \
  return int_cmp(l.priority, r.priority); \
} \
\
DEFINE_HEAP(queue_name##_heap, queue_name##_datum_t, queue_name##_cmp) \
\
struct queue_name##_t \
{ \
  queue_name##_heap_t* heap; \
}; \
\
static inline queue_name##_t* queue_name##_new() \
{ \
  queue_name##_t* queue = polymec_malloc(sizeof(queue_name##_t)); \
  queue->heap = queue_name##_heap_new(); \
  return queue; \
} \
\
static inline void queue_name##_free(queue_name##_t* queue) \
{ \
  queue_name##_heap_free(queue->heap); \
  polymec_free(queue); \
} \
\
static inline void queue_name##_datum_dtor(queue_name##_datum_t datum) \
{ \
  if (datum.dtor != NULL) \
    datum.dtor(datum.value); \
} \
\
static inline void queue_name##_push_with_dtor(queue_name##_t* queue, queue_name##_element_t value, int priority, void (*dtor)(element e)) \
{ \
  queue_name##_datum_t datum = {.priority = priority, .value = value, .dtor = dtor}; \
  queue_name##_heap_push_with_dtor(queue->heap, datum, queue_name##_datum_dtor); \
} \
\
static inline void queue_name##_push(queue_name##_t* queue, queue_name##_element_t value, int priority) \
{ \
  queue_name##_datum_t datum = {.priority = priority, .value = value}; \
  queue_name##_heap_push(queue->heap, datum); \
} \
\
static inline void queue_name##_pop(queue_name##_t* queue) \
{ \
  queue_name##_heap_pop(queue->heap); \
} \
static inline int queue_name##_highest(queue_name##_t* queue) \
{ \
  return queue_name##_heap_front(queue->heap).priority; \
} \
static inline queue_name##_element_t queue_name##_front(queue_name##_t* queue) \
{ \
  return queue_name##_heap_front(queue->heap).value; \
} \
static inline bool queue_name##_empty(queue_name##_t* queue) \
{ \
  return queue_name##_heap_empty(queue->heap); \
} \
static inline void queue_name##_clear(queue_name##_t* queue) \
{ \
  queue_name##_heap_clear(queue->heap); \
} \
\

// Define some priority_queues.
DEFINE_PRIORITY_QUEUE(int_priority_queue, int)
DEFINE_PRIORITY_QUEUE(index_priority_queue, index_t)

#endif
