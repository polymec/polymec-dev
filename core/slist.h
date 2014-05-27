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

#ifndef POLYMEC_SLIST_H
#define POLYMEC_SLIST_H

#include "core/polymec.h"
#include "core/comparators.h"

// An slist is a singly-linked list that stores homogeneous types.
// One defines an slist using
// DEFINE_SLIST(list_name, element)
//
// Interface for a type x_slist_t (with datum x) defined with 
// DEFINE_SLIST(x_slist, x):
//
// x_slist_t* x_slist_new() - Creates a new, empty slist.
// void x_slist_free(slist_t* list) - Destroys the list.
// x_slist_node_t* x_slist_find(x_slist_t* list, x value, x_slist_cmp comparator) - Returns the node at which a value appears in the list.
// void x_slist_insert(x_slist_t* list, x value, x_slist_node_t* node) - Inserts an x into the list in front of the given node.
// void x_slist_insert_with_dtor(x_slist_t* list, x value, x_slist_node_t* node, destructor dtor) - Inserts an x into the list, using dtor to destroy it when finished.
// void x_slist_append(x_slist_t* list, x value) - Appends an x to the end of the list.
// void x_slist_append_with_dtor(x_slist_t* list, x value, destructor dtor) - Appends an x to the end of the list, using dtor to destroy when finished.
// void x_slist_push(x_slist_t* list, x value) - Inserts an x at the front of the list.
// void x_slist_push_with_dtor(x_slist_t* list, x value, x_slist_dtor dtor) - Inserts an x at the front of the list with a destructor.
// x x_slist_pop(x_slist_t* list, x_slist_dtor* dtor) - Removes an x from the front of the list, returning it and its destructor (if dtor != NULL).
// void x_slist_remove(x_slist_t* list, x_slist_node_t* node) - Removes a node from the list.
// bool x_slist_next(x_slist_t* list, x_slist_node_t** pos, x* value) - Allows the traversal of the linked list.
// bool x_slist_empty(x_slist_t* list) - Returns true if empty, false otherwise.
// void x_slist_clear(x_slist_t* list) - Clears the given list, making it empty.

#define DEFINE_SLIST(list_name, element) \
typedef struct list_name##_node_t list_name##_node_t; \
typedef void (*list_name##_dtor)(element); \
struct list_name##_node_t \
{ \
  element value; \
  list_name##_dtor dtor; \
  list_name##_node_t* next; \
}; \
\
typedef struct list_name##_t list_name##_t; \
struct list_name##_t \
{ \
  list_name##_node_t* front; \
  list_name##_node_t* back; \
  int size; \
}; \
\
typedef int (*list_name##_comparator)(element, element); \
\
static inline list_name##_t* list_name##_new() \
{ \
  list_name##_t* list = polymec_malloc(sizeof(list_name##_t)); \
  list->front = list->back = NULL; \
  list->size = 0; \
  return list; \
} \
\
static inline void list_name##_free(list_name##_t* list) \
{ \
  list_name##_node_t* n; \
  while (list->front != NULL) \
  { \
    n = list->front; \
    list->front = n->next; \
    if (n->dtor != NULL) \
      (n->dtor)(n->value); \
    polymec_free(n); \
  } \
  polymec_free(list); \
} \
\
static inline list_name##_node_t* list_name##_find(list_name##_t* list, element value, list_name##_comparator comparator) \
{ \
  list_name##_node_t* n = list->front; \
  while ((n != NULL) && (comparator(n->value, value) != 0)) \
    n = n->next; \
  return n; \
} \
static inline void list_name##_insert_with_dtor(list_name##_t* list, element value, list_name##_dtor dtor, list_name##_node_t* node) \
{ \
  ASSERT(node != NULL); \
  list_name##_node_t* n = polymec_malloc(sizeof(list_name##_node_t)); \
  n->value = value; \
  n->next = NULL; \
  if (list->front == NULL) \
  { \
    ASSERT(node == NULL); \
    list->front = n; \
    list->back = n; \
    list->size = 1; \
  } \
  else if (list->front == list->back) \
  { \
    ASSERT(node == list->front); \
    n->next = list->front; \
    list->front = n; \
    list->back = n->next; \
    list->size = 2; \
  } \
  else \
  { \
    list_name##_node_t* p = list->front; \
    while (p->next != node) \
    { \
      ASSERT(p != NULL); \
      p = p->next; \
    } \
    n->next = p->next; \
    p->next = n; \
    list->size += 1; \
  } \
} \
\
static inline void list_name##_insert(list_name##_t* list, element value, list_name##_node_t* node) \
{ \
  list_name##_insert_with_dtor(list, value, NULL, node); \
} \
static inline void list_name##_append_with_dtor(list_name##_t* list, element value, list_name##_dtor dtor) \
{ \
  list_name##_node_t* n = polymec_malloc(sizeof(list_name##_node_t)); \
  n->value = value; \
  n->dtor = dtor; \
  n->next = NULL; \
  if (list->back == NULL) \
  { \
    ASSERT(list->front == NULL); \
    list->front = n; \
    list->back = n; \
  } \
  else \
  { \
    list->back->next = n; \
    list->back = n; \
  } \
  list->size += 1; \
} \
static inline void list_name##_append(list_name##_t* list, element value) \
{ \
  list_name##_append_with_dtor(list, value, NULL); \
} \
\
static inline void list_name##_push_with_dtor(list_name##_t* list, element value, list_name##_dtor dtor) \
{ \
  list_name##_insert_with_dtor(list, value, dtor, list->front); \
} \
\
static inline void list_name##_push(list_name##_t* list, element value) \
{ \
  list_name##_push_with_dtor(list, value, NULL); \
} \
\
static inline element list_name##_pop(list_name##_t* list, list_name##_dtor* dtor) \
{ \
  ASSERT(list->size > 0); \
  list_name##_node_t* node = list->front; \
  list->front = list->front->next; \
  if (list->size == 1) \
    list->back = NULL; \
  list->size -= 1; \
  element val = node->value; \
  if (dtor != NULL) \
    *dtor = node->dtor; \
  polymec_free(node); \
  return val; \
} \
\
static inline void list_name##_remove(list_name##_t* list, list_name##_node_t* node) \
{ \
  if (list->front == node) \
  { \
    list_name##_dtor dtor; \
    element e = list_name##_pop(list, &dtor); \
    if (dtor != NULL) \
      dtor(e); \
    return; \
  } \
  list_name##_node_t* p = list->front; \
  while ((p->next != node) && (p != NULL)) \
    p = p->next; \
  if (p != NULL) \
  { \
    p->next = node->next; \
    if (node->dtor != NULL) \
      (node->dtor)(node->value); \
    polymec_free(node); \
    list->size -= 1; \
  } \
} \
static inline bool list_name##_next(list_name##_t* list, list_name##_node_t** node, element* value) \
{ \
  if (*node == NULL) \
    *node = list->front; \
  else \
    *node = (*node)->next; \
  if (*node != NULL) \
    *value = (*node)->value; \
  return (*node != NULL); \
} \
\
static inline bool list_name##_empty(list_name##_t* list) \
{ \
  return (list->size == 0); \
} \
\
static inline void list_name##_clear(list_name##_t* list) \
{ \
  while (list->front != NULL) \
    list_name##_pop(list, NULL); \
} \

// Define some basic slist types.
DEFINE_SLIST(int_slist, int)
DEFINE_SLIST(long_slist, long)
DEFINE_SLIST(real_slist, real_t)
DEFINE_SLIST(string_slist, char*)
DEFINE_SLIST(ptr_slist, void*)

#endif
