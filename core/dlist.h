// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_DLIST_H
#define POLYMEC_DLIST_H

#include "core/polymec.h"
#include "core/comparators.h"

/// \addtogroup core core
///@{

/// \def DEFINE_DLIST(list_name, element)
/// Defines a doubly-linked list for a given element type. The following
/// interface is defined for a list with name `x_dlist`.
/// * `x_dlist_t* x_dlist_new()` - Creates a new, empty dlist.
/// * `void x_dlist_free(dlist_t* list)` - Destroys the list.
/// * `x_dlist_node_t* x_dlist_find(x_dlist_t* list, x value, x_dlist_cmp comparator)` - Returns the node at which a value appears in the list.
/// * `void x_dlist_insert(x_dlist_t* list, x value, x_dlist_node_t* node)` - Inserts an x into the list in front of the given node.
/// * `void x_dlist_insert_with_dtor(x_dlist_t* list, x value, x_dlist_node_t* node, destructor dtor)` - Inserts an x into the list, using dtor to destroy it when finished.
/// * `void x_dlist_append(x_dlist_t* list, x value)` - Appends an x to the end of the list.
/// * `void x_dlist_append_with_dtor(x_dlist_t* list, x value, destructor dtor)` - Appends an x to the end of the list, using dtor to destroy when finished.
/// * `void x_dlist_push(x_dlist_t* list, x value)` - Inserts an x at the front of the list.
/// * `void x_dlist_push_with_dtor(x_dlist_t* list, x value, x_dlist_dtor dtor)` - Inserts an x at the front of the list with a destructor.
/// * `x x_dlist_pop(x_dlist_t* list, x_dlist_dtor* dtor)` - Removes an x from the front of the list, returning it and its destructor (if dtor != NULL).
/// * `void x_dlist_remove(x_dlist_t* list, x_dlist_node_t* node)` - Removes a node from the list.
/// * `bool x_dlist_next(x_dlist_t* list, x_dlist_node_t** pos, x* value)` - Allows the traversal of the linked list.
/// * `bool x_dlist_prev(x_dlist_t* list, x_dlist_node_t** pos, x* value)` - Allows a reversed traversal of the linked list.
/// * `bool x_dlist_empty(x_dlist_t* list)` - Returns true if empty, false otherwise.
/// * `void x_dlist_clear(x_dlist_t* list)` - Clears the given list, making it empty.
///
/// Data members for a list `list`:
/// * `list->front` - The node (with members `value`, `dtor`, `prev`, and `next` at the front of the list.
/// * `list->back` - The node at the back of the list.
/// * `list->size` - The number of elements in the list.
///
/// \param list_name The name of the doubly-linked list.
/// \param element The data type stored in the list.
#define DEFINE_DLIST(list_name, element) \
typedef struct list_name##_node_t list_name##_node_t; \
typedef void (*list_name##_dtor)(element); \
struct list_name##_node_t \
{ \
  element value; \
  list_name##_dtor dtor; \
  list_name##_node_t* prev; \
  list_name##_node_t* next; \
}; \
\
typedef struct list_name##_t list_name##_t; \
struct list_name##_t \
{ \
  list_name##_node_t* front; \
  list_name##_node_t* back; \
  size_t size; \
}; \
\
typedef int (*list_name##_comparator)(element, element); \
\
static inline list_name##_t* list_name##_new() \
{ \
  list_name##_t* list = (list_name##_t*)polymec_malloc(sizeof(list_name##_t)); \
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
  list_name##_node_t* n = (list_name##_node_t*)polymec_malloc(sizeof(list_name##_node_t)); \
  n->value = value; \
  n->next = NULL; \
  n->prev = NULL; \
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
    n->prev = NULL; \
    list->front = n; \
    list->back = n->next; \
    list->size = 2; \
  } \
  else \
  { \
    list_name##_node_t* p = list->front; \
    list_name##_node_t* p_prev = NULL; \
    while (p->next != node) \
    { \
      ASSERT(p != NULL); \
      p_prev = p; \
      p = p->next; \
    } \
    n->next = p->next; \
    p->next = n; \
    p->prev = p_prev; \
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
  list_name##_node_t* n = (list_name##_node_t*)polymec_malloc(sizeof(list_name##_node_t)); \
  n->value = value; \
  n->dtor = dtor; \
  n->next = NULL; \
  if (list->back == NULL) \
  { \
    ASSERT(list->front == NULL); \
    n->prev = NULL; \
    list->front = n; \
    list->back = n; \
  } \
  else \
  { \
    list->back->next = n; \
    n->prev = list->back; \
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
  if (list->size == 0) \
    return; \
  if (list->front == node) \
  { \
    list_name##_dtor dtor; \
    element e = list_name##_pop(list, &dtor); \
    if (dtor != NULL) \
      dtor(e); \
    return; \
  } \
  node->prev->next = node->next; \
  if (node->next != NULL) \
    node->next->prev = node->prev; \
  if (node->dtor != NULL) \
    (node->dtor)(node->value); \
  polymec_free(node); \
  list->size -= 1; \
} \
\
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
static inline bool list_name##_prev(list_name##_t* list, list_name##_node_t** node, element* value) \
{ \
  if (*node == NULL) \
    *node = list->back; \
  else \
    *node = (*node)->prev; \
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

///@}

// Define some basic dlist types.
DEFINE_DLIST(int_dlist, int)
DEFINE_DLIST(long_dlist, long)
DEFINE_DLIST(real_dlist, real_t)
DEFINE_DLIST(string_dlist, char*)
DEFINE_DLIST(ptr_dlist, void*)

#endif
