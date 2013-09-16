// Copyright 2012-2013 Jeffrey Johnson.
// 
// This file is part of Polymec, and is licensed under the Apache License, 
// Version 2.0 (the "License"); you may not use this file except in 
// compliance with the License. You may may find the text of the license in 
// the LICENSE file at the top-level source directory, or obtain a copy of 
// it at
// 
// http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef POLYMEC_SLIST_H
#define POLYMEC_SLIST_H

#include "core/polymec.h"
#include "arena/proto.h"
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
// void x_slist_insert(x_slist_t* list, x value, x_slist_node_t* node) - Inserts an x into the list.
// void x_slist_insert_with_dtor(x_slist_t* list, x value, x_slist_node_t* node, destructor dtor) - Inserts an x into the list, using dtor to destroy it when finished.
// void x_slist_append(x_slist_t* list, x value) - Appends an x to the end of the list.
// void x_slist_append_with_dtor(x_slist_t* list, x value, destructor dtor) - Appends an x to the end of the list, using dtor to destroy when finished.
// x x_slist_pop(x_slist_t* list, x_slist_dtor* dtor) - Removes an x from the front of the list, returning it and its destructor (if dtor != NULL).
// void x_slist_remove(x_slist_t* list, x_slist_node_t* node) - Removes a node from the list.
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
  ARENA* arena; \
}; \
\
typedef int (*list_name##_comparator)(element, element); \
\
static inline list_name##_t* list_name##_new() \
{ \
  list_name##_t* list = (list_name##_t*)malloc(sizeof(list_name##_t)); \
  list->front = list->back = NULL; \
  list->size = 0; \
  list->arena = NULL; \
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
    free(n); \
  } \
  free(list); \
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
  list_name##_node_t* n = (list_name##_node_t*)malloc(sizeof(list_name##_node_t)); \
  n->value = value; \
  n->next = NULL; \
  if (list->front == NULL) \
  { \
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
  list_name##_node_t* n = (list_name##_node_t*)malloc(sizeof(list_name##_node_t)); \
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
  free(node); \
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
  } \
  list_name##_node_t* p = list->front; \
  while ((p->next != node) && (p != NULL)) \
    p = p->next; \
  if (p != NULL) \
  { \
    p->next = node->next; \
    if (node->dtor != NULL) \
      (node->dtor)(node->value); \
    free(node); \
    list->size -= 1; \
  } \
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
DEFINE_SLIST(double_slist, double)
DEFINE_SLIST(string_slist, char*)
DEFINE_SLIST(ptr_slist, void*)

#endif
