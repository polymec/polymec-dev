#ifndef POLYMEC_SLIST_H
#define POLYMEC_SLIST_H

#include <stdlib.h>
#include <string.h>
#include "core/polymec.h"
#include "arena/proto.h"
#include "core/comparators.h"

// An slist is a singly-linked list that stores homogeneous types.
// One defines an slist using
// DEFINE_SLIST(list_name, element)
//
// Interface for a type x_slist_t (with datum x) defined with 
// DEFINE_SLIST(x_slist, x, dtor):
//
// x_slist_t* x_slist_new() - Creates a new, empty slist.
// void x_slist_free(slist_t* list) - Destroys the list.
// x_slist_node_t* x_slist_find(x_slist_t* list, x value, x_slist_cmp comparator) - Returns the node at which a value appears in the list.
// void x_slist_insert(x_slist_t* list, x value, x_slist_node_t* node) - Inserts an x into the list.
// void x_slist_append(x_slist_t* list, x value) - Appends an x to the end of the list.
// x x_slist_pop(x_slist_t* list) - Removes an x from the front of the list, returning it.
// void x_slist_remove(x_slist_t* list, x_slist_node_t* node) - Removes a node from the list.

#define DEFINE_SLIST(list_name, element) \
typedef struct list_name##_node_t list_name##_node_t; \
struct list_name##_node_t \
{ \
  element value; \
  list_name##_node_t* next; \
}; \
\
typedef struct list_name##_t list_name##_t; \
struct list_name##_t \
{ \
  list_name##_node_t* front; \
  list_name##_node_t* back; \
  int size; \
  bool owns_data; \
  ARENA* arena; \
}; \
\
typedef int (*list_name##_comparator)(element, element); \
\
static inline list_name##_t* list_name##_new() \
{ \
  list_name##_t* list = malloc(sizeof(list_name##_t)); \
  list->front = list->back = NULL; \
  list->size = 0; \
  list->owns_data = true; \
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
static inline void list_name##_insert(list_name##_t* list, element value, list_name##_node_t* node) \
{ \
  ASSERT(node != NULL); \
  list_name##_node_t* n = malloc(sizeof(list_name##_node_t)); \
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
static inline void list_name##_append(list_name##_t* list, element value) \
{ \
  list_name##_node_t* n = malloc(sizeof(list_name##_node_t)); \
  n->value = value; \
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
\
static inline element list_name##_pop(list_name##_t* list) \
{ \
  ASSERT(list->size > 0); \
  list_name##_node_t* node = list->front; \
  list->front = list->front->next; \
  if (list->size == 1) \
    list->back = NULL; \
  list->size -= 1; \
  element val = node->value; \
  free(node); \
  return val; \
} \
\
static inline void list_name##_remove(list_name##_t* list, list_name##_node_t* node) \
{ \
  if (list->front == node) \
    list_name##_pop(list); \
  list_name##_node_t* p = list->front; \
  while ((p->next != node) && (p != NULL)) \
    p = p->next; \
  if (p != NULL) \
  { \
    p->next = node->next; \
    free(node); \
    list->size -= 1; \
  } \
} \

// Define some basic slist types.
DEFINE_SLIST(int_slist, int)
DEFINE_SLIST(double_slist, double)
DEFINE_SLIST(string_slist, char*)
DEFINE_SLIST(ptr_slist, void*)

#endif
