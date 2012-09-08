#ifndef ARBI_LIST_H
#define ARBI_LIST_H

#include "arbi.h"
#include "tommylist.h"

#ifdef __cplusplus
extern "C" {
#endif

// list_t is a doubly-linked list class based off of the tommyds list class.
// All elements in the list have the same type.
typedef struct list_t list_t;
typedef tommy_node list_node_t;
typedef tommy_compare_func list_comp_func;
typedef void (*list_element_dtor)(void*);

// Creates a new empty list with the given element destructor.
list_t* list_new(list_element_dtor dtor);

// Frees resources used by the given list.
void list_free(list_t* list);

// Clears the contents of the list.
void list_clear(list_t* list);

// Returns the head of the list.
list_node_t* list_head(list_t* list);

// Returns the tail of the list.
list_node_t* list_tail(list_t* list);

// Returns the value associated with the list node.
void* list_node_value(list_node_t* node);

// Returns the size of the list.
int list_size(list_t* list);

// Sorts the elements within the list using the given comparator using 
// a stable, N log N sort.
void list_sort(list_t* list, list_comp_func comparator);

// Inserts the given value at the head of the list.
void list_prepend(list_t* list, void* value);

// Inserts the given value at the end of the list.
void list_append(list_t* list, void* value);

// Inserts the given value before the given node in the list.
void list_insert_before(list_t* list, void* value, list_node_t* node);

// Returns the node in the list containing the given value, or NULL if 
// the list does not contain it.
list_node_t* list_find(list_t* list, void* value, list_comp_func comp);

// Removes the given node from the list.
void list_delete(list_t* list, list_node_t* node);

#ifdef __cplusplus
}
#endif

#endif
