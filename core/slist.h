#ifndef ARBI_SLIST_H
#define ARBI_SLIST_H

#include "arbi.h"

#ifdef __cplusplus
extern "C" {
#endif

// slist_t is a singly-linked list class that stores homogeneous types.
typedef struct slist_node_t slist_node_t;
struct slist_node_t
{
  void* value;
  slist_node_t* next;
};

typedef struct slist_t slist_t;
typedef int (*slist_cmp)(void*, void*);
typedef void (*slist_value_dtor)(void*);

slist_t* slist_new(slist_value_dtor dtor);
void slist_free(slist_t* list);
slist_node_t* slist_front(slist_t* list);
slist_node_t* slist_back(slist_t* list);
int slist_size(slist_t* list);
slist_node_t* slist_find(slist_t* list, void* value, slist_cmp comparator);
void slist_insert(slist_t* list, void* value, slist_node_t* node);
void slist_append(slist_t* list, void* value);
void* slist_pop(slist_t* list);
void slist_remove(slist_t* list, slist_node_t* node);

#ifdef __cplusplus
}
#endif

#endif
