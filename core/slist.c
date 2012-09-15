#include "core/slist.h"
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

struct slist_t
{
  slist_node_t* front;
  slist_node_t* back;
  int size;
  slist_value_dtor dtor;
};

slist_t* slist_new(slist_value_dtor dtor)
{
  slist_t* list = malloc(sizeof(slist_t));
  list->front = list->back = NULL;
  list->size = 0;
  list->dtor = dtor;
  return list;
}

void slist_free(slist_t* list)
{
  slist_node_t* n;
  while (list->front != NULL)
  {
    n = list->front;
    list->front = n->next;
    if (list->dtor)
      list->dtor(n->value);
    free(n);
  }
}

slist_node_t* slist_front(slist_t* list)
{
  return list->front;
}

slist_node_t* slist_back(slist_t* list)
{
  return list->back;
}

int slist_size(slist_t* list)
{
  return list->size;
}

slist_node_t* slist_find(slist_t* list, void* value, slist_cmp comparator)
{
  slist_node_t* n = list->front;
  while ((n != NULL) && (comparator(n->value, value) != 0))
    n = n->next;
  return n;
}

void slist_insert(slist_t* list, void* value, slist_node_t* node)
{
  ASSERT(node != NULL);

  slist_node_t* n = malloc(sizeof(slist_node_t));
  n->value = value;
  n->next = NULL;

  if (list->front == NULL)
  {
    list->front = n;
    list->back = n;
    list->size = 1;
  }
  else if (list->front == list->back)
  {
    ASSERT(node == list->front);
    n->next = list->front;
    list->front = n;
    list->back = n->next;
    list->size = 2;
  }
  else
  {
    slist_node_t* p = list->front;
    while (p->next != node)
    {
      ASSERT(p != NULL); // make sure node is in this list!
      p = p->next;
    }
    n->next = p->next;
    p->next = n;
    list->size += 1;
  }
}

void slist_append(slist_t* list, void* value)
{
  slist_node_t* n = malloc(sizeof(slist_node_t));
  n->value = value;
  n->next = NULL;
  if (list->back == NULL)
  {
    ASSERT(list->front == NULL);
    list->front = n;
    list->back = n;
  }
  else
  {
    list->back->next = n;
    list->back = n;
  }
}

void* slist_pop(slist_t* list)
{
  ASSERT(list->size > 0);
  slist_node_t* node = list->front;
  list->front = list->front->next;
  if (list->size == 1)
    list->back = NULL;
  list->size -= 1;
  void* val = node->value;
  free(node);
  return val;
}

void slist_remove(slist_t* list, slist_node_t* node)
{
  if (list->front == node)
    slist_pop(list);
  slist_node_t* p = list->front;
  while ((p->next != node) && (p != NULL))
    p = p->next;
  if (p != NULL)
  {
    p->next = node->next;
    if (list->dtor)
      list->dtor(node->value);
    free(node);
    list->size -= 1;
  }
}

#ifdef __cplusplus
}
#endif

