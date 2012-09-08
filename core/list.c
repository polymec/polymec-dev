#include "core/list.h"

#ifdef __cplusplus
extern "C" {
#endif

struct list_t
{
  tommy_list list;
  int size;
  list_element_dtor dtor;
};

typedef struct 
{
  tommy_node node;
  void* value;
} list_element_t;

list_t* list_new(list_element_dtor dtor)
{
  list_t* l = malloc(sizeof(list_t));
  tommy_list_init(&l->list);
  l->size = 0;
  l->dtor = dtor;
  return l;
}

void list_free(list_t* list)
{
  list_clear(list);
  free(list);
}

void list_clear(list_t* list)
{
  tommy_node* i = tommy_list_head(&list->list);
  while (i != NULL)
  {
    tommy_node* i_next = i->next;
    list_element_t* elem = i->data;
    if (list->dtor != NULL)
      list->dtor(elem->value);
    free(elem);
    i = i_next;
  }
  list->size = 0;
}

list_node_t* list_head(list_t* list)
{
  return tommy_list_head(&list->list);
}

list_node_t* list_tail(list_t* list)
{
  return tommy_list_tail(&list->list);
}

void* list_node_value(list_node_t* node)
{
  list_element_t* elem = i->data;
  return elem->value;
}

int list_size(list_t* list)
{
  return list->size;
}

void list_sort(list_t* list, list_comp_func comparator)
{
  tommy_list_sort(&list->list, comparator);
}

void list_prepend(list_t* list, void* value)
{
  list_element_t* node = malloc(sizeof(list_element_t));
  node->value = value;
  tommy_list_insert_head(&list->list, &node->node, node);
}

void list_append(list_t* list, void* value)
{
  list_element_t* node = malloc(sizeof(list_element_t));
  node->value = value;
  tommy_list_insert_head(&list->list, &node->node, node);
}

void list_insert_before(list_t* list, void* value, list_node_t* node)
{
  tommy_node* i = tommy_list_head(&list->list);
  while ((i != NULL) && (i != node))
    i = i->next;
  if (i == NULL) return;
  list_element_t* n = malloc(sizeof(list_element_t));
  n->value = value;
  tommy_list_insert_head(&i, &n->node, n);
}

list_node_t* list_find(list_t* list, void* value, list_comp_func comp)
{
  while (i != NULL)
  {
    list_element_t* elem = i->data;
    if (comp(value, elem->value) == 0)
      return i;
    i = i->next;
  }
  return NULL;
}

void list_delete(list_t* list, list_node_t* node)
{
  tommy_node* i = tommy_list_head(&list->list);
  while ((i != NULL) && (i != node))
    i = i->next;
  if (i == NULL) return;
  tommy_list_remove_existing(&list->list, node);

  // Deallocate the node.
  list_element_t* elem = i->data;
  if (list->dtor != NULL)
    list->dtor(elem->value);
  free(elem);
}

#ifdef __cplusplus
}
#endif

