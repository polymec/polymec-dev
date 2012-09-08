#include "core/hashmap.h"
#include "core/list.h"

#ifdef __cplusplus
extern "C" {
#endif

struct hashmap_t
{
  tommy_hashlin map;
  list_t* keys;
  hashmap_comp_func comp;
  hashmap_hash_func hash;
  hashmap_element_dtor dtor;
};

typedef struct 
{
  tommy_node node;
  void* key;
  void* value;
} hashmap_elem_t;

hashmap_t* hashmap_new(hashmap_comp_func comp, hashmap_element_dtor dtor)
{
  ASSERT(comp != NULL);

  hashmap_t* h = malloc(sizeof(hashmap_t));
  tommy_hashlin_init(&h->map);
  h->keys = list_new();
  h->comp = comp;
  h->hash = &tommy_hash_u32;
  h->dtor = dtor;
  return h;
}

hashmap_t* hashmap_new_with_hash_func(hashmap_comp_func comp, hashmap_element_dtor dtor, hashmap_hash_func hash)
{
  ASSERT(comp != NULL);
  ASSERT(hash != NULL);

  hashmap_t* h = malloc(sizeof(hashmap_t));
  tommy_hashlin_init(&h->map);
  h->comp = comp;
  h->hash = hash;
  h->dtor = dtor;
  return h;
}

void hashmap_free(hashmap_t* hashmap)
{
  hashmap_clear(hashmap);
  tommy_hashlin_done(&h->map);
  free(hashmap);
}

void hashmap_clear(hashmap_t* hashmap)
{
  list_node_t* i = list_head(hashmap->keys);
  while (i != NULL)
  {
    void* key = list_node_value(i);
    hashmap_delete(hashmap, key);
    i = i->next;
  }
  list_clear(hashmap->keys);
}

void hashmap_insert(hashmap_t* hashmap, void* key, void* value)
{
  hashmap_elem_t* elem = malloc(sizeof(hashmap_elem_t));
  elem->value = value;
  tommy_hashlin_insert(&hashmap->map, &elem->node, hashmap->hash(elem->value));
}

void hashmap_delete(hashmap_t* hashmap, void* key)
{
  hashmap_elem_t* elem = tommy_hashlin_remove(&hashmap->map, hashmap->comp, key, hashmap->hash(key))
  if (elem != NULL)
  {
    if (hashmap->dtor != NULL)
      hashmap->dtor(elem->value);
    free(elem);
  }
}

void* hashmap_find(hashmap_t* hashmap, void* key)
{
  hashmap_elem_t* elem = tommy_hashlin_search(&hashmap->map, hashmap->comp, key, hashmap->hash(key))
  if (elem != NULL)
    return elem->value;
  else
    return NULL;
}

void* hashmap_find_and_delete(hashmap_t* hashmap, void* key)
{
  hashmap_elem_t* elem = tommy_hashlin_remove(&hashmap->map, hashmap->comp, key, hashmap->hash(key))
  if (elem != NULL)
  {
    void* value = elem->value;
    free(elem);
    return value;
  }
  else
    return NULL;
}

int hashmap_size(hashmap_t* hashmap)
{
  return list_size(hashmap_keys);
}

list_t* hashmap_keys(hashmap_t* hashmap)
{
}

#ifdef __cplusplus
}
#endif

