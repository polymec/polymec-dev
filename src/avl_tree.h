#ifndef ARBI_AVL_TREE_H
#define ARBI_AVL_TREE_H

#include "arbi.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct avl_node_t avl_node_t;
struct avl_node_t
{
  avl_node_t* left;
  avl_node_t* right;
  int depth;
  void* attribute;
};

typedef struct avl_tree_t avl_tree_t;
typedef int (*avl_tree_attribute_cmp)(void*, void*);
typedef void (*avl_tree_attribute_dtor)(void*);

avl_tree_t* avl_tree_new(avl_tree_attribute_cmp cmp, avl_tree_attribute_dtor dtor);
void        avl_tree_free(avl_tree_t* tree);
void        avl_tree_clear(avl_tree_t* tree);
avl_node_t* avl_tree_find(avl_tree_t* tree, void* datum);
void        avl_tree_insert(avl_tree_t* tree, void* datum);
void        avl_tree_delete(avl_tree_t* tree, void* datum);

#ifdef __cplusplus
}
#endif

#endif
