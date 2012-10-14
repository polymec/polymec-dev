#ifndef ARBI_AVL_TREE_H
#define ARBI_AVL_TREE_H

#include <stdlib.h>
#include "core/arbi.h"
#include "core/comparators.h"
#include "arena/proto.h"

// An AVL tree is a balanced binary tree that can be used to implement other data structures.
// One defines tree and node types using
// DEFINE_AVL_TREE(tree_name, element, comparator, destructor)
//
// Interface for a type x_tree_t (with node type x_tree_node_t and datum x) defined with 
// DEFINE_AVL_TREE(x_tree, x, x_comparator, x_destructor):
// 
// x_tree_t* x_tree_new() - Creates a new empty tree.
// void x_tree_free(x_tree_t* tree) - Destroys the tree.
// void x_tree_clear(x_tree_t* tree) - Empties the tree.
// x_tree_node_t* x_tree_find(x_tree_t* tree, x datum) - Finds the node containing the datum.
// void x_tree_insert(x_tree_t* tree, x datum) - Inserts a datum into the tree.
// void x_tree_delete(x_tree_t* tree, x_tree_node_t* node) - Deletes the node from the tree.
// void x_tree_node_visit(x_tree_node_t* node, x_tree_node_visitor visit, void*) - Visits the given node and its subtree.
// int x_tree_size(x_tree_t* tree) - Returns the number of nodes in the tree. Computed, not stored.

// This maximum function is used to keep AVL trees ordered.
static int avl_tree_max(int x, int y)
{
  return ((x > y) ? x : y);
}

#define DEFINE_AVL_TREE(tree_name, element, comparator, destructor) \
typedef struct tree_name##_node_t tree_name##_node_t; \
typedef void (*tree_name##_node_visitor)(tree_name##_node_t*, void*); \
struct tree_name##_node_t \
{ \
  tree_name##_node_t* left; \
  tree_name##_node_t* right; \
  int depth; \
  element attribute; \
}; \
\
typedef struct tree_name##_t tree_name##_t; \
typedef void (*tree_name##_dtor)(element); \
struct tree_name##_t \
{ \
  tree_name##_node_t* root; \
  ARENA* arena; \
  tree_name##_dtor dtor; \
}; \
\
static inline tree_name##_t* tree_name##_new() \
{ \
  tree_name##_t* tree = malloc(sizeof(tree_name##_t)); \
  tree->arena = NULL; \
  tree->root = NULL; \
  tree->dtor = destructor; \
  return tree; \
} \
\
static inline void tree_name##_clear_node(tree_name##_node_t* node, tree_name##_dtor dtor) \
{ \
  if (node != NULL) \
  { \
    tree_name##_clear_node(node->left, dtor); \
    tree_name##_clear_node(node->right, dtor); \
    if (dtor != NULL) \
      dtor(node->attribute); \
    free(node); \
  } \
} \
\
static inline void tree_name##_clear(tree_name##_t* tree) \
{ \
  tree_name##_clear_node(tree->root, tree->dtor); \
  tree->root = NULL; \
} \
\
static inline void tree_name##_free(tree_name##_t* tree) \
{ \
  tree_name##_clear(tree); \
  free(tree); \
} \
\
static inline tree_name##_node_t* tree_name##_find_node(tree_name##_node_t* node, element datum) \
{ \
  if (node == NULL) \
    return NULL; \
  int result = comparator(datum, node->attribute); \
  if (result == 0) \
    return node; \
  else if (result < 0) \
    return tree_name##_find_node(node->left, datum); \
  else \
    return tree_name##_find_node(node->right, datum); \
} \
\
static inline tree_name##_node_t* tree_name##_find(tree_name##_t* tree, element datum) \
{ \
  return tree_name##_find_node(tree->root, datum); \
} \
\
static inline int tree_name##_node_depth(tree_name##_node_t* node) \
{ \
  if (node == NULL) \
    return -1; \
  else return node->depth; \
} \
\
static inline tree_name##_node_t* tree_name##_node_single_rotate_with_left(tree_name##_node_t* node) \
{ \
  tree_name##_node_t* n = node->left; \
  ASSERT(n != NULL); \
  node->left = n->right; \
  n->right = node; \
  node->depth = avl_tree_max(tree_name##_node_depth(node->left), tree_name##_node_depth(node->right)) + 1; \
  n->depth = avl_tree_max(tree_name##_node_depth(n->left), node->depth) + 1; \
  return n; \
} \
\
static inline tree_name##_node_t* tree_name##_node_single_rotate_with_right(tree_name##_node_t* node) \
{ \
  tree_name##_node_t* n = node->right; \
  ASSERT(n != NULL); \
  node->right = n->left; \
  n->left = node; \
  node->depth = avl_tree_max(tree_name##_node_depth(node->left), tree_name##_node_depth(node->right)) + 1; \
  n->depth = avl_tree_max(tree_name##_node_depth(n->right), node->depth) + 1; \
  return n; \
} \
\
static inline tree_name##_node_t* tree_name##_node_double_rotate_with_left(tree_name##_node_t* node) \
{ \
  node->left = tree_name##_node_single_rotate_with_right(node->left); \
  return tree_name##_node_single_rotate_with_left(node); \
} \
\
static inline tree_name##_node_t* tree_name##_node_double_rotate_with_right(tree_name##_node_t* node) \
{ \
  node->right = tree_name##_node_single_rotate_with_left(node->right); \
  return tree_name##_node_single_rotate_with_right(node); \
} \
\
static inline tree_name##_node_t* tree_name##_insert_node(tree_name##_node_t* node, element datum) \
{ \
  if (node == NULL) \
  { \
    node = malloc(sizeof(tree_name##_node_t)); \
    node->left = NULL; \
    node->right = NULL; \
    node->attribute = datum; \
    node->depth = 0; \
    return node; \
  } \
  else \
  { \
    int result = comparator(datum, node->attribute); \
    if (result < 0) \
    { \
      node->left = tree_name##_insert_node(node->left, datum); \
      if ((tree_name##_node_depth(node->left) - tree_name##_node_depth(node->right)) == 2) \
      { \
        int result2 = comparator(datum, node->left->attribute); \
        if (result2 < 0) \
          node = tree_name##_node_single_rotate_with_left(node); \
        else \
          node = tree_name##_node_double_rotate_with_left(node); \
      } \
    } \
    else if (result > 0) \
    { \
      node->right = tree_name##_insert_node(node->right, datum); \
      if ((tree_name##_node_depth(node->right) - tree_name##_node_depth(node->left)) == 2) \
      { \
        int result2 = comparator(datum, node->right->attribute); \
        if (result2 > 0) \
          node = tree_name##_node_single_rotate_with_right(node); \
        else \
          node = tree_name##_node_double_rotate_with_right(node); \
      } \
    } \
    node->depth = avl_tree_max(tree_name##_node_depth(node->left), tree_name##_node_depth(node->right)) + 1; \
    return node; \
  } \
} \
\
static inline void tree_name##_insert(tree_name##_t* tree, element datum) \
{ \
  tree->root = tree_name##_insert_node(tree->root, datum); \
} \
\
static inline tree_name##_node_t* tree_name##_find_node_parent(tree_name##_node_t* root, \
                                             tree_name##_node_t* node) \
{ \
  ASSERT(node != NULL); \
  int result = comparator(node->attribute, root->attribute); \
  if (result == 0) \
  { \
    return NULL; \
  } \
  else if (result < 0) \
  { \
    if (root->left == node) \
      return root; \
    else \
      return tree_name##_find_node_parent(root->left, node); \
  } \
  else \
  { \
    if (root->right == node) \
      return root; \
    else \
      return tree_name##_find_node_parent(root->right, node); \
  } \
} \
\
static inline void tree_name##_delete(tree_name##_t* tree, tree_name##_node_t* node) \
{ \
  ASSERT(node != NULL); \
  tree_name##_node_t* parent = tree_name##_find_node_parent(tree->root, node); \
  if ((node->left == NULL) || (node->right == NULL)) \
  { \
    if ((node->left == NULL) && (node->right == NULL)) \
    { \
      if (parent != NULL) \
      { \
        if (parent->left == node) \
          parent->left = NULL; \
        else \
          parent->right = NULL; \
      } \
    } \
    else \
    { \
      tree_name##_node_t* child = (node->left != NULL) ? node->left : node->right; \
      if (parent != NULL) \
      { \
        if (parent->left == node) \
          parent->left = child; \
        else \
          parent->right = child; \
      } \
    } \
  } \
  else \
  { \
    tree_name##_node_t* pred = node->left; \
    tree_name##_node_t* pred_parent = node; \
    tree_name##_node_t* succ = node->right; \
    tree_name##_node_t* succ_parent = node; \
    while (pred->right != NULL) \
    { \
      pred_parent = pred; \
      pred = pred->right; \
    } \
    while (succ->left != NULL) \
    { \
      succ_parent = succ; \
      succ = succ->left; \
    } \
    tree_name##_node_t* replacement; \
    tree_name##_node_t* replacement_parent; \
    int imbalance = succ->depth - pred->depth; \
    if (imbalance < 0) \
    { \
      replacement = succ; \
      replacement_parent = succ_parent; \
    } \
    else \
    { \
      replacement = pred; \
      replacement_parent = pred_parent; \
    } \
    if (replacement->left != NULL) \
    { \
      replacement_parent->right = replacement->left; \
    } \
    parent->left = replacement; \
  } \
  if (tree->dtor != NULL) \
    tree->dtor(node->attribute); \
  free(node); \
} \
\
static inline void tree_name##_node_visit(tree_name##_node_t* node, tree_name##_node_visitor visit, void* arg) \
{ \
  if (node == NULL) return; \
  if (node->left != NULL) \
    tree_name##_node_visit(node->left, visit, arg); \
  visit(node, arg); \
  if (node->right != NULL) \
    tree_name##_node_visit(node->right, visit, arg); \
} \
\
static inline void tree_name##_count_node(tree_name##_node_t* node, void* arg) \
{ \
  int* count = (int*)arg; \
  (*count)++; \
} \
\
static inline int tree_name##_size(tree_name##_t* tree) \
{ \
  int count = 0; \
  tree_name##_node_visit(tree->root, &tree_name##_count_node, (void*)&count); \
  return count; \
} \
\


// Define some avl_trees.
DEFINE_AVL_TREE(int_tree, int, int_cmp, NULL)
DEFINE_AVL_TREE(double_tree, double, double_cmp, NULL)


// Old stuff.

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
void        avl_tree_delete(avl_tree_t* tree, avl_node_t* node);
avl_node_t* avl_tree_root(avl_tree_t* tree);

// Visit a node and its subtree.
typedef void (*avl_node_visitor)(avl_node_t*, void*);
void        avl_node_visit(avl_node_t* node, avl_node_visitor visit, void*);

// Some specialized AVL tree constructors.
avl_tree_t* int_avl_tree_new();

#ifdef __cplusplus
}
#endif

#endif
