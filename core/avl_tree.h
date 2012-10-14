#ifndef ARBI_AVL_TREE_H
#define ARBI_AVL_TREE_H

#include "core/arbi.h"
#include "core/comparators.h"
#include "arena/proto.h"

#ifdef MAX(x,y)
#undef MAX
#define MAX(x,y) ((x > y) ? x : y)
#endif

#define DEFINE_AVL_TREE(tree_name, tree_node_name, element, comparator, destructor) \
typedef struct tree_node_name##_t tree_node_name##_t; \
struct tree_node_name##_t \
{ \
  tree_node_name##_t* left; \
  tree_node_name##_t* right; \
  int depth; \
  element attribute; \
}; \
\
typedef struct tree_name##_t tree_name##_t; \
typedef void (*tree_name##_dtor)(element); \
struct tree_name##_t \
{ \
  tree_node_name##_t* root; \
  ARENA* arena; \
  tree_name##_dtor dtor; \
}; \
\
tree_name##_t* tree_name##_new() \
{ \
  tree_name##_t* tree = malloc(sizeof(tree_name##_t)); \
  tree->arena = NULL; \
  tree->root = NULL; \
  tree->dtor = destructor; \
  return tree; \
} \
\
static inline void tree_name##_clear_node(tree_node_name##_t* node, tree_name##_dtor dtor) \
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
static inline tree_node_name##_t* tree_name##_find_node(tree_node_name##_t* node, element datum) \
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
static inline tree_node_name##_t* tree_name##_find(tree_name##_t* tree, element datum) \
{ \
  return tree_name##_find_node(tree->root, datum); \
} \
\
static inline int tree_node_name##_depth(tree_node_name##_t* node) \
{ \
  if (node == NULL) \
    return -1; \
  else return node->depth; \
} \
\
static inline tree_node_name##_t* tree_node_name##_single_rotate_with_left(tree_node_name##_t* node) \
{ \
  tree_node_name##_t* n = node->left; \
  ASSERT(n != NULL); \
  node->left = n->right; \
  n->right = node; \
  node->depth = MAX(tree_node_name##_depth(node->left), tree_node_name##_depth(node->right)) + 1; \
  n->depth = MAX(tree_node_name##_depth(n->left), node->depth) + 1; \
  return n; \
} \
\
static inline tree_node_name##_t* tree_node_name##_single_rotate_with_right(tree_node_name##_t* node) \
{ \
  tree_node_name##_t* n = node->right; \
  ASSERT(n != NULL); \
  node->right = n->left; \
  n->left = node; \
  node->depth = MAX(tree_node_name##_depth(node->left), tree_node_name##_depth(node->right)) + 1; \
  n->depth = MAX(tree_node_name##_depth(n->right), node->depth) + 1; \
  return n; \
} \
\
static inline tree_node_name##_t* tree_node_name##double_rotate_with_left(tree_node_name##_t* node) \
{ \
  node->left = tree_node_name##single_rotate_with_right(node->left); \
  return tree_node_name##single_rotate_with_left(node); \
} \
\
static inline tree_node_name##_t* tree_node_name##double_rotate_with_right(tree_node_name##_t* node) \
{ \
  node->right = tree_node_name##single_rotate_with_left(node->right); \
  return tree_node_name##single_rotate_with_right(node); \
} \
\
static inline tree_node_name##_t* tree_name##_insert_node(tree_node_name##_t* node, element datum) \
{ \
  if (node == NULL) \
  { \
    node = malloc(sizeof(tree_node_name##_t)); \
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
      if ((tree_node_name##_depth(node->left) - tree_node_name##_depth(node->right)) == 2) \
      { \
        int result2 = comparator(datum, node->left->attribute); \
        if (result2 < 0) \
          node = tree_node_name##single_rotate_with_left(node); \
        else \
          node = tree_node_name##double_rotate_with_left(node); \
      } \
    } \
    else if (result > 0) \
    { \
      node->right = tree_name##_insert_node(node->right, datum); \
      if ((tree_node_name##_depth(node->right) - tree_node_name##_depth(node->left)) == 2) \
      { \
        int result2 = comparator(datum, node->right->attribute); \
        if (result2 > 0) \
          node = tree_node_name##single_rotate_with_right(node); \
        else \
          node = tree_node_name##double_rotate_with_right(node); \
      } \
    } \
    node->depth = MAX(tree_node_name##_depth(node->left), tree_node_name##_depth(node->right)) + 1; \
    return node; \
  } \
} \
\
static inline void tree_name##_insert(tree_name##_t* tree, element datum) \
{ \
  tree->root = tree_name##_insert_node(tree->root, datum); \
} \
\
static inline tree_node_name##_t* tree_name##_find_node_parent(tree_node_name##_t* root, \
                                             tree_node_name##_t* node) \
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
      return tree_name##_find_node(root->right, node); \
  } \
} \
\
static inline void tree_name##_delete(tree_name##_t* tree, tree_node_name##_t* node) \
{ \
  ASSERT(node != NULL); \
  tree_node_name##_t* parent = tree_name##_find_node_parent(tree->root, node); \
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
      tree_node_name##_t* child = (node->left != NULL) ? node->left : node->right; \
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
    tree_node_name##_t* pred = node->left; \
    tree_node_name##_t* pred_parent = node; \
    tree_node_name##_t* succ = node->right; \
    tree_node_name##_t* succ_parent = node; \
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
    tree_node_name##_t* replacement; \
    tree_node_name##_t* replacement_parent; \
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
  if (dtor != NULL) \
    dtor(node->attribute); \
  free(node); \
} \
\
static inline void tree_node_name##_visit(tree_node_name##_t* node, tree_node_name##_visitor visit, element arg) \
{ \
  if (node == NULL) return; \
  if (node->left != NULL) \
    tree_node_name##_visit(node->left, visit, arg); \
  visit(node, arg); \
  if (node->right != NULL) \
    tree_node_name##_visit(node->right, visit, arg); \
} \
\
tree_node_name##_t* tree_name##_root(tree_name##_t* tree) \
{ \
  return tree->root; \
} \


// Define some avl_trees.
DEFINE_AVL_TREE(int_tree, int_tree_node, int, int_cmp, NULL)
DEFINE_AVL_TREE(double_tree, double_tree_node, double, double_cmp, NULL)


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
