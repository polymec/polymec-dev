#include <stdlib.h>
#include "avl_tree.h"

#ifdef MAX
#undef MAX
#endif
#define MAX(a, b) ((a > b) ? a : b)

#ifdef __cplusplus
extern "C" {
#endif

struct avl_tree_t
{
  avl_node_t* root;
  avl_tree_attribute_cmp cmp;
  avl_tree_attribute_dtor dtor;
};

//------------------------------------------------------------------------
avl_tree_t* avl_tree_new(avl_tree_attribute_cmp cmp, avl_tree_attribute_dtor dtor)
{
  ASSERT(cmp != NULL);
  avl_tree_t* tree = malloc(sizeof(avl_tree_t));
  tree->root = NULL;
  tree->cmp = cmp;
  tree->dtor = dtor;
  return tree;
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
static void avl_tree_clear_node(avl_node_t* node, avl_tree_attribute_dtor dtor)
{
  if (node != NULL)
  {
    avl_tree_clear_node(node->left, dtor);
    avl_tree_clear_node(node->right, dtor);
    if (dtor != NULL)
      dtor(node->attribute);
    free(node);
  }
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
void avl_tree_free(avl_tree_t* tree)
{
  avl_tree_clear(tree);
  free(tree);
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
void avl_tree_clear(avl_tree_t* tree)
{
  avl_tree_clear_node(tree->root, tree->dtor);
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
static avl_node_t* avl_tree_find_node(avl_node_t* node, void* datum, avl_tree_attribute_cmp cmp)
{
  if (node == NULL)
    return NULL;
  int result = cmp(datum, node->attribute);
  if (result == 0)
    return node;
  else if (result < 0)
    return avl_tree_find_node(node->left, datum, cmp);
  else
    return avl_tree_find_node(node->right, datum, cmp);
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
avl_node_t* avl_tree_find(avl_tree_t* tree, void* datum)
{
  return avl_tree_find_node(tree->root, datum, tree->cmp);
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
static avl_node_t* single_rotate_with_left(avl_node_t* node)
{
  avl_node_t* n = node->left;
  node->left = n->right;
  n->right = node;

  node->depth = MAX(node->left->depth, node->right->depth) + 1;
  n->depth = MAX(n->left->depth, node->depth) + 1;
  return n;
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
static avl_node_t* single_rotate_with_right(avl_node_t* node)
{
  avl_node_t* n = node->right;
  node->right = n->left;
  n->left = node;

  node->depth = MAX(node->left->depth, node->right->depth) + 1;
  n->depth = MAX(n->right->depth, node->depth) + 1;
  return n;
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
static avl_node_t* double_rotate_with_left(avl_node_t* node)
{
  node->left = single_rotate_with_right(node->left);
  return single_rotate_with_left(node);
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
static avl_node_t* double_rotate_with_right(avl_node_t* node)
{
  node->right = single_rotate_with_left(node->right);
  return single_rotate_with_right(node);
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
static avl_node_t* avl_tree_insert_node(avl_node_t* node, void* datum, 
                                        avl_tree_attribute_cmp cmp, int depth)
{
  if (node == NULL)
  {
    node = malloc(sizeof(avl_node_t));
    node->left = NULL;
    node->right = NULL;
    node->attribute = datum;
    node->depth = depth;
    return node;
  }
  else 
  {
    int result = cmp(datum, node->attribute);
    if (result < 0)
    {
      node->left = avl_tree_insert_node(node->left, datum, cmp, depth+1);
      if ((node->left->depth - node->right->depth) == 2)
      {
        int result2 = cmp(datum, node->left->attribute);
        if (result2 < 0)
          node = single_rotate_with_left(node);
        else
          node = double_rotate_with_left(node);
      }
    }
    else if (result > 0)
    {
      node->right = avl_tree_insert_node(node->right, datum, cmp, depth+1);
      if ((node->right->depth - node->left->depth) == 2)
      {
        int result2 = cmp(datum, node->right->attribute);
        if (result2 < 0)
          node = single_rotate_with_right(node);
        else
          node = double_rotate_with_right(node);
      }
    }
    // Otherwise the tree already contains this datum -- do nothing.
    return node;
  }
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
void avl_tree_insert(avl_tree_t* tree, void* datum)
{
  avl_tree_insert_node(tree, datum, tree->cmp, 0);
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
void avl_tree_delete(avl_tree_t* tree, void* datum)
{
  // FIXME: Unimplemented for now!
}
//------------------------------------------------------------------------

#ifdef __cplusplus
}
#endif

