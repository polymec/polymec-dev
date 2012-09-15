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

avl_tree_t* avl_tree_new(avl_tree_attribute_cmp cmp, avl_tree_attribute_dtor dtor)
{
  ASSERT(cmp != NULL);
  avl_tree_t* tree = malloc(sizeof(avl_tree_t));
  tree->root = NULL;
  tree->cmp = cmp;
  tree->dtor = dtor;
  return tree;
}

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

void avl_tree_free(avl_tree_t* tree)
{
  avl_tree_clear(tree);
  free(tree);
}

void avl_tree_clear(avl_tree_t* tree)
{
  avl_tree_clear_node(tree->root, tree->dtor);
  tree->root = NULL;
}

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

avl_node_t* avl_tree_find(avl_tree_t* tree, void* datum)
{
  return avl_tree_find_node(tree->root, datum, tree->cmp);
}

static inline int avl_node_depth(avl_node_t* node)
{
  if (node == NULL)
    return -1;
  else return node->depth;
}

static avl_node_t* single_rotate_with_left(avl_node_t* node)
{
  avl_node_t* n = node->left;
  ASSERT(n != NULL);
  node->left = n->right;
  n->right = node;

  node->depth = MAX(avl_node_depth(node->left), avl_node_depth(node->right)) + 1;
  n->depth = MAX(avl_node_depth(n->left), node->depth) + 1;
  return n;
}

static avl_node_t* single_rotate_with_right(avl_node_t* node)
{
  avl_node_t* n = node->right;
  ASSERT(n != NULL);
  node->right = n->left;
  n->left = node;

  node->depth = MAX(avl_node_depth(node->left), avl_node_depth(node->right)) + 1;
  n->depth = MAX(avl_node_depth(n->right), node->depth) + 1;
  return n;
}

static avl_node_t* double_rotate_with_left(avl_node_t* node)
{
  node->left = single_rotate_with_right(node->left);
  return single_rotate_with_left(node);
}

static avl_node_t* double_rotate_with_right(avl_node_t* node)
{
  node->right = single_rotate_with_left(node->right);
  return single_rotate_with_right(node);
}

static avl_node_t* avl_tree_insert_node(avl_node_t* node, void* datum, 
                                        avl_tree_attribute_cmp cmp)
{
  if (node == NULL)
  {
    node = malloc(sizeof(avl_node_t));
    node->left = NULL;
    node->right = NULL;
    node->attribute = datum;
    node->depth = 0;
    return node;
  }
  else 
  {
    int result = cmp(datum, node->attribute);
    if (result < 0)
    {
      node->left = avl_tree_insert_node(node->left, datum, cmp);
      if ((avl_node_depth(node->left) - avl_node_depth(node->right)) == 2)
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
      node->right = avl_tree_insert_node(node->right, datum, cmp);

      if ((avl_node_depth(node->right) - avl_node_depth(node->left)) == 2)
      {
        int result2 = cmp(datum, node->right->attribute);
        if (result2 > 0)
          node = single_rotate_with_right(node);
        else
          node = double_rotate_with_right(node);
      }
    }
    // Otherwise the tree already contains this datum -- do nothing.
    node->depth = MAX(avl_node_depth(node->left), avl_node_depth(node->right)) + 1;
    return node;
  }
}

void avl_tree_insert(avl_tree_t* tree, void* datum)
{
  tree->root = avl_tree_insert_node(tree->root, datum, tree->cmp);
}

static avl_node_t* avl_tree_find_node_parent(avl_node_t* root, 
                                             avl_node_t* node,
                                             avl_tree_attribute_cmp cmp)
{
  ASSERT(node != NULL);
  int result = cmp(node->attribute, root->attribute);
  if (result == 0)
  {
    return NULL;
  }
  else if (result < 0)
  {
    if (root->left == node) 
      return root;
    else 
      return avl_tree_find_node_parent(root->left, node, cmp);
  }
  else
  {
    if (root->right == node)
      return root;
    else
      return avl_tree_find_node(root->right, node, cmp);
  }
}

void avl_tree_delete(avl_tree_t* tree, avl_node_t* node)
{
  ASSERT(node != NULL);

  // Find the parent node.
  avl_node_t* parent = avl_tree_find_node_parent(tree->root, node, tree->cmp);

  if ((node->left == NULL) || (node->right == NULL))
  {
    if ((node->left == NULL) && (node->right == NULL)) // leaf!
    {
      // Unhook the node.
      if (parent != NULL)
      {
        if (parent->left == node)
          parent->left = NULL;
        else
          parent->right = NULL;
      }
    }
    else 
    {
      avl_node_t* child = (node->left != NULL) ? node->left : node->right;
      // Rewire the tree.
      if (parent != NULL)
      {
        if (parent->left == node)
          parent->left = child;
        else
          parent->right = child;
      }
    }
  }
  else
  {
    // We must do some gymnastics, since the node has two children.
    // Replace the node with either its in-order predecessor (pred) or 
    // successor (succ) depending on balancing.
    avl_node_t* pred = node->left;
    avl_node_t* pred_parent = node;
    avl_node_t* succ = node->right;
    avl_node_t* succ_parent = node;
    while (pred->right != NULL)
    {
      pred_parent = pred;
      pred = pred->right;
    }
    while (succ->left != NULL)
    {
      succ_parent = succ;
      succ = succ->left;
    }

    // Which one do we use as a replacement?
    avl_node_t* replacement;
    avl_node_t* replacement_parent;
    int imbalance = succ->depth - pred->depth;
    if (imbalance < 0) // left-heavy 
    {
      replacement = succ;
      replacement_parent = succ_parent;
    }
    else // right-heavy (or balanced)
    {
      replacement = pred;
      replacement_parent = pred_parent;
    }

    // Replace the node.
    if (replacement->left != NULL)
    {
      replacement_parent->right = replacement->left;
    }
    parent->left = replacement;

  }

  // Kill the node. Nobody loves you, node!
  if (tree->dtor)
    tree->dtor(node->attribute);
  free(node);
}

void avl_node_visit(avl_node_t* node, avl_node_visitor visit, void* arg)
{
  if (node == NULL) return;
  if (node->left != NULL)
    avl_node_visit(node->left, visit, arg);
  visit(node, arg);
  if (node->right != NULL)
    avl_node_visit(node->right, visit, arg);
}

static int avl_tree_int_cmp(void* left, void* right)
{
  int l = (int)left, r = (int)right;
  return (l < r) ? -1 : (l > r) ? 1 : 0;
}

avl_node_t* avl_tree_root(avl_tree_t* tree)
{
  return tree->root;
}

avl_tree_t* int_avl_tree_new()
{
  return avl_tree_new(&avl_tree_int_cmp, NULL);
}

#ifdef __cplusplus
}
#endif

