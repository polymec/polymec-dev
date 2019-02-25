// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/octree.h"
#include "core/timer.h"

// A bucket PR (point region) octree.

typedef enum
{
  OCTREE_BRANCH_NODE,
  OCTREE_LEAF_NODE
} octree_node_type_t;

typedef struct octree_node_t octree_node_t;

// Branch node.
typedef struct 
{
  int depth;
  octree_node_t* children[8]; 
} octree_branch_node_t;

// Leaf node.
typedef struct
{
  point_t point;
} octree_leaf_node_t;

// Basic octree node -- can be either a branch node or a leaf node.
// We cast between this type and the actual node types below.
struct octree_node_t
{
  octree_node_type_t type;
  int index; // leaf or branch index
  int parent; // index of parent branch (-1 if root)
  union 
  {
    octree_branch_node_t branch_node;
    octree_leaf_node_t leaf_node;
  };
};

struct octree_t 
{
  octree_node_t* root;  // The root node.
  bbox_t bbox; // The bounding box for the tree.
  int num_leaves; // Total number of leaf nodes (points) in the tree.
  int num_branches; // Total number of branch nodes in the tree.
};

static octree_node_t* branch_new(octree_t* tree, int parent, int depth)
{
  octree_node_t* node = polymec_malloc(sizeof(octree_node_t));
  node->type = OCTREE_BRANCH_NODE;
  node->index = tree->num_branches;
  node->parent = parent;
  node->branch_node.depth = depth;
  memset(node->branch_node.children, 0, 8*sizeof(octree_node_t*));
  ++(tree->num_branches);
  return node;
}

// This creates a new leaf node with a single occupant.
static octree_node_t* leaf_new(octree_t* tree, point_t* point, int parent)
{
  octree_node_t* node = polymec_malloc(sizeof(octree_node_t));
  node->type = OCTREE_LEAF_NODE;
  node->index = tree->num_leaves;
  node->parent = parent;
  node->leaf_node.point = *point;
  ++(tree->num_leaves);
  return node;
}

// This returns the slot within a branch into which the given point fits given 
// a rectangular domain centered at the given center point.
static inline int find_slot(point_t* center, point_t* point)
{
  return ((point->x >= center->x) << 2) + 
         ((point->y >= center->y) << 1) + 
          (point->z >= center->z);
}

octree_t* octree_new(bbox_t* bounding_box)
{
  ASSERT(bounding_box->x1 < bounding_box->x2);
  ASSERT(bounding_box->y1 < bounding_box->y2);
  ASSERT(bounding_box->z1 < bounding_box->z2);

  octree_t* tree = polymec_malloc(sizeof(octree_t));
  tree->root = NULL;
  tree->bbox = *bounding_box;
  tree->num_leaves = 0;
  tree->num_branches = 0;
  return tree;
}

void octree_free(octree_t* tree)
{
  octree_clear(tree);
  polymec_free(tree);
}

void octree_insert(octree_t* tree, point_t* point, int index)
{
  START_FUNCTION_TIMER();
  if (tree->root == NULL) // Empty tree
  {
    octree_node_t* node = leaf_new(tree, point, -1);
    tree->root = node;
    STOP_FUNCTION_TIMER();
    return;
  }
  else if (tree->root->type == OCTREE_LEAF_NODE)
  {
    point_t center = {.x = 0.5 * (tree->bbox.x1 + tree->bbox.x2),
                      .y = 0.5 * (tree->bbox.y1 + tree->bbox.y2),
                      .z = 0.5 * (tree->bbox.z1 + tree->bbox.z2)};

    // The tree consists of a single node.
    octree_node_t* root = tree->root;

    // Does the given point already exist here?
    if (points_coincide(&root->leaf_node.point, point)) 
    {
      STOP_FUNCTION_TIMER();
      return;
    }

    // We need to create a branch node here and place the existing 
    // leaf under it.
    octree_node_t* node = root;
    tree->root = branch_new(tree, -1, 0);
    int slot = find_slot(&center, point);
    tree->root->branch_node.children[slot] = node;
    node->parent = tree->root->index;
  }
  
  // Now we proceed with the normal logic, given that the root node 
  // is a branch node.
  ASSERT(tree->root->type == OCTREE_BRANCH_NODE);
  octree_node_t* node = tree->root;
  point_t center = {.x = 0.5 * (tree->bbox.x1 + tree->bbox.x2),
                    .y = 0.5 * (tree->bbox.y1 + tree->bbox.y2),
                    .z = 0.5 * (tree->bbox.z1 + tree->bbox.z2)};
  real_t lx = tree->bbox.x2 - tree->bbox.x1;
  real_t ly = tree->bbox.y2 - tree->bbox.y1;
  real_t lz = tree->bbox.z2 - tree->bbox.z1;
  int slot = find_slot(&center, point);
  static real_t xf[] = {-0.25, -0.25, -0.25, -0.25, +0.25, +0.25, +0.25, +0.25};
  static real_t yf[] = {-0.25, -0.25, +0.25, +0.25, -0.25, -0.25, +0.25, +0.25};
  static real_t zf[] = {-0.25, +0.25, -0.25, +0.25, -0.25, +0.25, -0.25, +0.25};
  int depth = 0;

  // Locate a branch node in the tree on which we can hang a new leaf.
  while ((node->branch_node.children[slot] != NULL) && 
         (node->branch_node.children[slot]->type == OCTREE_BRANCH_NODE))
  {
    node = node->branch_node.children[slot];
    center.x += xf[slot] * lx;
    lx *= 0.5;
    center.y += yf[slot] * ly;
    ly *= 0.5;
    center.z += zf[slot] * lz;
    lz *= 0.5;
    slot = find_slot(&center, point);
    ++depth;
  }

  // Is there a leaf in our spot already?
  octree_node_t* leaf = node->branch_node.children[slot];
  if (leaf == NULL)
  {
    // No leaf here, so we create a new one!
    leaf = leaf_new(tree, point, node->index);
    node->branch_node.children[slot] = leaf;
  }
  else
  {
    // There's a leaf here. Is it identical to the point we're adding?
    if (points_coincide(&leaf->leaf_node.point, point)) 
    {
      STOP_FUNCTION_TIMER();
      return;
    }
    else
    {
      // We have to make a new branch that holds the existing leaf and our 
      // new one. Subdivide our octant until point and leaf->leaf_node.point
      // no longer occupy the same octant (slot).
      int old_point_slot = slot, new_point_slot = slot; 
      do
      {
        // Make a new branch node and point node at it.
        node->branch_node.children[new_point_slot] = branch_new(tree, node->index, depth+1);
        node = node->branch_node.children[new_point_slot];

        // Construct the geometry for the new octant.
        center.x += xf[new_point_slot] * lx;
        lx *= 0.5;
        center.y += yf[new_point_slot] * ly;
        ly *= 0.5;
        center.z += zf[new_point_slot] * lz;
        lz *= 0.5;

        // Compute the slots for the existing and new points.
        new_point_slot = find_slot(&center, point);
        old_point_slot = find_slot(&center, &(leaf->leaf_node.point));
        ++depth;
      }
      while (new_point_slot == old_point_slot);

      // Everything in its right place.
      node->branch_node.children[old_point_slot] = leaf;
      octree_node_t* new_leaf = leaf_new(tree, point, node->index);
      node->branch_node.children[new_point_slot] = new_leaf;
    }
  }
  STOP_FUNCTION_TIMER();
}

int octree_size(octree_t* tree)
{
  return tree->num_leaves;
}

static void node_clear(octree_node_t* node)
{
  if (node == NULL) 
    return;
  else if (node->type == OCTREE_LEAF_NODE)
    polymec_free(node);
  else 
  {
    ASSERT(node->type == OCTREE_BRANCH_NODE);
    for (int i = 0; i < 8; ++i)
      node_clear(node->branch_node.children[i]);
    polymec_free(node);
  }
}

void octree_clear(octree_t* tree)
{
  node_clear(tree->root);
  tree->root = NULL;
  tree->num_leaves = 0;
  tree->num_branches = 0;
}

void* octree_new_leaf_array(octree_t* tree, size_t data_size)
{
  return polymec_calloc(tree->num_leaves, data_size);
}

void* octree_new_branch_array(octree_t* tree, size_t data_size)
{
  return polymec_calloc(tree->num_branches, data_size);
}

static int octree_delete_node(octree_node_t* node, 
                              octree_node_t* parent, 
                              int slot, 
                              int index)
{
  if (node != NULL)
  {
    if (node->type == OCTREE_BRANCH_NODE)
    {
      bool has_children = false;
      int num_deleted = 0;
      for (int s = 0; s < 8; ++s)
      {
        if (num_deleted == 0)
        {
          num_deleted = octree_delete_node(node->branch_node.children[s], 
                                           node, s, index);
        }
        if (node->branch_node.children[s] != NULL)
          has_children = true;
      }
      // Prune this branch node if we need to.
      if (!has_children && (parent != NULL))
      {
        polymec_free(node);
        parent->branch_node.children[slot] = NULL;
      }
      return num_deleted;
    }
    else
    {
      if (node->index == index)
      {
        polymec_free(node);
        if (parent != NULL)
          parent->branch_node.children[slot] = NULL;
        return 1;
      }
    }
  }
  return 0;
}

void octree_delete(octree_t* tree, int index)
{
  START_FUNCTION_TIMER();
  if (tree->root == NULL) 
    return;
  else if (tree->root->type == OCTREE_LEAF_NODE)
  {
    polymec_free(tree->root);
    tree->root = NULL;
    tree->num_leaves = 0;
  }
  else
  {
    // Try to delete the leaf with the given index.
    int num_deleted = octree_delete_node(tree->root, NULL, -1, index);
    tree->num_leaves -= num_deleted;
  }
  STOP_FUNCTION_TIMER();
}

static void octree_visit_pre(octree_node_t* node, 
                             void* context,
                             bool (*visit_branch)(void* context, int depth, int branch_index, int parent_branch_index),
                             void (*visit_leaf)(void* context, int leaf_index, point_t* point, int parent_branch_index))
{
  if (node != NULL)
  {
    if (node->type == OCTREE_BRANCH_NODE)
    {
      for (int slot = 0; slot < 8; ++slot)
        octree_visit_pre(node->branch_node.children[slot], context, visit_branch, visit_leaf);
      if (visit_branch != NULL)
        visit_branch(context, node->branch_node.depth, node->index, node->parent);
    }
    else
    {
      if (visit_leaf != NULL)
        visit_leaf(context, node->index, &(node->leaf_node.point), node->parent);
    }
  }
}

static void octree_visit_post(octree_node_t* node, 
                              void* context,
                              bool (*visit_branch)(void* context, int depth, int branch_index, int parent_branch_index),
                              void (*visit_leaf)(void* context, int leaf_index, point_t* point, int parent_branch_index))
{
  if (node != NULL)
  {
    if (node->type == OCTREE_BRANCH_NODE)
    {
      bool traverse_children = false;
      if (visit_branch != NULL)
        traverse_children = visit_branch(context, node->branch_node.depth, node->index, node->parent);
      if (traverse_children)
      {
        for (int slot = 0; slot < 8; ++slot)
          octree_visit_post(node->branch_node.children[slot], context, visit_branch, visit_leaf);
      }
    }
    else
    {
      if (visit_leaf != NULL)
        visit_leaf(context, node->index, &(node->leaf_node.point), node->parent);
    }
  }
}

void octree_visit(octree_t* tree, 
                  octree_traversal_t order,
                  void* context,
                  bool (*visit_branch)(void* context, int depth, int branch_index, int parent_branch_index),
                  void (*visit_leaf)(void* context, int leaf_index, point_t* point, int parent_branch_index))
{
  START_FUNCTION_TIMER();
  if (order == OCTREE_PRE)
    octree_visit_pre(tree->root, context, visit_branch, visit_leaf);
  else // if (order == OCTREE_POST)
    octree_visit_post(tree->root, context, visit_branch, visit_leaf);
  STOP_FUNCTION_TIMER();
}

