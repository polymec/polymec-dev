// Copyright 2012-2013 Jeffrey Johnson.
// 
// This file is part of Polymec, and is licensed under the Apache License, 
// Version 2.0 (the "License"); you may not use this file except in 
// compliance with the License. You may may find the text of the license in 
// the LICENSE file at the top-level source directory, or obtain a copy of 
// it at
// 
// http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "core/kd_tree.h"

#ifdef SQ
#undef SQ
#endif
#define SQ(x) ((x)*(x))
#define SQ_DIST(x, y) (SQ(x[0]-y[0]) + SQ(x[1]-y[1]) + SQ(x[2]-y[2]))

// A kd-tree with k = 3.

// This, the basic structure in the tree, is basically a 3d-tree node.
typedef struct kd_tree_node_t kd_tree_node_t;
struct kd_tree_node_t
{
  double pos[3];  // Coordinates of point.
  int dir;        // Direction of indexing (0, 1, or 2)
  int index;      // Identifying index of point.

  kd_tree_node_t *left, *right; // Children.
};

static kd_tree_node_t* node_new(double* pos, int index, int dir)
{
  ASSERT(dir >= 0);
  ASSERT(dir < 3);
  kd_tree_node_t* node = malloc(sizeof(kd_tree_node_t));
  node->pos[0] = pos[0]; node->pos[1] = pos[1]; node->pos[2] = pos[2];
  node->index = index;
  node->dir = dir;
  node->left = node->right = NULL;
  return node;
}

static void node_free(kd_tree_node_t* node)
{
  free(node);
}

// This defines a bounding 3-rectangle for points in the kd_tree.
typedef struct 
{
  double min[3], max[3];
} kd_tree_rect_t;

static kd_tree_rect_t* rect_new(double* min, double* max)
{
  kd_tree_rect_t* rect = malloc(sizeof(kd_tree_rect_t));
  rect->min[0] = min[0]; rect->min[1] = min[1]; rect->min[2] = min[2];
  rect->max[0] = max[0]; rect->max[1] = max[1]; rect->max[2] = max[2];
  return rect;
}

static void rect_free(kd_tree_rect_t* rect)
{
  free(rect);
}

static double rect_square_dist(kd_tree_rect_t* rect, double* pos)
{
  double r2 = 0.0;
  for (int i = 0; i < 3; ++i)
  {
    if (pos[i] < rect->min[i])
      r2 += SQ(rect->min[i] - pos[i]);
    else if (pos[i] > rect->max[i])
      r2 += SQ(rect->max[i] - pos[i]);
  }
  return r2;
}

static void rect_extend(kd_tree_rect_t* rect, double* pos)
{
  for (int i = 0; i < 3; ++i)
  {
    if (pos[i] < rect->min[i])
      rect->min[i] = pos[i];
    if (pos[i] > rect->max[i])
      rect->max[i] = pos[i];
  }
}

struct kd_tree_t 
{
  kd_tree_node_t* root; // Root node of 3d-tree.
  kd_tree_rect_t* rect; // Containing rectangle.
  int size;               // Number of points.
};

static void node_insert(kd_tree_node_t** node_ptr, double* pos, int index, int dir)
{
  if (*node_ptr == NULL)
  {
    *node_ptr = node_new(pos, index, dir);
    return;
  }
  kd_tree_node_t* node = *node_ptr;
  int new_dir = (node->dir + 1) % 3;
  if (pos[node->dir] < node->pos[node->dir])
    node_insert(&(*node_ptr)->left, pos, index, new_dir);
  else
    node_insert(&(*node_ptr)->right, pos, index, new_dir);
}

static void kd_tree_insert(kd_tree_t* tree, point_t* point, int index)
{
  double pos[3];
  pos[0] = point->x, pos[1] = point->y, pos[2] = point->z;
  node_insert(&tree->root, pos, index, 0);
  ++(tree->size);
  if (tree->rect == NULL)
    tree->rect = rect_new(pos, pos);
  else
    rect_extend(tree->rect, pos);
}

kd_tree_t* kd_tree_new(point_t* points, int num_points)
{
  kd_tree_t* tree = malloc(sizeof(kd_tree_t));
  tree->root = NULL;
  tree->rect = NULL;
  tree->size = 0;

  for (int i = 0; i < num_points; ++i)
    kd_tree_insert(tree, &points[i], i);
  return tree;
}


int kd_tree_size(kd_tree_t* tree)
{
  return tree->size;
}

static void node_clear(kd_tree_node_t* node)
{
  if (node == NULL) return;

  node_clear(node->left);
  node_clear(node->right);
  node_free(node);
}

static void kd_tree_clear(kd_tree_t* tree)
{
  node_clear(tree->root);
  tree->root = NULL;
  if (tree->rect != NULL)
  {
    rect_free(tree->rect);
    tree->rect = NULL;
  }
  tree->size = 0;
}

void kd_tree_free(kd_tree_t* tree)
{
  kd_tree_clear(tree);
  free(tree);
}

static void find_nearest(kd_tree_node_t* node, double* pos, kd_tree_node_t** result, double* r2, kd_tree_rect_t* rect)
{
  // Go left or right?
  int dir = node->dir;
  double which = pos[dir] - node->pos[dir];
  kd_tree_node_t *near_subtree, *far_subtree;
  double *near_coord, *far_coord;
  if (which <= 0.0)
  {
    near_subtree = node->left;
    far_subtree = node->right;
    near_coord = rect->max + dir;
    far_coord = rect->min + dir;
  }
  else
  {
    near_subtree = node->right;
    far_subtree = node->left;
    near_coord = rect->min + dir;
    far_coord = rect->max + dir;
  }

  if (near_subtree != NULL)
  {
    // Bisect and recurse.
    double coord = *near_coord;
    *near_coord = node->pos[dir];
    find_nearest(near_subtree, pos, result, r2, rect);
    *near_coord = coord;
  }

  // Distance from point to current node.
  double my_r2 = SQ_DIST(node->pos, pos);
  if (my_r2 < *r2)
  {
    *result = node;
    *r2 = my_r2;
  }

  if (far_subtree != NULL)
  {
    // Bisect and recurse (if needed).
    double coord = *far_coord;
    *far_coord = node->pos[dir];
    if (rect_square_dist(rect, pos) < *r2)
      find_nearest(far_subtree, pos, result, r2, rect);
    *far_coord = coord;
  }
}

int kd_tree_nearest(kd_tree_t* tree, point_t* point)
{
  if (tree->root == NULL)
    return -1;

  // Start with the root.
  kd_tree_node_t* node = tree->root;
  double pos[3];
  pos[0] = point->x, pos[1] = point->y, pos[2] = point->z;
  double r2 = SQ_DIST(pos, node->pos);
  kd_tree_rect_t rect;
  rect.min[0] = tree->rect->min[0]; rect.max[0] = tree->rect->max[0];
  rect.min[1] = tree->rect->min[1]; rect.max[1] = tree->rect->max[1];
  rect.min[2] = tree->rect->min[2]; rect.max[2] = tree->rect->max[2];
  
  // Search recursively for the closest node.
  find_nearest(tree->root, pos, &node, &r2, &rect);
  return node->index;
}

static void find_nearest_n(kd_tree_node_t* node, 
                           double* pos, 
                           int n, 
                           int* neighbors, 
                           double* square_distances, 
                           kd_tree_rect_t* rect)
{
  // Go left or right?
  int dir = node->dir;
  double which = pos[dir] - node->pos[dir];
  kd_tree_node_t *near_subtree, *far_subtree;
  double *near_coord, *far_coord;
  if (which <= 0.0)
  {
    near_subtree = node->left;
    far_subtree = node->right;
    near_coord = rect->max + dir;
    far_coord = rect->min + dir;
  }
  else
  {
    near_subtree = node->right;
    far_subtree = node->left;
    near_coord = rect->min + dir;
    far_coord = rect->max + dir;
  }

  if (near_subtree != NULL)
  {
    // Bisect and recurse.
    double coord = *near_coord;
    *near_coord = node->pos[dir];
    find_nearest_n(near_subtree, pos, n, neighbors, square_distances, rect);
    *near_coord = coord;
  }

  // Get the square distance from point to the current node, and insert 
  // the current node into the list of n nearest neighbors if it's closer 
  // than any of the current nearest neighbors.
  double my_r2 = SQ_DIST(node->pos, pos);
  if (my_r2 < square_distances[n-1])
  {
    int i = n-1;
    while ((i > 0) && (my_r2 < square_distances[i])) --i;
    if (my_r2 > square_distances[i])
      ++i; // Back up one step to where we'll insert the new neighbor.
printf("%g vs %g: %d slotted at %d\n", my_r2, square_distances[i], node->index, i);
    for (int j = n-1; j > i; --j)
    {
      neighbors[j] = neighbors[j-1];
      square_distances[j] = square_distances[j-1];
    }
    neighbors[i] = node->index;
    square_distances[i] = my_r2;
  }

  if (far_subtree != NULL)
  {
    // Bisect and recurse (if needed).
    double coord = *far_coord;
    *far_coord = node->pos[dir];
    if (rect_square_dist(rect, pos) < square_distances[n-1])
      find_nearest_n(far_subtree, pos, n, neighbors, square_distances, rect);
    *far_coord = coord;
  }
}

void kd_tree_nearest_n(kd_tree_t* tree, point_t* point, int n, int* neighbors)
{
  ASSERT(n > 0);
  ASSERT(neighbors != NULL);
  for (int i = 0; i < n; ++i)
    neighbors[n] = -1;

  if (tree->root == NULL)
    return;

  double square_distances[n];
  for (int i = 0; i < n; ++i)
    square_distances[i] = FLT_MAX;

  // Start with the root.
  double pos[3];
  pos[0] = point->x, pos[1] = point->y, pos[2] = point->z;
  kd_tree_rect_t rect;
  rect.min[0] = tree->rect->min[0]; rect.max[0] = tree->rect->max[0];
  rect.min[1] = tree->rect->min[1]; rect.max[1] = tree->rect->max[1];
  rect.min[2] = tree->rect->min[2]; rect.max[2] = tree->rect->max[2];
  
  // Search recursively for the closest nodes.
  printf("tree is %p\n", tree);
  find_nearest_n(tree->root, pos, n, neighbors, square_distances, &rect);
  printf("tree is now %p\n", tree);
  ASSERT((neighbors[n-1] >= 0) || ((tree->size < n) && (neighbors[n-1] == -1)));
}

static void find_within_radius(kd_tree_node_t* node, 
                               double* pos, 
                               double radius,
                               kd_tree_rect_t* rect,
                               int_slist_t* results)
{
  // Go left or right?
  int dir = node->dir;
  double which = pos[dir] - node->pos[dir];
  kd_tree_node_t *near_subtree, *far_subtree;
  double *near_coord, *far_coord;
  if (which <= 0.0)
  {
    near_subtree = node->left;
    far_subtree = node->right;
    near_coord = rect->max + dir;
    far_coord = rect->min + dir;
  }
  else
  {
    near_subtree = node->right;
    far_subtree = node->left;
    near_coord = rect->min + dir;
    far_coord = rect->max + dir;
  }

  if (near_subtree != NULL)
  {
    // Bisect and recurse.
    double coord = *near_coord;
    *near_coord = node->pos[dir];
    find_within_radius(near_subtree, pos, radius, rect, results);
    *near_coord = coord;
  }

  // Distance from point to current node.
  double my_r2 = SQ_DIST(node->pos, pos);
  if (my_r2 < radius*radius)
    int_slist_append(results, node->index);

  if (far_subtree != NULL)
  {
    // Bisect and recurse (if needed).
    double coord = *far_coord;
    *far_coord = node->pos[dir];
    if (rect_square_dist(rect, pos) < radius*radius)
      find_within_radius(far_subtree, pos, radius, rect, results);
    *far_coord = coord;
  }
}

int_slist_t* kd_tree_within_radius(kd_tree_t* tree, 
                                     point_t* point, 
                                     double radius)
{
  if (tree->root == NULL)
    return NULL;

  // Start with the root.
  kd_tree_node_t* node = tree->root;
  double pos[3];
  pos[0] = point->x, pos[1] = point->y, pos[2] = point->z;
  double r2 = SQ_DIST(pos, node->pos);
  kd_tree_rect_t rect;
  int_slist_t* results = int_slist_new();
  rect.min[0] = tree->rect->min[0]; rect.max[0] = tree->rect->max[0];
  rect.min[1] = tree->rect->min[1]; rect.max[1] = tree->rect->max[1];
  rect.min[2] = tree->rect->min[2]; rect.max[2] = tree->rect->max[2];
  
  // Search recursively for the closest node.
  find_within_radius(tree->root, pos, radius, &rect, results);
  return results;
}

kd_tree_pos_t kd_tree_start(kd_tree_t* tree)
{
  kd_tree_pos_t pos = {.node = tree->root };
  return pos; 
}

bool kd_tree_next(kd_tree_t* tree, kd_tree_pos_t* pos, int* index, double* coords)
{
  ASSERT(pos != NULL);
  ASSERT(index != NULL);
  ASSERT(coords != NULL);

  kd_tree_node_t* node = pos->node;
  if (node == NULL)
    return false;

  if (node->left != NULL)
  {
    pos->node = node->left;
    return kd_tree_next(tree, pos, index, coords);
  }
  *index = node->index;
  coords[0] = node->pos[0];
  coords[1] = node->pos[1];
  coords[2] = node->pos[2];
  if (node->right != NULL)
  {
    pos->node = node->right;
    return kd_tree_next(tree, pos, index, coords);
  }
  return true;
}


