// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/kd_tree.h"

// Convenience functions.
static inline real_t square(real_t x)
{
  return x * x;
}

static inline real_t square_dist(real_t* x, real_t* y)
{
  return square(x[0]-y[0]) + square(x[1]-y[1]) + square(x[2]-y[2]);
}

// A kd-tree with k = 3.

// This, the basic structure in the tree, is basically a 3d-tree node.
typedef struct kd_tree_node_t kd_tree_node_t;
struct kd_tree_node_t
{
  real_t pos[3];  // Coordinates of point.
  int dir;        // Direction of indexing (0, 1, or 2)
  int index;      // Identifying index of point.

  kd_tree_node_t *left, *right; // Children.
};

static kd_tree_node_t* node_new(real_t* pos, int index, int dir)
{
  ASSERT(dir >= 0);
  ASSERT(dir < 3);
  kd_tree_node_t* node = polymec_malloc(sizeof(kd_tree_node_t));
  node->pos[0] = pos[0]; node->pos[1] = pos[1]; node->pos[2] = pos[2];
  node->index = index;
  node->dir = dir;
  node->left = node->right = NULL;
  return node;
}

static void node_free(kd_tree_node_t* node)
{
  polymec_free(node);
}

// This defines a bounding 3-rectangle for points in the kd_tree.
typedef struct
{
  real_t min[3], max[3];
} kd_tree_rect_t;

static kd_tree_rect_t* rect_new(real_t* min, real_t* max)
{
  kd_tree_rect_t* rect = polymec_malloc(sizeof(kd_tree_rect_t));
  rect->min[0] = min[0]; rect->min[1] = min[1]; rect->min[2] = min[2];
  rect->max[0] = max[0]; rect->max[1] = max[1]; rect->max[2] = max[2];
  return rect;
}

static void rect_free(kd_tree_rect_t* rect)
{
  polymec_free(rect);
}

// Returns the square of the distance from the rectange to the vector
// with components given by pos.
static real_t rect_square_dist(kd_tree_rect_t* rect, real_t* pos)
{
  real_t r2 = 0.0;
  for (int i = 0; i < 3; ++i)
  {
    if (pos[i] < rect->min[i])
      r2 += square(rect->min[i] - pos[i]);
    else if (pos[i] > rect->max[i])
      r2 += square(rect->max[i] - pos[i]);
  }
  return r2;
}

// Returns true if any point within the given rectangle satisfies the given
// predicate with the coordinates pos, false if not.
static bool rect_satisfies_predicate(kd_tree_rect_t* rect,
                                     real_t* pos,
                                     bool (*predicate)(void* context, point_t* x, point_t* y),
                                     void* context)
{
  // Find the coordinates of the closest point in rect to pos and
  // store them in xi.
  real_t xi[3] = {pos[0], pos[1], pos[2]};
  for (int i = 0; i < 3; ++i)
  {
    if (pos[i] < rect->min[i])
      xi[i] = rect->min[i];
    else if (pos[i] > rect->max[i])
      xi[i] = rect->max[i];
  }

  // Does this closest point satisfy the predicate?
  point_t *x = (point_t*)pos, *y = (point_t*)xi;
  return predicate(context, x, y);
}

static void rect_extend(kd_tree_rect_t* rect, real_t* pos)
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
  size_t size;          // Number of points.

  // Point cache for traversal.
  real_array_t* point_cache;
};

static void node_insert(kd_tree_node_t** node_ptr, real_t* pos, int index, int dir)
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
  real_t pos[3];
  pos[0] = point->x; pos[1] = point->y; pos[2] = point->z;
  node_insert(&tree->root, pos, index, 0);
  ++(tree->size);
  if (tree->rect == NULL)
    tree->rect = rect_new(pos, pos);
  else
    rect_extend(tree->rect, pos);

  // Update the point cache.
  real_array_append(tree->point_cache, point->x);
  real_array_append(tree->point_cache, point->y);
  real_array_append(tree->point_cache, point->z);
}

kd_tree_t* kd_tree_new(point_t* points, size_t num_points)
{
  kd_tree_t* tree = polymec_malloc(sizeof(kd_tree_t));
  tree->root = NULL;
  tree->rect = NULL;
  tree->size = 0;
  tree->point_cache = real_array_new();

  for (int i = 0; i < num_points; ++i)
    kd_tree_insert(tree, &points[i], i);
  return tree;
}

size_t kd_tree_size(kd_tree_t* tree)
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

  // Clear the point cache.
  real_array_clear(tree->point_cache);
}

void kd_tree_free(kd_tree_t* tree)
{
  kd_tree_clear(tree);
  real_array_free(tree->point_cache);
  polymec_free(tree);
}

static void find_nearest(kd_tree_node_t* node, real_t* pos, kd_tree_node_t** result, real_t* r2, kd_tree_rect_t* rect)
{
  // Go left or right?
  int dir = node->dir;
  real_t which = pos[dir] - node->pos[dir];
  kd_tree_node_t *near_subtree, *far_subtree;
  real_t *near_coord, *far_coord;
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
    real_t coord = *near_coord;
    *near_coord = node->pos[dir];
    find_nearest(near_subtree, pos, result, r2, rect);
    *near_coord = coord;
  }

  // Distance from point to current node.
  real_t my_r2 = square_dist(node->pos, pos);
  if (my_r2 < *r2)
  {
    *result = node;
    *r2 = my_r2;
  }

  if (far_subtree != NULL)
  {
    // Bisect and recurse (if needed).
    real_t coord = *far_coord;
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
  real_t pos[3];
  pos[0] = point->x; pos[1] = point->y; pos[2] = point->z;
  real_t r2 = square_dist(pos, node->pos);
  kd_tree_rect_t rect;
  rect.min[0] = tree->rect->min[0]; rect.max[0] = tree->rect->max[0];
  rect.min[1] = tree->rect->min[1]; rect.max[1] = tree->rect->max[1];
  rect.min[2] = tree->rect->min[2]; rect.max[2] = tree->rect->max[2];

  // Search recursively for the closest node.
  find_nearest(tree->root, pos, &node, &r2, &rect);
  return node->index;
}

static void find_nearest_n(kd_tree_node_t* node,
                           real_t* pos,
                           int n,
                           int* neighbors,
                           real_t* square_distances,
                           kd_tree_rect_t* rect)
{
  // Go left or right?
  int dir = node->dir;
  real_t which = pos[dir] - node->pos[dir];
  kd_tree_node_t *near_subtree, *far_subtree;
  real_t *near_coord, *far_coord;
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
    real_t coord = *near_coord;
    *near_coord = node->pos[dir];
    find_nearest_n(near_subtree, pos, n, neighbors, square_distances, rect);
    *near_coord = coord;
  }

  // Get the square distance from point to the current node, and insert
  // the current node into the list of n nearest neighbors if it's closer
  // than any of the current nearest neighbors.
  real_t my_r2 = square_dist(node->pos, pos);
  if (my_r2 < square_distances[n-1])
  {
    int i = n-1;
    while ((i > 0) && (my_r2 < square_distances[i])) --i;
    if (my_r2 > square_distances[i])
      ++i; // Back up one step to where we'll insert the new neighbor.
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
    real_t coord = *far_coord;
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
    neighbors[i] = -1;

  if (tree->root == NULL)
    return;

  real_t square_distances[n];
  for (int i = 0; i < n; ++i)
    square_distances[i] = REAL_MAX;

  // Start with the root.
  real_t pos[3];
  pos[0] = point->x; pos[1] = point->y; pos[2] = point->z;
  kd_tree_rect_t rect;
  rect.min[0] = tree->rect->min[0]; rect.max[0] = tree->rect->max[0];
  rect.min[1] = tree->rect->min[1]; rect.max[1] = tree->rect->max[1];
  rect.min[2] = tree->rect->min[2]; rect.max[2] = tree->rect->max[2];

  // Search recursively for the closest nodes.
  find_nearest_n(tree->root, pos, n, neighbors, square_distances, &rect);
  ASSERT((neighbors[n-1] >= 0) || ((tree->size < n) && (neighbors[n-1] == -1)));
}

static void find_within_radius(kd_tree_node_t* node,
                               real_t* pos,
                               real_t radius,
                               kd_tree_rect_t* rect,
                               int_array_t* results)
{
  // Go left or right?
  int dir = node->dir;
  real_t which = pos[dir] - node->pos[dir];
  kd_tree_node_t *near_subtree, *far_subtree;
  real_t *near_coord, *far_coord;
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
    real_t coord = *near_coord;
    *near_coord = node->pos[dir];
    find_within_radius(near_subtree, pos, radius, rect, results);
    *near_coord = coord;
  }

  // Distance from point to current node.
  real_t my_r2 = square_dist(node->pos, pos);
  if (my_r2 < radius*radius)
    int_array_append(results, node->index);

  if (far_subtree != NULL)
  {
    // Bisect and recurse (if needed).
    real_t coord = *far_coord;
    *far_coord = node->pos[dir];
    if (rect_square_dist(rect, pos) < radius*radius)
      find_within_radius(far_subtree, pos, radius, rect, results);
    *far_coord = coord;
  }
}

int_array_t* kd_tree_within_radius(kd_tree_t* tree,
                                   point_t* point,
                                   real_t radius)
{
  int_array_t* results = int_array_new();
  if (tree->root == NULL)
    return results;

  real_t pos[3];
  pos[0] = point->x; pos[1] = point->y; pos[2] = point->z;
  kd_tree_rect_t rect;
  rect.min[0] = tree->rect->min[0]; rect.max[0] = tree->rect->max[0];
  rect.min[1] = tree->rect->min[1]; rect.max[1] = tree->rect->max[1];
  rect.min[2] = tree->rect->min[2]; rect.max[2] = tree->rect->max[2];

  // Search recursively for the closest node.
  find_within_radius(tree->root, pos, radius, &rect, results);
  return results;
}

static void find_for_predicate(kd_tree_node_t* node,
                               real_t* pos,
                               bool (*pred)(void* context, point_t* x, point_t* y),
                               void* context,
                               kd_tree_rect_t* rect,
                               int_array_t* results)
{
  // Go left or right?
  int dir = node->dir;
  real_t which = pos[dir] - node->pos[dir];
  kd_tree_node_t *near_subtree, *far_subtree;
  real_t *near_coord, *far_coord;
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
    real_t coord = *near_coord;
    *near_coord = node->pos[dir];
    find_for_predicate(near_subtree, pos, pred, context, rect, results);
    *near_coord = coord;
  }

  // Does the current node satisfy the predicate?
  point_t x = {pos[0], pos[1], pos[2]};
  point_t y = {node->pos[0], node->pos[1], node->pos[2]};
  if (pred(context, &x, &y))
    int_array_append(results, node->index);

  if (far_subtree != NULL)
  {
    // Bisect and recurse (if needed).
    real_t coord = *far_coord;
    *far_coord = node->pos[dir];
    if (rect_satisfies_predicate(rect, pos, pred, context))
      find_for_predicate(far_subtree, pos, pred, context, rect, results);
    *far_coord = coord;
  }
}

int_array_t* kd_tree_for_predicate(kd_tree_t* tree,
                                   point_t* point,
                                   bool (*pred)(void* context, point_t* point, point_t* x),
                                   void* context)
{
  int_array_t* results = int_array_new();
  if (tree->root == NULL)
    return results;

  real_t pos[3];
  pos[0] = point->x; pos[1] = point->y; pos[2] = point->z;
  kd_tree_rect_t rect;
  rect.min[0] = tree->rect->min[0]; rect.max[0] = tree->rect->max[0];
  rect.min[1] = tree->rect->min[1]; rect.max[1] = tree->rect->max[1];
  rect.min[2] = tree->rect->min[2]; rect.max[2] = tree->rect->max[2];

  // Search recursively for the closest node.
  find_for_predicate(tree->root, pos, pred, context, &rect, results);
  return results;
}

bool kd_tree_next(kd_tree_t* tree, int* pos, int* index, point_t* coords)
{
  ASSERT(*pos >= 0);
  ASSERT(index != NULL);
  ASSERT(coords != NULL);

  int offset = 3 * (*pos);
  if (offset >= tree->point_cache->size)
    return false;

  *index = *pos;
  coords->x = tree->point_cache->data[offset];
  coords->y = tree->point_cache->data[offset+1];
  coords->z = tree->point_cache->data[offset+2];
  ++(*pos);

  return true;
}

exchanger_t* kd_tree_find_ghost_points(kd_tree_t* tree, MPI_Comm comm, real_t R_max)
{
  ASSERT(R_max >= 0.0);
  exchanger_t* ex = exchanger_new(comm);
#if POLYMEC_HAVE_MPI

  // Create a bounding box containing all the local points in the tree.
  bbox_t bbox = {.x1 = REAL_MAX, .x2 = -REAL_MAX,
                 .y1 = REAL_MAX, .y2 = -REAL_MAX,
                 .z1 = REAL_MAX, .z2 = -REAL_MAX};
  {
    int pos = 0, index;
    point_t x;
    while (kd_tree_next(tree, &pos, &index, &x))
    {
      if (!bbox_contains(&bbox, &x))
        bbox_grow(&bbox, &x);
    }
  }

  // Grow our bounding box outward by the maximum radius.
  bbox.x1 -= R_max;
  bbox.x2 += R_max;
  bbox.y1 -= R_max;
  bbox.y2 += R_max;
  bbox.z1 -= R_max;
  bbox.z2 += R_max;

  // Find out which processes we interact with by intersecting bounding boxes.
  int num_neighbor_procs;
  int* neighbor_procs = bbox_intersecting_processes(&bbox, comm, &num_neighbor_procs);

  int rank;
  MPI_Comm_rank(comm, &rank);

  // Get bounding boxes from all of these processes.
  bbox_t bboxes[num_neighbor_procs];
  {
    real_t data[num_neighbor_procs][6];
    MPI_Request requests[2*num_neighbor_procs];
    for (int p = 0; p < num_neighbor_procs; ++p)
    {
      int proc = neighbor_procs[p];
      int err = MPI_Irecv(data[p], 6, MPI_REAL_T, proc, 0, comm, &requests[p]);
      if (err != MPI_SUCCESS)
        polymec_error("%d: Could not receive bounding box data from %d", rank, proc);
    }

    // Now send asynchronously.
    real_t my_data[6] = {bbox.x1, bbox.x2, bbox.y1, bbox.y2, bbox.z1, bbox.z2};
    for (int p = 0; p < num_neighbor_procs; ++p)
    {
      int proc = neighbor_procs[p];
      int err = MPI_Isend(my_data, 6, MPI_REAL_T, proc, 0, comm, &requests[num_neighbor_procs + p]);
      if (err != MPI_SUCCESS)
        polymec_error("%d: Could not send point numbers to %d", rank, proc);
    }

    // Wait for messages to fly.
    MPI_Status statuses[2*num_neighbor_procs];
    MPI_Waitall(2*num_neighbor_procs, requests, statuses);

    // Copy the data into place.
    for (int p = 0; p < num_neighbor_procs; ++p)
    {
      bboxes[p].x1 = data[p][0];
      bboxes[p].x2 = data[p][1];
      bboxes[p].y1 = data[p][2];
      bboxes[p].y2 = data[p][3];
      bboxes[p].z1 = data[p][4];
      bboxes[p].z2 = data[p][5];
    }
  }

  // Now send each neighboring process the points that fall within its
  // bounding box.
  {
    real_array_t* points_sent[num_neighbor_procs];
    MPI_Request requests[2*num_neighbor_procs];
    MPI_Status statuses[2*num_neighbor_procs];

    // Post receives for the number of points we expect.
    int num_points_received[num_neighbor_procs];
    for (int p = 0; p < num_neighbor_procs; ++p)
    {
      int proc = neighbor_procs[p];
      int err = MPI_Irecv(&num_points_received[p], 1, MPI_INT, proc, 0, comm, &requests[p]);
      if (err != MPI_SUCCESS)
        polymec_error("%d: Could not receive point numbers from %d", rank, proc);
    }

    for (int p = 0; p < num_neighbor_procs; ++p)
    {
      int_array_t* sent_point_indices = int_array_new();
      points_sent[p] = real_array_new();

      // Make a list of our points that fall into process p's
      // bounding box, and populate the send info in our exchanger.
      int pos = 0, index;
      point_t x;
      while (kd_tree_next(tree, &pos, &index, &x))
      {
        if (bbox_contains(&bboxes[p], &x))
        {
          real_array_append(points_sent[p], x.x);
          real_array_append(points_sent[p], x.y);
          real_array_append(points_sent[p], x.z);
          int_array_append(sent_point_indices, index);
        }
      }

      // Set the send information in the exchanger.
      int proc = neighbor_procs[p];

      if (sent_point_indices->size > 0)
        exchanger_set_send(ex, proc, sent_point_indices->data, (int)sent_point_indices->size, false);
      int_array_release_data_and_free(sent_point_indices);

      // Now send the number of points to process p.
      int num_points_sent = (int)(points_sent[p]->size/3);
      int err = MPI_Isend(&num_points_sent, 1, MPI_INT,
                          proc, 0, comm, &requests[num_neighbor_procs + p]);
      if (err != MPI_SUCCESS)
        polymec_error("%d: Could not send point numbers to %d", rank, proc);
    }

    // Finish the sends.
    MPI_Waitall(2*num_neighbor_procs, requests, statuses);

    // Post receives for the actual point coordinates.
    real_array_t* points_received[num_neighbor_procs];
    for (int p = 0; p < num_neighbor_procs; ++p)
    {
      points_received[p] = real_array_new();
      real_array_resize(points_received[p], 3*num_points_received[p]);

      int proc = neighbor_procs[p];
      int err = MPI_Irecv(points_received[p]->data, (int)points_received[p]->size,
                          MPI_REAL_T, proc, 0, comm, &requests[p]);
      if (err != MPI_SUCCESS)
        polymec_error("%d: Could not receive point data from %d", rank, proc);
    }

    // Now send the point coordinates.
    for (int p = 0; p < num_neighbor_procs; ++p)
    {
      // Now send the number of points to process p.
      int proc = neighbor_procs[p];
      int err = MPI_Isend(points_sent[p]->data, (int)points_sent[p]->size, MPI_REAL_T,
                          proc, 0, comm, &requests[num_neighbor_procs + p]);
      if (err != MPI_SUCCESS)
        polymec_error("%d: Could not send point numbers to %d", rank, proc);
    }

    // Finish the sends.
    MPI_Waitall(2*num_neighbor_procs, requests, statuses);

    // Add the points to our tree and populate the receive entries in
    // our exchanger.
    for (int p = 0; p < num_neighbor_procs; ++p)
    {
      int_array_t* received_point_indices = int_array_new();
      real_array_t* point_data = points_received[p];
      for (int i = 0; i < point_data->size/3; ++i)
      {
        point_t x = {.x = point_data->data[3*i],
                     .y = point_data->data[3*i+1],
                     .z = point_data->data[3*i+2]};
        int ghost_index = (int)tree->size;
        kd_tree_insert(tree, &x, ghost_index);
        int_array_append(received_point_indices, ghost_index);
      }

      // Set the send information in the exchanger.
      int proc = neighbor_procs[p];
      if (received_point_indices->size > 0)
        exchanger_set_receive(ex, proc, received_point_indices->data, (int)received_point_indices->size, false);
      int_array_release_data_and_free(received_point_indices);
    }

    // Clean up.
    for (int p = 0; p < num_neighbor_procs; ++p)
    {
      real_array_free(points_sent[p]);
      real_array_free(points_received[p]);
    }
  }

  // Clean up.
  polymec_free(neighbor_procs);

  return ex;
#else
  return ex;
#endif
}

