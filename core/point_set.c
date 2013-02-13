#include "core/point_set.h"

#ifdef __cplusplus
extern "C" {
#endif

#ifdef SQ
#undef SQ
#endif
#define SQ(x) ((x)*(x))
#define SQ_DIST(x, y) (SQ(x[0]-y[0]) + SQ(x[1]-y[1]) + SQ(x[2]-y[2]))

// The point set is an embellished kd-tree with k = 3.

// This, the basic structure in the point set, is basically a 3d-tree node.
typedef struct point_set_node_t point_set_node_t;
struct point_set_node_t
{
  double pos[3];  // Coordinates of point.
  int dir;        // Direction of indexing (0, 1, or 2)
  int index;      // Identifying index of point.

  point_set_node_t *left, *right; // Children.
};

static point_set_node_t* node_new(double* pos, int index, int dir)
{
  ASSERT(dir >= 0);
  ASSERT(dir < 3);
  point_set_node_t* node = malloc(sizeof(point_set_node_t));
  node->pos[0] = pos[0]; node->pos[1] = pos[1]; node->pos[2] = pos[2];
  node->index = index;
  node->dir = dir;
  node->left = node->right = NULL;
  return node;
}

static void node_free(point_set_node_t* node)
{
  free(node);
}

// This defines a bounding 3-rectangle for points in the point set.
typedef struct 
{
  double min[3], max[3];
} point_set_rect_t;

static point_set_rect_t* rect_new(double* min, double* max)
{
  point_set_rect_t* rect = malloc(sizeof(point_set_rect_t));
  rect->min[0] = min[0]; rect->min[1] = min[1]; rect->min[2] = min[2];
  rect->max[0] = max[0]; rect->max[1] = max[1]; rect->max[2] = max[2];
  return rect;
}

static void rect_free(point_set_rect_t* rect)
{
  free(rect);
}

static double rect_square_dist(point_set_rect_t* rect, double* pos)
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

static void rect_extend(point_set_rect_t* rect, double* pos)
{
  for (int i = 0; i < 3; ++i)
  {
    if (pos[i] < rect->min[i])
      rect->min[i] = pos[i];
    if (pos[i] > rect->max[i])
      rect->max[i] = pos[i];
  }
}

struct point_set_t 
{
  point_set_node_t* root; // Root node of 3d-tree.
  point_set_rect_t* rect; // Containing rectangle.
  int size;               // Number of points.
};

point_set_t* point_set_new()
{
  point_set_t* pset = malloc(sizeof(point_set_t));
  pset->root = NULL;
  pset->rect = NULL;
  pset->size = 0;
  return pset;
}

void point_set_free(point_set_t* pset)
{
  point_set_clear(pset);
  free(pset);
}

static void node_insert(point_set_node_t** node_ptr, double* pos, int index, int dir)
{
  if (*node_ptr == NULL)
  {
    *node_ptr = node_new(pos, index, dir);
    return;
  }
  point_set_node_t* node = *node_ptr;
  int new_dir = (node->dir + 1) % 3;
  if (pos[node->dir] < node->pos[node->dir])
    node_insert(&(*node_ptr)->left, pos, index, new_dir);
  else
    node_insert(&(*node_ptr)->right, pos, index, new_dir);
}

void point_set_insert(point_set_t* pset, point_t* point, int index)
{
  double pos[3];
  pos[0] = point->x, pos[1] = point->y, pos[2] = point->z;
  node_insert(&pset->root, pos, index, 0);
  ++(pset->size);
  if (pset->rect == NULL)
    pset->rect = rect_new(pos, pos);
  else
    rect_extend(pset->rect, pos);
}

int point_set_size(point_set_t* pset)
{
  return pset->size;
}

static void node_clear(point_set_node_t* node)
{
  if (node == NULL) return;

  node_clear(node->left);
  node_clear(node->right);
  node_free(node);
}

void point_set_clear(point_set_t* pset)
{
  node_clear(pset->root);
  pset->root = NULL;
  if (pset->rect != NULL)
  {
    rect_free(pset->rect);
    pset->rect = NULL;
  }
  pset->size = 0;
}

static void find_nearest(point_set_node_t* node, double* pos, point_set_node_t** result, double* r2, point_set_rect_t* rect)
{
  // Go left or right?
  int dir = node->dir;
  double which = pos[dir] - node->pos[dir];
  point_set_node_t *near_subtree, *far_subtree;
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

int point_set_nearest(point_set_t* pset, point_t* point)
{
  if (pset->root == NULL)
    return -1;

  // Start with the root.
  point_set_node_t* node = pset->root;
  double pos[3];
  pos[0] = point->x, pos[1] = point->y, pos[2] = point->z;
  double r2 = SQ_DIST(pos, node->pos);
  point_set_rect_t rect;
  rect.min[0] = pset->rect->min[0]; rect.max[0] = pset->rect->max[0];
  rect.min[1] = pset->rect->min[1]; rect.max[1] = pset->rect->max[1];
  rect.min[2] = pset->rect->min[2]; rect.max[2] = pset->rect->max[2];
  
  // Search recursively for the closest node.
  find_nearest(pset->root, pos, &node, &r2, &rect);
  return node->index;
}

point_set_pos_t point_set_start(point_set_t* pset)
{
  point_set_pos_t pos = {.node = pset->root };
  return pos; 
}

bool point_set_next(point_set_t* pset, point_set_pos_t* pos, int* index, double* coords)
{
  ASSERT(pos != NULL);
  ASSERT(index != NULL);
  ASSERT(coords != NULL);

  point_set_node_t* node = pos->node;
  if (node == NULL)
    return false;

  if (node->left != NULL)
  {
    pos->node = node->left;
    return point_set_next(pset, pos, index, coords);
  }
  *index = node->index;
  coords[0] = node->pos[0];
  coords[1] = node->pos[1];
  coords[2] = node->pos[2];
  if (node->right != NULL)
  {
    pos->node = node->right;
    return point_set_next(pset, pos, index, coords);
  }
  return true;
}


#ifdef __cplusplus
}
#endif

