#include "core/octree.h"

#ifdef SQ
#undef SQ
#endif
#define SQ(x) ((x)*(x))
#define SQ_DIST(x, y) (SQ(x[0]-y[0]) + SQ(x[1]-y[1]) + SQ(x[2]-y[2]))

// A 3D octree.

typedef enum
{
  OCTREE_BRANCH_NODE,
  OCTREE_LEAF_NODE
} octree_node_type_t;

// Basic octree node -- can be either a branch node or a leaf node.
// We cast between this type and the actual node types below.
typedef struct 
{
  octree_node_type_t type;
} octree_node_t;

// Branch node.
typedef struct 
{
  // Present for compatibility with octree_node_t.
  octree_node_type_t type;
  
  // Children.
  octree_node_t* children[2][2][2];
} octree_branch_node_t;

// Leaf node.
typedef struct
{
  // Present for compatibility with octree_node_t.
  octree_node_type_t type;
  
  // Index of point.
  int index;
}

static octree_node_t* branch_new()
{
  octree_branch_node_t* branch = malloc(sizeof(octree_branch_node_t));
  branch->type = OCTREE_BRANCH_NODE;
  memset(branch->children, 0, 8*sizeof(octree_node_t*));
  return (octree_node_t*)branch;
}

static octree_node_t* leaf_new(int index)
{
  octree_leaf_node_t* leaf = malloc(sizeof(octree_leaf_node_t));
  leaf->type = OCTREE_LEAF_NODE;
  leaf->index = index;
  return (octree_node_t*)leaf;
}

static void node_delete(octree_node_t** node)
{
  if (*node == NULL) return;

  if ((*node)->type == OCTREE_BRANCH_NODE)
  {
    octree_branch_node_t* branch = *node;
    for (int i = 0; i < 2; ++i)
    {
      for (int j = 0; j < 2; ++j)
      {
        for (int k = 0; k < 2; ++k)
          node_delete(&branch->children[i][j][k]);
      }
    }
    free(branch);
  }
  else
  {
    octree_leaf_node_t* leaf = *node;
    free(leaf);
  }
  *node = NULL;
}

struct octree_t 
{
  octree_node_t* root;  // The root node.
  double size; // The size of the cubic block of space covered by the tree.
  int num_points; // Total number of points in the tree.
  int num_bins; // The number of bins on the side of the space.
};

octree_t* octree_new(point_t* center, double size, int num_bins)
{
  octree_t* tree = malloc(sizeof(octree_t));
  tree->center = *center;
  tree->size = size;
  tree->root = NULL;
  tree->rect = NULL;
  tree->size = 0;
  return tree;
}

void octree_free(octree_t* tree)
{
  octree_clear(tree);
  free(tree);
}

static void node_insert(octree_node_t** node_ptr, double* pos, int index)
{
  if (*node_ptr == NULL)
  {
    *node_ptr = node_new(pos, index);
    return;
  }
  octree_node_t* node = *node_ptr;
  for (int d = 0; d < 3; ++d)
  {
    if (pos[node->dir] < node->pos[node->dir])
    node_insert(&(*node_ptr)->left, pos, index, new_dir);
  else
    node_insert(&(*node_ptr)->right, pos, index, new_dir);
}

void octree_insert(octree_t* tree, point_t* point, int index)
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

int octree_size(octree_t* tree)
{
  return tree->size;
}

static void node_clear(octree_node_t* node)
{
  if (node == NULL) return;

  node_clear(node->left);
  node_clear(node->right);
  node_free(node);
}

void octree_clear(octree_t* tree)
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

static void find_nearest(octree_node_t* node, double* pos, octree_node_t** result, double* r2, octree_rect_t* rect)
{
  // Go left or right?
  int dir = node->dir;
  double which = pos[dir] - node->pos[dir];
  octree_node_t *near_subtree, *far_subtree;
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

int octree_nearest(octree_t* tree, point_t* point)
{
  if (tree->root == NULL)
    return -1;

  // Start with the root.
  octree_node_t* node = tree->root;
  double pos[3];
  pos[0] = point->x, pos[1] = point->y, pos[2] = point->z;
  double r2 = SQ_DIST(pos, node->pos);
  octree_rect_t rect;
  rect.min[0] = tree->rect->min[0]; rect.max[0] = tree->rect->max[0];
  rect.min[1] = tree->rect->min[1]; rect.max[1] = tree->rect->max[1];
  rect.min[2] = tree->rect->min[2]; rect.max[2] = tree->rect->max[2];
  
  // Search recursively for the closest node.
  find_nearest(tree->root, pos, &node, &r2, &rect);
  return node->index;
}

static void find_within_radius(octree_node_t* node, 
                               double* pos, 
                               double radius,
                               octree_rect_t* rect,
                               int_slist_t* results)
{
  // Go left or right?
  int dir = node->dir;
  double which = pos[dir] - node->pos[dir];
  octree_node_t *near_subtree, *far_subtree;
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

int_slist_t* octree_within_radius(octree_t* tree, 
                                     point_t* point, 
                                     double radius)
{
  if (tree->root == NULL)
    return NULL;

  // Start with the root.
  octree_node_t* node = tree->root;
  double pos[3];
  pos[0] = point->x, pos[1] = point->y, pos[2] = point->z;
  double r2 = SQ_DIST(pos, node->pos);
  octree_rect_t rect;
  int_slist_t* results = int_slist_new();
  rect.min[0] = tree->rect->min[0]; rect.max[0] = tree->rect->max[0];
  rect.min[1] = tree->rect->min[1]; rect.max[1] = tree->rect->max[1];
  rect.min[2] = tree->rect->min[2]; rect.max[2] = tree->rect->max[2];
  
  // Search recursively for the closest node.
  find_within_radius(tree->root, pos, radius, &rect, results);
  return results;
}

octree_pos_t octree_start(octree_t* tree)
{
  octree_pos_t pos = {.node = tree->root };
  return pos; 
}

bool octree_next(octree_t* tree, octree_pos_t* pos, int* index, double* coords)
{
  ASSERT(pos != NULL);
  ASSERT(index != NULL);
  ASSERT(coords != NULL);

  octree_node_t* node = pos->node;
  if (node == NULL)
    return false;

  if (node->left != NULL)
  {
    pos->node = node->left;
    return octree_next(tree, pos, index, coords);
  }
  *index = node->index;
  coords[0] = node->pos[0];
  coords[1] = node->pos[1];
  coords[2] = node->pos[2];
  if (node->right != NULL)
  {
    pos->node = node->right;
    return octree_next(tree, pos, index, coords);
  }
  return true;
}


