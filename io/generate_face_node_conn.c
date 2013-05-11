#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "core/point.h"
#include "core/slist.h"
#include "core/avl_tree.h"
#include "core/mesh.h"
#include "io/generate_face_node_conn.h"

// Traverses the given points of a polygonal facet along their convex
// hull, writing their indices to indices in order.
static void traverse_convex_hull(double* points, int num_points, int* indices, int* count)
{
  *count = 0;

  // Find the "lowest" point in the set.
  double ymin = FLT_MAX;
  int index0 = -1;
  for (int p = 0; p < num_points; ++p)
  {
    if (ymin > points[2*p+1])
    {
      ymin = points[2*p+1];
      index0 = p;
    }
  }

  // We start with this point and a horizontal angle.
  double theta_prev = 0.0;
  indices[(*count)++] = index0;

  // Now start gift wrapping.
  int i = index0;
  do 
  {
    double dtheta_min = 2.0*M_PI;
    int j_min = -1;
    for (int j = 0; j < num_points; ++j)
    {
      if (j != i)
      {
        double dx = points[2*j] - points[2*i],
               dy = points[2*j+1] - points[2*i+1];
        double theta = atan2(dy, dx);
        double dtheta = theta - theta_prev;
        if (dtheta < 0.0)
          dtheta += 2.0*M_PI;
        if (dtheta_min > dtheta)
        {
          dtheta_min = dtheta;
          j_min = j;
        }
      }
    }
    if (j_min != index0)
      indices[(*count)++] = j_min;
    theta_prev += dtheta_min;
    i = j_min;
  }
  while (i != index0);

  // The convex hull should be a polygon unless the input points 
  // don't form a polygon.
  ASSERT((num_points <= 2) || 
         ((num_points > 2) && (*count > 2)));
}

void generate_face_node_conn(mesh_t* mesh,
                             int** face_nodes,
                             int* face_node_offsets)
{
  // Figure out face-node connectivity. We do this by computing centers
  // for all the cells and then using them to define face normals, 
  // from which node orderings can be determining using a convex hull 
  // determination algorithm (gift wrapping).

  // Make a list of all nodes attached to faces (unordered).
  int num_cells = mesh->num_cells;
  int num_faces = mesh->num_faces;
  int* nodes_for_face[num_faces];
  face_node_offsets[0] = 0;
  for (int f = 0; f < num_faces; ++f)
  {
    int rays = 0, ne = mesh->faces[f].num_edges;
    for (int e = 0; e < ne; ++e)
    {
      if (mesh->faces[f].edges[e]->node2 == NULL)
        ++rays;
    }
    int num_nodes = ne - rays;
    ASSERT(num_nodes >= 3);
    face_node_offsets[f+1] = face_node_offsets[f] + num_nodes;
    nodes_for_face[f] = malloc(num_nodes*sizeof(int));
  }

  int_avl_tree_t* fnodes = int_avl_tree_new();
  for (int f = 0; f < num_faces; ++f)
  {
    int counter = 0, ne = mesh->faces[f].num_edges;
    for (int e = 0; e < ne; ++e)
    {
      edge_t* edge = mesh->faces[f].edges[e];
      int node1_id = edge->node1 - &mesh->nodes[0];
      if (int_avl_tree_find(fnodes, node1_id) == NULL)
      {
        nodes_for_face[f][counter++] = node1_id;
        int_avl_tree_insert(fnodes, node1_id);
      }
      if (edge->node2 != NULL)
      {
        int node2_id = edge->node2 - &mesh->nodes[0];
        if (int_avl_tree_find(fnodes, node2_id) == NULL)
        {
          nodes_for_face[f][counter++] = node2_id;
          int_avl_tree_insert(fnodes, node2_id);
        }
      }
    }
    ASSERT(counter == face_node_offsets[f+1] - face_node_offsets[f]);
    int_avl_tree_clear(fnodes);
  }
  int_avl_tree_free(fnodes);

  // Compute cell centers from face nodes.
  point_t cell_centers[num_cells];
  memset(cell_centers, 0, num_cells*sizeof(point_t));
  int_avl_tree_t* cell_nodes = int_avl_tree_new();
  for (int c = 0; c < num_cells; ++c)
  {
    int num_nodes = 0;
    for (int f = 0; f < mesh->cells[c].num_faces; ++f)
    {
      int num_face_nodes = face_node_offsets[f+1] - face_node_offsets[f];
      for (int n = 0; n < num_face_nodes; ++n)
      {
        int node_id = nodes_for_face[f][n];
        if (int_avl_tree_find(cell_nodes, node_id) == NULL)
        {
          int_avl_tree_insert(cell_nodes, node_id);
          node_t* node = &mesh->nodes[nodes_for_face[f][n]];
          cell_centers[c].x += node->x;
          cell_centers[c].y += node->y;
          cell_centers[c].z += node->z;
          ++num_nodes;
        }
      }
    }
    cell_centers[c].x /= num_nodes;
    cell_centers[c].y /= num_nodes;
    cell_centers[c].z /= num_nodes;
    int_avl_tree_clear(cell_nodes);
  }
  int_avl_tree_free(cell_nodes);

  int_slist_t* all_face_nodes_list = int_slist_new();
  for (int f = 0; f < mesh->num_faces; ++f)
  {
    // Compute the normal vector for the face, pointing outward from 
    // its first cell.
    int nn = face_node_offsets[f+1] - face_node_offsets[f];
    int* nodes = nodes_for_face[f];
    ASSERT(nn >= 3);
    point_t face_center = {.x = 0.0, .y = 0.0, .z = 0.0};
    for (int n = 0; n < nn; ++n)
    {
      node_t* node = &mesh->nodes[nodes[n]];
      face_center.x += node->x;
      face_center.y += node->y;
      face_center.z += node->z;
    }
    face_center.x /= nn;
    face_center.y /= nn;
    face_center.z /= nn;

    // Construct vectors v1, v2, and v3, where v1 is the vector pointing from the 
    // face center to the first face node, v2 is a vector pointing from the face 
    // center to any other face, node, and v3 is their cross product.
    vector_t v1;
    v1.x = mesh->nodes[nodes[0]].x - face_center.x;
    v1.y = mesh->nodes[nodes[0]].y - face_center.y;
    v1.z = mesh->nodes[nodes[0]].z - face_center.z;

    vector_t v2, normal;
    double normal_mag;
    for (int n = 1; n < nn; ++n)
    {
      v2.x = mesh->nodes[nodes[n]].x - face_center.x;
      v2.y = mesh->nodes[nodes[n]].y - face_center.y;
      v2.z = mesh->nodes[nodes[n]].z - face_center.z;

      // normal = v1 x v2.
      vector_cross(&v1, &v2, &normal);
      normal_mag = sqrt(vector_dot(&normal, &normal));
      if (normal_mag > 1e-14) break;
    }
    ASSERT(normal_mag > 1e-14);
    normal.x /= normal_mag; normal.y /= normal_mag; normal.z /= normal_mag;

    vector_t v3;
    point_t cell_center;
    int cell1 = mesh->faces[f].cell1 - &mesh->cells[0];
    cell_center.x = cell_centers[cell1].x;
    cell_center.y = cell_centers[cell1].y;
    cell_center.z = cell_centers[cell1].z;
    v3.x = face_center.x - cell_center.x;
    v3.y = face_center.y - cell_center.y;
    v3.z = face_center.z - cell_center.z;
    if (vector_dot(&normal, &v3) < 0.0)
    {
      normal.x *= -1.0; normal.y *= -1.0; normal.z *= -1.0;
    }

    // Now project the coordinates of the face's nodes to the plane
    // with the given normal and centered about the face center.
    double points[2*nn]; // NOTE: planar coordinates (2D)
    vector_t e1, e2; // Basis vectors in the plane.
    double v1_mag = sqrt(vector_dot(&v1, &v1));
    e1.x = v1.x / v1_mag;
    e1.y = v1.y / v1_mag;
    e1.z = v1.z / v1_mag;

    // e2 = normal x e1.
    vector_cross(&normal, &e1, &e2);
    for (int p = 0; p < nn; ++p)
    {
      // v = node center - cell center.
      vector_t v;
      v.x = mesh->nodes[nodes[p]].x - cell_center.x;
      v.y = mesh->nodes[nodes[p]].y - cell_center.y;
      v.z = mesh->nodes[nodes[p]].z - cell_center.z;

      // Compute the perpendicular component of the point
      // with location v:
      // v_perp = v - (n o v)n.
      vector_t v_perp;
      double nov = vector_dot(&normal, &v);
      v_perp.x = v.x - nov * normal.x;
      v_perp.y = v.y - nov * normal.y;
      v_perp.z = v.z - nov * normal.z;

      // Project it to the plane.
      points[2*p]   = vector_dot(&v_perp, &e1);
      points[2*p+1] = vector_dot(&v_perp, &e2);
    }

    // Find the node order by traversing the convex hull of 
    // the points within the plane, appending them to face_nodes.
    int indices[nn], count;
    traverse_convex_hull(points, nn, indices, &count);
    face_node_offsets[f+1] = face_node_offsets[f] + nn;
    for (int n = 0; n < count; ++n)
      int_slist_append(all_face_nodes_list, nodes_for_face[f][indices[n]]);
    free(nodes_for_face[f]);
  }

  // Write the connectivity information.
  ASSERT(all_face_nodes_list->size == face_node_offsets[mesh->num_faces]);
  *face_nodes = malloc(sizeof(int)*all_face_nodes_list->size);
  int nf = all_face_nodes_list->size;
  for (int i = 0; i < nf; ++i)
    (*face_nodes)[i] = int_slist_pop(all_face_nodes_list);

  int_slist_free(all_face_nodes_list);
}

