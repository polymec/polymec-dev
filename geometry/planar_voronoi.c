#include <stdlib.h>
#include <string.h>
#include "geometry/voronoi.h"
#include "core/mesh.h"
#include "core/slist.h"
#include "core/avl_tree.h"

#define MIN(a, b) ((a > b) ? b : a)
#define MAX(a, b) ((a > b) ? a : b)

#ifdef __cplusplus
extern "C" {
#endif

// This goes here for protection from C++!
// Because Hang Si has messed with some of Shewchuk's predicates and
// included them with his own Tetgen library, we need to rename some of 
// the symbols therein to prevent duplicate symbols from confusing the 
// linker. Gross.
#define DTRILIBRARY 
#define REDUCED 
#define CDT_ONLY 
#define ANSI_DECLARATORS 
#define REAL double 
#define VOID void
#define exactinit triangle_exactinit 
#define fast_expansion_sum_zeroelem triangle_fast_expansion_sum_zeroelem 
#define scale_expansion_zeroelim triangle_scale_expansion_zeroelim 
#define estimate triangle_estimate 
#define orient3dadapt triangle_orient3dadapt 
#define orient3d triangle_orient3d 
#define incircle triangle_incircle
#include "triangle.h"

// This AVL node visitor appends the tree's data to a tag.
static void append_to_tag(avl_node_t* node, void* p)
{
  int* tag_p = (int*)p;
  *tag_p = (int)node->attribute;
  ++tag_p;
}

mesh_t* voronoi_plane(point_t* points, int num_points, 
                      point_t* ghost_points, int num_ghost_points)
{
  ASSERT(points != NULL);
  ASSERT(num_points >= 2);
  ASSERT(num_ghost_points >= 0);

  // Now we use Triangle to obtain a Voronoi graph.
  struct triangulateio in, delaunay, voro;

  // Define input points.
  in.numberofpoints = num_points;
  in.pointlist = malloc(num_points*sizeof(double));
  for (int i = 0; i < num_points; ++i)
  {
    in.pointlist[2*i+0] = points[i].x;
    in.pointlist[2*i+1] = points[i].y;
  }

  // No point attributes or markers.
  in.numberofpointattributes = 0;
  in.pointattributelist = NULL; 
  in.pointmarkerlist = NULL; // No point markers.

  // No segments or holes.
  in.numberofsegments = 0;
  in.segmentlist = NULL;
  in.segmentmarkerlist = NULL;
  in.numberofholes = 0;
  in.holelist = NULL;

  // No regions.
  in.numberofregions = 0;
  in.regionlist = NULL;

  // Set up the structure for the triangulation.
  delaunay.pointlist = NULL;
  delaunay.pointattributelist = NULL;
  delaunay.pointmarkerlist = NULL;
  delaunay.trianglelist = NULL;
  delaunay.triangleattributelist = NULL;
  delaunay.neighborlist = NULL;
  delaunay.segmentlist = NULL;
  delaunay.segmentmarkerlist = NULL;
  delaunay.edgelist = NULL;
  delaunay.edgemarkerlist = NULL;
  delaunay.holelist = NULL;

  // Set up the structure that will hold the Voronoi graph.
  voro.pointlist = NULL;
  voro.edgelist = NULL;
  voro.normlist = NULL;

  // Do the triangulation. Switches pass to triangle are:
  // -Q : Quiet (shaddap!), no output on the terminal except errors.
  // -z : Indices are all numbered from zero.
  // -e : Generates edges and places them in out.edgelist.
  // -D : Forces a true delaunay triangulation, which makes the Voronoi graph valid.
  // -v : Generates the Voronoi graph, storing it in voro.
  triangulate((char*)"QzeDv", &in, &delaunay, &voro);

  // Make sure we got something.
  if (delaunay.numberoftriangles == 0)
    arbi_error("voronoi_plane: delaunay triangulation produced 0 triangles!");
  if (delaunay.numberofpoints != num_points)
  {
    char err[1024];
    snprintf(err, 1024, "voronoi_plane: delaunay triangulation produced %d triangles\n(%d generating points given)", 
             delaunay.numberofpoints, num_points);
    arbi_error(err);
  }
  if (voro.numberofpoints <= 0)
    arbi_error("voronoi_plane: Error occurred generating Voronoi graph.");

  // Now assemble the mesh.
  int num_cells = num_points;
  int num_ghost_cells = num_ghost_points;
  int num_nodes = 2*voro.numberofpoints;
  int horiz_edges = 2*delaunay.numberofedges;
  int vert_edges = voro.numberofpoints;
  int num_edges = horiz_edges + vert_edges;
  int horiz_faces = delaunay.numberofedges;
  int vert_faces = 2*num_points;
  int num_faces = horiz_faces + vert_faces;
  mesh_t* mesh = mesh_new(num_cells, num_ghost_cells,
                          num_faces, num_edges, num_nodes);

  // Node coordinates.
  for (int n = 0; n < num_nodes/2; ++n)
  {
    mesh->nodes[2*n].x = voro.pointlist[2*n+0];
    mesh->nodes[2*n].y = voro.pointlist[2*n+1];
    mesh->nodes[2*n].z = -0.5;
    mesh->nodes[2*n+1].x = voro.pointlist[2*n+0];
    mesh->nodes[2*n+1].y = voro.pointlist[2*n+1];
    mesh->nodes[2*n+1].z = 0.5;
  }

  // Edge <-> node connectivity.

  // Horizontal edges: 
  // mesh->edges[2*i]   <-- ith bottom edge
  // mesh->edges[2*i+1] <-- ith top edge
  // NOTE: We keep track of "outer" edges and tag them accordingly.
  avl_tree_t* outer_edges = int_avl_tree_new();
  int num_outer_edges = 0;
  for (int i = 0; i < delaunay.numberofedges; ++i)
  {
    mesh->edges[2*i].node1 = &mesh->nodes[2*voro.edgelist[2*i]];
    mesh->edges[2*i+1].node1 = &mesh->nodes[2*voro.edgelist[2*i]+1];
    int n2 = voro.edgelist[2*i+1]; // -1 if ghost
    if (n2 == -1)
    {
      int e1 = 2*i, e2 = 2*i+1;
      avl_tree_insert(outer_edges, (void*)e1);
      avl_tree_insert(outer_edges, (void*)e2);
      ++num_outer_edges;
      mesh->edges[2*i].node2 = NULL;
      mesh->edges[2*i+1].node2 = NULL;
    }
    else
    {
      mesh->edges[2*i].node2 = &mesh->nodes[2*n2];
      mesh->edges[2*i+1].node2 = &mesh->nodes[2*n2+1];
    }
  }
  // Vertical edges.
  // mesh->edges[horiz_edges+2*i] <-- ith vertical edge
  for (int i = 0; i < num_nodes; ++i)
  {
    mesh->edges[horiz_edges+2*i].node1 = &mesh->nodes[2*i];
    mesh->edges[horiz_edges+2*i].node2 = &mesh->nodes[2*i+1];
  }
  // Tag the outer edges as such.
  if (num_outer_edges > 0)
  {
    int* outer_edge_tag = mesh_create_tag(mesh->edge_tags, "outer_edges", num_outer_edges);
    avl_node_t* root = avl_tree_root(outer_edges);
    int* tag_p = outer_edge_tag;
    avl_node_visit(root, &append_to_tag, (void*)tag_p);

    // Stick the outward rays in a "rays" property on this tag.
    double* rays = malloc(3*num_outer_edges*sizeof(double));
    mesh_tag_set_property(mesh->edge_tags, "outer_edges", "rays", rays, &free);
    for (int i = 0; i < num_outer_edges; ++i)
    {
      int j = outer_edge_tag[i];
      ASSERT(voro.edgelist[2*i+1] == -1);
      rays[3*i+0] = voro.normlist[2*j];
      rays[3*i+1] = voro.normlist[2*j+1];
      rays[3*i+2] = 0.0;
    }
  }

  // Face <-> edge/cell connectivity.

  // Horizontal faces:
  // mesh->faces[i] <-- ith horizontal face
  avl_tree_t* outer_cells = int_avl_tree_new();
  int num_outer_cells = 0;
  for (int i = 0; i < delaunay.numberofedges; ++i)
  {
    // Hook up this face to its cells.
    int cell1 = delaunay.edgelist[2*i];
    int cell2 = delaunay.edgelist[2*i+1];
    mesh->faces[i].cell1 = &mesh->cells[cell1];
    mesh->cells[cell1].num_faces++;
    mesh->faces[i].cell2 = &mesh->cells[cell2];
    mesh->cells[cell2].num_faces++;

    // Horizontal faces are quadrilaterals.
    mesh->faces[i].num_edges = 4;
    mesh->faces[i].edges = malloc(4*sizeof(edge_t*));
    int e1 = 2*i, e2 = 2*i+1; // Horizontal edges
    int e3 = horiz_edges + voro.edgelist[2*i], // Vertical edges
        e4 = horiz_edges + voro.edgelist[2*i+1];
    mesh->faces[i].edges[0] = &mesh->edges[e1];
    mesh->faces[i].edges[1] = &mesh->edges[e2];
    mesh->faces[i].edges[2] = &mesh->edges[e3];
    mesh->faces[i].edges[3] = &mesh->edges[e4];

    // If one of these edges is an outer edge, then both cells that share 
    // the face are outer cells.
    if ((avl_tree_find(outer_edges, (void*)e1) != NULL) || 
        (avl_tree_find(outer_edges, (void*)e2) != NULL) || 
        (avl_tree_find(outer_edges, (void*)e3) != NULL) || 
        (avl_tree_find(outer_edges, (void*)e4) != NULL))
    {
      if (avl_tree_find(outer_cells, (void*)cell1) == NULL)
      {
        avl_tree_insert(outer_cells, (void*)cell1);
        ++num_outer_cells;
      }
      if (avl_tree_find(outer_cells, (void*)cell2) == NULL)
      {
        avl_tree_insert(outer_cells, (void*)cell2);
        ++num_outer_cells;
      }
    }
  }

  // Hook up the horizontal faces to cells, leaving room for top and 
  // bottom faces.
  for (int i = 0; i < mesh->num_cells; ++i)
  {
    mesh->cells[i].num_faces += 2;
    mesh->cells[i].faces = malloc(mesh->cells[i].num_faces*sizeof(face_t*));
  }
  int cell_next_face[mesh->num_cells];
  memset(cell_next_face, 0, sizeof(int)*mesh->num_cells);
  for (int i = 0; i < delaunay.numberofedges; ++i)
  {
    int cell1 = delaunay.edgelist[2*i];
    int cell2 = delaunay.edgelist[2*i+1];
    mesh->cells[cell1].faces[cell_next_face[cell1]++] = &mesh->faces[i];
    mesh->cells[cell2].faces[cell_next_face[cell2]++] = &mesh->faces[i];
  }

  // Vertical faces (arbitrary polygons).
  // mesh->faces[horiz_faces+2*i]   <-- bottom face for ith cell
  // mesh->faces[horiz_faces+2*i+1] <-- top face for ith cell
  int vface_num_edges[2*delaunay.numberofedges],
      vface_next_edge[2*delaunay.numberofedges];
  memset(vface_num_edges, 0, sizeof(int)*delaunay.numberofedges);
  memset(vface_next_edge, 0, sizeof(int)*delaunay.numberofedges);
  for (int i = 0; i < delaunay.numberofedges; ++i)
  {
    vface_num_edges[2*voro.edgelist[2*i]]++;
    vface_num_edges[2*voro.edgelist[2*i]+1]++;
    vface_num_edges[2*voro.edgelist[2*i+1]]++;
    vface_num_edges[2*voro.edgelist[2*i+1]+1]++;
  }
  // Hook up the top and bottom faces to cells.
  for (int i = 0; i < mesh->num_cells; ++i)
  {
    int bottom_face = horiz_faces + 2*i;
    int top_face = horiz_faces + 2*i + 1;
    mesh->cells[i].faces[cell_next_face[i]++] = &mesh->faces[bottom_face];
    mesh->cells[i].faces[cell_next_face[i]++] = &mesh->faces[top_face];
    mesh->faces[bottom_face].cell1 = &mesh->cells[i];
    mesh->faces[bottom_face].cell2 = NULL;
    mesh->faces[top_face].cell1 = &mesh->cells[i];
    mesh->faces[top_face].cell2 = NULL;
  }
  for (int i = 0; i < delaunay.numberofedges; ++i)
  {
    int edge1 = voro.edgelist[2*i], edge2 = voro.edgelist[2*i+1];

    // Bottom face - edge1
    if (mesh->faces[horiz_faces+2*edge1].edges == NULL)
      mesh->faces[horiz_faces+2*edge1].edges = malloc(vface_num_edges[edge1]*sizeof(edge_t*));
    mesh->faces[horiz_faces+2*edge1].edges[vface_next_edge[2*edge1]++] = &mesh->edges[2*edge1];

    // Top face - edge1
    if (mesh->faces[horiz_faces+2*edge1+1].edges == NULL)
      mesh->faces[horiz_faces+2*edge1+1].edges = malloc(vface_num_edges[edge1]*sizeof(edge_t*));
    mesh->faces[horiz_faces+2*edge1+1].edges[vface_next_edge[2*edge1+1]++] = &mesh->edges[2*edge1+1];

    // Bottom face - edge2
    if (mesh->faces[horiz_faces+2*edge2].edges == NULL)
      mesh->faces[horiz_faces+2*edge2].edges = malloc(vface_num_edges[edge2]*sizeof(edge_t*));
    mesh->faces[horiz_faces+2*edge2].edges[vface_next_edge[2*edge1]++] = &mesh->edges[2*edge1];

    // Top face - edge2
    if (mesh->faces[horiz_faces+2*edge2+1].edges == NULL)
      mesh->faces[horiz_faces+2*edge2+1].edges = malloc(vface_num_edges[edge2]*sizeof(edge_t*));
    mesh->faces[horiz_faces+2*edge2+1].edges[vface_next_edge[2*edge1+1]++] = &mesh->edges[2*edge1+1];
  }

  // Tag the outer cells as such.
  ASSERT(num_outer_cells > 0);
  int* outer_cell_tag = mesh_create_tag(mesh->cell_tags, "outer_cells", num_outer_cells);
  avl_node_t* root = avl_tree_root(outer_cells);
  int* tag_p = outer_cell_tag;
  avl_node_visit(root, &append_to_tag, (void*)tag_p);

  // Finally, we create properties on the outer_edges and outer_cells tags 
  // that associate one with the other.
  slist_t* outer_cell_edges = slist_new(NULL);
  for (int i = 0; i < num_outer_cells; ++i)
  {
    int num_edges = 0;
    slist_node_t* pos = slist_back(outer_cell_edges);
    for (int f = 0; f < mesh->cells[i].num_faces; ++f)
    {
      for (int e = 0; e < mesh->cells[i].faces[f]->num_edges; ++e)
      {
        int edgeid = (int)(mesh->faces[f].edges[e] - &mesh->edges[0]);
        if (avl_tree_find(outer_edges, (void*)edgeid) != NULL)
        {
          slist_append(outer_cell_edges, (void*)edgeid);
          ++num_edges;
        }
      }
    }
    pos = pos->next;
    slist_insert(outer_cell_edges, pos, (void*)num_edges);
  }

  // Add 'outer_edges' as a property of the outer_cells.
  int* oce = malloc(slist_size(outer_cell_edges)*sizeof(double));
  mesh_tag_set_property(mesh->edge_tags, "outer_cells", "outer_edges", outer_cell_edges, &free);
  int offset = 0;
  for (slist_node_t* n = slist_front(outer_cell_edges); n != NULL;)
  {
    // Read the number of edges for the cell.
    int num_edges = (int)n->value;
    n = n->next;
    oce[offset++] = num_edges;
    for (int e = 0; e < num_edges; ++e)
    {
      oce[offset++] = (int)n->value;
      n = n->next;
    }
  }

  // Clean up.
  slist_free(outer_cell_edges);
  avl_tree_free(outer_cells);
  avl_tree_free(outer_edges);
  free(in.pointlist);
  free(delaunay.pointlist);
  free(delaunay.trianglelist);
  free(delaunay.edgelist);
  free(voro.pointlist);
  free(voro.edgelist);
  free(voro.normlist);

  return mesh;
}

mesh_t* voronoi_stack(mesh_t* plane, int num_planes, double z1, double z2)
{
  int num_points = 0; // FIXME
  int num_ghost_points = 0;
  int num_cells = num_planes * num_points;
  int num_ghost_cells = num_planes * num_ghost_points;
  int nodes_per_plane = 0;
  int num_nodes = (num_planes + 1) * nodes_per_plane;
  int horiz_edges_per_plane = 0;
  int vert_edges_per_plane = 0;
  int num_edges = num_planes * horiz_edges_per_plane + 
                  (num_planes - 1) * vert_edges_per_plane;
  int horiz_faces_per_plane = 0;
  int vert_faces_per_plane = num_points;
  int num_faces = num_planes * horiz_faces_per_plane + 
                  (num_planes + 1) * vert_faces_per_plane;
  mesh_t* mesh = mesh_new(num_cells, num_ghost_cells,
                          num_faces, num_edges, num_nodes);
  return mesh;
}

#ifdef __cplusplus
}
#endif

