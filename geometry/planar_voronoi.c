#include <stdlib.h>
#include "geometry/voronoi.h"
#include "core/mesh.h"

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

// This represents an edge in a planar Voronoi graph.
typedef struct planar_voronoi_cell_t planar_voronoi_cell_t;
typedef struct 
{
  node_t* node1;
  node_t* node2;
  planar_voronoi_cell_t* cell1;
  planar_voronoi_cell_t* cell2;
} planar_voronoi_edge_t;

struct planar_voronoi_cell_t
{
  int num_edges;
  planar_voronoi_edge_t** edges;
};

struct planar_voronoi_t
{
  int num_cells;
  planar_voronoi_cell_t* cells;
  int num_edges;
  planar_voronoi_edge_t* edges;
  int num_nodes;
  node_t* nodes;
};

static planar_voronoi_t* planar_voronoi_new(int num_cells, int num_edges, int num_nodes)
{
  planar_voronoi_t* v = malloc(sizeof(planar_voronoi_t));
  v->num_cells = num_cells;
  v->cells = malloc(num_cells*sizeof(planar_voronoi_cell_t));
  v->num_edges = num_edges;
  v->edges = malloc(num_edges*sizeof(planar_voronoi_edge_t));
  v->num_nodes = num_nodes;
  v->nodes = malloc(num_nodes*sizeof(node_t));
  return v;
}

static void planar_voronoi_free(planar_voronoi_t* v)
{
  free(v->cells);
  free(v->edges);
  free(v->nodes);
  free(v);
}

planar_voronoi_t* planar_voronoi_graph(point_t* points, int num_points, sp_func_t* F)
{
  ASSERT(points != NULL);
  ASSERT(num_points >= 2);
  ASSERT(F != NULL);
  ASSERT(!sp_func_is_homogeneous(F));

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
  // -D : Forces a true Delauney triangulation, which makes the Voronoi graph valid.
  // -v : Generates the Voronoi graph, storing it in voro.
  triangulate((char*)"QzeDv", &in, &delaunay, &voro);

  // Make sure we got something.
  if (delaunay.numberoftriangles == 0)
    arbi_error("planar_voronoi_new: Delauney triangulation produced 0 triangles!");
  if (delaunay.numberofpoints != num_points)
  {
    char err[1024];
    snprintf(err, 1024, "planar_voronoi_new: Delauney triangulation produced %d triangles\n(%d generating points given)", 
             delaunay.numberofpoints, num_points);
    arbi_error(err);
  }
  if (voro.numberofpoints <= 0)
    arbi_error("planar_voronoi_new: Error occurred generating Voronoi graph.");

  // Now assemble the planar Voronoi graph.
  int num_cells = num_points;
  int num_edges = delaunay.numberofedges;
  int num_nodes = voro.numberofpoints;

  // Count up the semi-infinite edge rays, add corresponding boundary nodes,
  // and count the number of boundary edges we need.
  int num_bnodes = 0, num_bedges = 0;
  int cell_bedge[num_cells];
  memset(cell_needs_edge, -1, num_cells*sizeof(int));
  for (int e = 0; e < num_edges; ++e)
  {
    if (voro.edgelist[2*e+1] == -1)
    {
      ++num_bnodes;
      int cell1 = delaunay.edgelist[2*e],
          cell2 = delaunay.edgelist[2*e+1];
      if (!cell_needs_edge[cell1] && !cell_needs_edge[cell2])
      {
        cell_bedge[cell1] = cell_bedge[cell2] = num_bedges;
        ++num_bedges;
      }
    }
  }

  // Construct the graph.
  int num_ghost_cells = num_bedges;
  planar_voronoi_t* v = planar_voronoi_new(num_cells + num_bedges, 
                                           num_edges + num_bedges, 
                                           num_nodes + num_bnodes);
  for (int n = 0; n < num_nodes; ++n)
  {
    v->nodes[n].x = voro.pointlist[2*n+0];
    v->nodes[n].y = voro.pointlist[2*n+1];
    v->nodes[n].z = 0.0;
  }
  for (int n = num_nodes; n < num_nodes + num_bnodes; ++n)
  {
    // Trace back from the given node to the boundary along the
    // gradient of the implicit function F.
    // FIXME
    v->nodes[n].x = 0.0, v->nodes[n].y = 0.0;
  }

  int bnode = 0;
  for (int e = 0; e < num_edges; ++e)
  {
    // Hook up the nodes.
    v->edges[e].node1 = &(v->nodes[voro.edgelist[2*e]]);
    int node2 = voro.edgelist[2*e+1];
    if (node2 == -1) // Semi-infinite ray!
    {
      // This cell is adjacent to the boundary.
      v->edges[e].node2 = num_nodes + bnode;
      v->edges[e].cell2 = NULL;
      ++bnode;
    }
    else
    {
      v->edges[e].node2 = &(v->nodes[node2]);
    }

    // Hook up the cells to the edges and vice versa.
    int cell1 = MIN(delaunay.edgelist[2*e], delaunay.edgelist[2*e+1]);
    int cell2 = MAX(delaunay.edgelist[2*e], delaunay.edgelist[2*e+1]);
    v->cells[cell1].edges[v->cells[cell1].num_edges] = &(v->edges[e]);
    ++(v->cells[cell1].num_edges);
    v->cells[cell2].edges[v->cells[cell2].num_edges] = &(v->edges[e]);
    ++(v->cells[cell2].num_edges);
  }

  // Add any needed boundary edges.
  for (int c = 0; c < num_cells; ++c)
  {
    int be = cell_bedges[c];
    if (be >= 0)
    {
      v->cells[c].edges[v->cells[cell1].num_edges] = &(v->edges[be]);
      ++(v->cells[c].num_edges);
      v->edges[be].cell1 = &(v->cells[c]);
      v->edges[be].cell2 = &(v->cells[be]);
    }
  }

  // Clean up our Triangle data structures.
  free(in.pointlist);
  free(delaunay.pointlist);
  free(delaunay.trianglelist);
  free(delaunay.edgelist);
  free(voro.pointlist);
  free(voro.edgelist);
  free(voro.normlist);

  return v;
}

mesh_t* voronoi_extrusion(planar_voronoi_t* planar_graph, int num_layers, double length)
{
  ASSERT(planar_graph != NULL);
  ASSERT(num_extruded_cells > 0);
  ASSERT(length > 0.0);
  arbi_not_implemented("extrusion");

  // Dispose of the planar graph.
  planar_voronoi_free(planar_graph);

  return NULL;
}

#ifdef __cplusplus
}
#endif

