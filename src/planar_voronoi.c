#include <stdlib.h>
#include "voronoi.h"
#include "mesh.h"

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

//------------------------------------------------------------------------
planar_voronoi_t* planar_voronoi_new(point_t* points, int num_points, bbox_t* bounding_box)
{
  ASSERT(points != NULL);
  ASSERT(num_points >= 2);
  ASSERT(bounding_box != NULL);

  // Now we use Triangle to obtain a Voronoi graph.
  struct triangulateio in, delaunay, voro;

  // Define input points. Include the 4 bounding box corners.
  in.numberofpoints = num_points + 4;
  in.pointlist = malloc(num_points*sizeof(double));
  in.pointlist[0] = bounding_box->x1;
  in.pointlist[1] = bounding_box->y1;
  in.pointlist[2] = bounding_box->x2;
  in.pointlist[3] = bounding_box->y1;
  in.pointlist[4] = bounding_box->x1;
  in.pointlist[5] = bounding_box->y2;
  in.pointlist[6] = bounding_box->x2;
  in.pointlist[7] = bounding_box->y2;
  for (int i = 4; i < num_points + 4; ++i)
  {
    in.pointlist[2*i+0] = points[i].x;
    in.pointlist[2*i+1] = points[i].y;
  }

  // No point attributes or markers.
  in.numberofpointattributes = 0;
  in.pointattributelist = 0; 
  in.pointmarkerlist = 0; // No point markers.

  // No segments or holes.
  in.numberofsegments = 0;
  in.segmentlist = 0;
  in.segmentmarkerlist = 0;
  in.numberofholes = 0;
  in.holelist = 0;

  // No regions.
  in.numberofregions = 0;
  in.regionlist = 0;

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
  // -c : Generates convex hull and places it in out.segmentlist.
  // -v : Generates the Voronoi graph, storing it in voro.
  triangulate((char*)"Qzecv", &in, &delaunay, &voro);

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

  planar_voronoi_t* v = NULL;//mesh_new(num_cells, num_ghost_cells, num_faces, num_edges, num_nodes);
  return v;
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
mesh_t* extrusion(planar_voronoi_t* planar_graph, int num_extruded_cells, double length)
{
  ASSERT(num_extruded_cells > 0);
  ASSERT(length > 0.0);
  arbi_not_implemented("extrusion");
  return NULL;
}
//------------------------------------------------------------------------

#ifdef __cplusplus
}
#endif

