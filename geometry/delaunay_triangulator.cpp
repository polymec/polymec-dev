#include "core/polymec.h"
#include <gc/gc.h>
#include "geometry/delaunay_triangulator.h"

#define TETLIBRARY
#include "tetgen.h"

extern "C" {

delaunay_triangulation_t* delaunay_triangulation_new(int num_tets, int num_vertices)
{
  delaunay_triangulation_t* t = (delaunay_triangulation_t*)malloc(sizeof(delaunay_triangulation_t));
  t->num_tets = num_tets;
  t->tets = (delaunay_tet_t*)malloc(num_tets*sizeof(delaunay_tet_t));
  t->num_vertices = num_vertices;
  t->vertices = (double*)malloc(3*num_vertices*sizeof(double));
  return t;
}

void delaunay_triangulation_free(delaunay_triangulation_t* triangulation)
{
  free(triangulation->tets);
  free(triangulation->vertices);
}

struct delaunay_triangulator_t 
{
};

delaunay_triangulator_t* delaunay_triangulator_new()
{
  delaunay_triangulator_t* t = (delaunay_triangulator_t*)GC_MALLOC(sizeof(delaunay_triangulator_t));
  return t;
}

delaunay_triangulation_t* 
delaunay_triangulator_triangulate(delaunay_triangulator_t* triangulator,
                                  double* points, int num_points)
{
  UNUSED_ARG(triangulator);
  ASSERT(points != NULL);
  ASSERT(num_points >= 4);

  // Set Tetgen's options. Real Delaunay triangulations only, please.
  tetgenio in;
  in.initialize();
  in.numberofpoints = num_points;
  in.pointlist = new double[sizeof(double)*3*num_points];
  memcpy(in.pointlist, points, sizeof(double)*3*num_points);

  // Tetrahedralize. Command line options are:
  // Q          -- Quiet, no output to terminal (shaddap, Tetgen!).
  // B, N, E, F -- Suppress the generation of boundary, node, edge, face files.
  // C          -- Perform a consistency check on the final mesh.
  tetgenio out;
  out.initialize();
  tetrahedralize((char*)"QBNEFC", &in, &out, NULL, NULL);
  ASSERT(out.numberofcorners == 4); // Linear tetrahedra only.

  // Copy stuff to a fresh triangulation object.
  delaunay_triangulation_t* t = delaunay_triangulation_new(out.numberoftetrahedra, 
                                                           out.numberofpoints);
  memcpy(t->vertices, out.pointlist, sizeof(double)*3*t->num_vertices);
  for (int i = 0; i < t->num_tets; ++i)
  {
    t->tets[i].vertices[0] = out.tetrahedronlist[4*i];
    t->tets[i].vertices[1] = out.tetrahedronlist[4*i+1];
    t->tets[i].vertices[2] = out.tetrahedronlist[4*i+2];
    t->tets[i].vertices[3] = out.tetrahedronlist[4*i+3];
  }
  
  return t;
}

}

