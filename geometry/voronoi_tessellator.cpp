// Welcome to tessellator.cpp, one of the only C++ files in the Polymec source 
// code. This code uses Tetgen, which is a C++ library for creating Delaunay 
// tetrahedralizations of domains. Since C++ and C are actually very different 
// languages, this code is insulated from the rest of Polymec in terms of 
// data types.

#include "core/polymec.h"
#include <gc/gc.h>
#include "geometry/voronoi_tessellator.h"

#define TETLIBRARY
#include "tetgen.h"

extern "C" {

void voronoi_tessellation_free(voronoi_tessellation_t* tessellation)
{
  for (int c = 0; c < tessellation->num_cells; ++c)
    free(tessellation->cells[c].faces);
  free(tessellation->cells);

  for (int f = 0; f < tessellation->num_faces; ++f)
    free(tessellation->faces[f].edges);
  free(tessellation->faces);

  free(tessellation->edges);
  free(tessellation->nodes);
  free(tessellation);
}

struct voronoi_tessellator_t
{
};

voronoi_tessellator_t* voronoi_tessellator_new()
{
  voronoi_tessellator_t* t = (voronoi_tessellator_t*)GC_MALLOC(sizeof(voronoi_tessellator_t));
  return t;
}

namespace 
{

voronoi_tessellation_t* voronoi_tessellation_new(int num_cells, int num_faces, int num_edges, int num_nodes)
{
  voronoi_tessellation_t* t = (voronoi_tessellation_t*)malloc(sizeof(voronoi_tessellation_t));
  t->num_cells = num_cells;
  t->cells = (voronoi_cell_t*)malloc(num_cells*sizeof(voronoi_cell_t));
  t->num_faces = num_faces;
  t->faces = (voronoi_face_t*)malloc(num_faces*sizeof(voronoi_face_t));
  t->num_edges = num_edges;
  t->edges = (voronoi_edge_t*)malloc(num_edges*sizeof(voronoi_edge_t));
  t->num_nodes = num_nodes;
  t->nodes = (double*)malloc(3*num_nodes*sizeof(double));
  return t;
}

}

voronoi_tessellation_t*
voronoi_tessellator_tessellate(voronoi_tessellator_t* tessellator,
                               double* points, int num_points)
{
  ASSERT(points != NULL);
  ASSERT(num_points >= 2);

  // Set Tetgen's options. We desire a Voronoi mesh.
  tetgenio in;
  in.initialize();
  in.numberofpoints = num_points;
  in.pointlist = new double[sizeof(double)*3*in.numberofpoints];
  memcpy(in.pointlist, points, sizeof(double)*3*in.numberofpoints);

  // Tetrahedralize. Command line options are:
  // Q          -- Quiet, no output to terminal (shaddap, Tetgen!).
  // v          -- Generate a Voronoi tessellation.
  // B, N, E, F -- Suppress the generation of boundary, node, edge, face files.
  // C          -- Perform a consistency check on the final mesh.
  tetgenio out;
  out.initialize();
  tetrahedralize((char*)"QvBNEFC", &in, &out, NULL, NULL);
  ASSERT(out.numberofvcells == num_points);

  // Copy stuff to a fresh tessellation object.
  voronoi_tessellation_t* t = voronoi_tessellation_new(out.numberofvcells, 
                                                       out.numberofvfacets, 
                                                       out.numberofvedges, 
                                                       out.numberofvpoints);
  memcpy(t->nodes, out.vpointlist, sizeof(double)*3*t->num_nodes);

  // Edge <-> node connectivity.
  for (int i = 0; i < t->num_edges; ++i)
  {
    t->edges[i].node1 = out.vedgelist[i].v1;
    t->edges[i].node2 = out.vedgelist[i].v2;
    if (t->edges[i].node2 == -1) // node2 is a "ghost"
    {
      t->edges[i].ray[0] = out.vedgelist[i].vnormal[0];
      t->edges[i].ray[1] = out.vedgelist[i].vnormal[1];
      t->edges[i].ray[2] = out.vedgelist[i].vnormal[2];
    }
  }

  // Face <-> edge/cell connectivity.
  for (int i = 0; i < t->num_faces; ++i)
  {
    t->faces[i].cell1 = out.vfacetlist[i].c1;
    t->faces[i].cell2 = out.vfacetlist[i].c2;
    int Ne = out.vfacetlist[i].elist[0];
    t->faces[i].num_edges = Ne;
    t->faces[i].edges = (int*)malloc(sizeof(int)*Ne);
    for (int j = 0; j < Ne; ++j)
      t->faces[i].edges[j] = out.vfacetlist[i].elist[j+1];
  }

  // Cell <-> face connectivity.
  for (int i = 0; i < t->num_cells; ++i)
  {
    int Nf = out.vcelllist[i][0];
    t->cells[i].num_faces = Nf;
    t->cells[i].faces = (int*)malloc(sizeof(int)*Nf);
    for (int f = 0; f < Nf; ++f)
      t->cells[i].faces[f] = out.vcelllist[i][f+1];
  }

  return t;
}

}

