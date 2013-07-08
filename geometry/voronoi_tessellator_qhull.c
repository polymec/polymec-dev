// This implementation of the Voronoi tessellator uses QHull.

#define qh_QHpointer (1)
#include "libqhull.h"

#include "core/polymec.h"
#include <gc/gc.h>
#include "geometry/voronoi_tessellator.h"

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
  // Temporary file where QHull data is stored.
  char qhull_filename[1024];
};

static void voronoi_tessellator_free(void* ctx, void* dummy)
{
  voronoi_tessellator_t* t = ctx;
  unlink(t->qhull_filename);
}

voronoi_tessellator_t* voronoi_tessellator_new()
{
  voronoi_tessellator_t* t = (voronoi_tessellator_t*)GC_MALLOC(sizeof(voronoi_tessellator_t));
  GC_register_finalizer(t, voronoi_tessellator_free, t, NULL, NULL);
  strcpy(t->qhull_filename, "QHULLXXXXXXXX");
  mkstemp(t->qhull_filename);
  return t;
}

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

voronoi_tessellation_t*
voronoi_tessellator_tessellate(voronoi_tessellator_t* tessellator,
                               double* points, int num_points)
{
  ASSERT(points != NULL);
  ASSERT(num_points >= 2);

  // Set up QHull to tessellate the points. Much of this was borrowed from 
  // user_eg.c in the QHull distribution.
  char command[128];
  sprintf(command, "qhull v Qbb");
  int status = qh_new_qhull(3, num_points, points, False,
                            command, NULL, NULL);
  if (status != 0)
    polymec_error("Could not create an instance of QHull for Voronoi tessellation.");

  // Determine the numbers of Voronoi cells and nodes.
  qh_findgood_all(qh facet_list);
  int num_cells = qh num_vertices - qh_setsize(qh del_vertices);
  int num_nodes = qh num_good;

  // Computes Voronoi centers for all facets. 
  qh_setvoronoi_all();
  facetT* facet;
  vertexT* vertex;

  // NOTE: For edges, see the QHull functions qh_facet3vertex, qh_nextridge3d.
  // NOTE: qh_facet3vertex returns a setT of vertices, and sets are described 
  // NOTE: in qset.h (qh_setsize(set) returns the set's size, for example).

  // Count the faces and edges.
  int num_faces = 0, num_edges = 0;
  FORALLfacets
  {
    facet->seen = false;
    ++num_faces;
    // FIXME
  }

  // Find the numbers of neighbors for each cell.
  int* num_neighbors = malloc(sizeof(int) * num_cells);

  int i = 0;
  FORALLvertices
  {
    qh_order_vertexneighbors(vertex);

    bool infinity_seen = false;
    facetT *neighbor, **neighbor_ptr;
    FOREACHneighbor_(vertex)
    {
      if (neighbor->upperdelaunay)
      {
        if (!infinity_seen)
        {
          infinity_seen = true;
//          num_neighbors[i]++;
        }
      }
      else
      {
        neighbor->seen = true;
        num_neighbors[i]++;
      }
    }
    ++i;
  }

  // Now we record the node coordinates.
  FORALLfacets
  {
    facet->seen = false;
  }

  i = 0;
  FORALLvertices
  {
    qh_order_vertexneighbors(vertex);
    bool infinity_seen = false;
    int index = qh_pointid(vertex->point);
    int nn = num_neighbors

    // Skip those cells with only one neighbor.
    // (This is an oddity of QHull, apparently.)
    if (nn == 1) continue;
    facetT *neighbor, **neighbor_ptr;

  }

  // Copy stuff to a fresh tessellation object.
  voronoi_tessellation_t* t = NULL;//voronoi_tessellation_new(out.numberofvcells, 
                                   //                    out.numberofvfacets, 
                                   //                    out.numberofvedges, 
                                   //                    out.numberofvpoints);
//  memcpy(t->nodes, out.vpointlist, sizeof(double)*3*t->num_nodes);

#if 0
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
#endif

  // Clean up.
#ifdef qh_NOmem
  qh_freeqhull( True);
#else
  qh_freeqhull( False);
  int curlong, totlong;
  qh_memfreeshort(&curlong, &totlong);
  if (curlong || totlong)
    fprintf(stderr, "qhull internal warning (main): did not free %d bytes of long memory(%d pieces)\n",
       totlong, curlong);
#endif
  return t;
}

