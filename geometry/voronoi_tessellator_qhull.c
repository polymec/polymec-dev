// This implementation of the Voronoi tessellator uses QHull.

#include "libqhull/libqhull.h"
#include "libqhull/mem.h"
#include "libqhull/qset.h"

#include "core/polymec.h"
#include "core/unordered_map.h"
#include "core/tuple.h"
#include "core/slist.h"
#include <gc/gc.h>
#include "geometry/voronoi_tessellator.h"

// We use fmemopen and open_memstream to communicate with QHull.
#ifdef __APPLE__
#include "fmemopen.h"
#endif
#include <stdio.h>

static char hidden_options[]=" d v H Qbb Qf Qg Qm Qr Qu Qv Qx Qz TR E V Fp Gt Q0 Q1 Q2 Q3 Q4 Q5 Q6 Q7 Q8 Q9 ";

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

voronoi_tessellation_t*
voronoi_tessellator_tessellate(voronoi_tessellator_t* tessellator,
                               double* points, int num_points)
{
  ASSERT(points != NULL);
  ASSERT(num_points >= 2);

  // Set up input and output streams to mimic stdin and stdout.
  char *input, *output;
  int input_size;
  FILE* fin = fmemopen(input, input_size, "r");
  size_t output_size;
  FILE* fout = open_memstream(&output, &output_size);

  // Initialize QHull.
  int argc = 5;
  char* argv[] = {"qvoronoi", "p", "Fv", "Fi", "Fo"};
  qh_init_A(fin, fout, stderr, argc, argv);
  int status = setjmp(qh errexit); 
  if (!status)
  {
    qh_option("voronoi  _bbound-last  _coplanar-keep", NULL, NULL);
    qh DELAUNAY = true;
    qh VORONOI = true;
    qh SCALElast = true;
    qh_checkflags(qh qhull_command, hidden_options);
    qh_initflags(qh qhull_command);
    bool ismalloc = false; // FIXME: ???
    int dim = 3;
    qh_init_B(points, num_points, dim+1, ismalloc);
    qh_qhull();
    qh_check_output();
    qh_produce_output();
    if (qh VERIFYoutput && !qh FORCEoutput && !qh STOPpoint && !qh STOPcone)
      qh_check_points();
    status = qh_ERRnone;
  }
  qh NOerrexit = true; // No more setjmp.

  fclose(fin);
  fclose(fout);

  // Clean up QHull.
#ifdef qh_NOmem
  qh_freeqhull(True);
#else
  qh_freeqhull(False);
  int curlong, totlong;
  qh_memfreeshort(&curlong, &totlong);
  if (curlong || totlong)
    fprintf(stderr, "qhull internal warning (main): did not free %d bytes of long memory(%d pieces)\n",
       totlong, curlong);
#endif

  // At this point, the output buffer contains the tessellation. We simply 
  // need to parse it into our own data.
  printf("%s\n", output);

  // Start assembling our tessellation.
  voronoi_tessellation_t* t = (voronoi_tessellation_t*)malloc(sizeof(voronoi_tessellation_t));
#if 0

  // Record the node coordinates.
  t->num_nodes = qh num_vertices;
  t->nodes = (double*)malloc(3*t->num_nodes*sizeof(double));
  {
    int n = 0;
    FORALLvertices
    {
      t->nodes[3*n]   = vertex->point[0];
      t->nodes[3*n+1] = vertex->point[1];
      t->nodes[3*n+2] = vertex->point[2];
      ++n;
    }
  }

  // Build the face->edge mapping and the edge->node mapping.
  int_ptr_unordered_map_t* edges_for_face = int_ptr_unordered_map_new();
  int_ptr_unordered_map_t* nodes_for_edge = int_ptr_unordered_map_new();
  FORALLfacets
  {
    facet->seen = false;
    int_slist_t* face_edges = int_slist_new();
    FOREACHridge_(facet->ridges)
    {
      int n1 = -1, n2 = -1;
      FOREACHvertex_(ridge->vertices)
      {
        if (n1 == -1)
          n1 = vertex->id;
        else
        {
          ASSERT(n2 == -1);
          n2 = vertex->id;
        }
      }
      int* node_pair = int_tuple_new(2);
      node_pair[0] = MIN(n1, n2);
      node_pair[1] = MAX(n1, n2);
      int_ptr_unordered_map_insert_with_v_dtor(nodes_for_edge, ridge->id, node_pair, DTOR(int_tuple_free));
      int_slist_t* face_edges = int_slist_new();
      int_slist_append(face_edges, ridge->id);
    }
    int_ptr_unordered_map_insert_with_v_dtor(edges_for_face, facet->id, face_edges, DTOR(int_slist_free));
  }
  t->num_edges = nodes_for_edge->size;
  t->edges = (voronoi_edge_t*)malloc(t->num_edges*sizeof(voronoi_edge_t));
  for (int e = 0; e < t->num_edges; ++e)
  {
    int* node_pair = *int_ptr_unordered_map_get(nodes_for_edge, e);
    t->edges[e].node1 = node_pair[0];
    t->edges[e].node2 = node_pair[1];
    if (t->edges[e].node2 == -1) // node2 is a "ghost"
    {
      // FIXME
//      t->edges[i].ray[0] = out.vedgelist[i].vnormal[0];
//      t->edges[i].ray[1] = out.vedgelist[i].vnormal[1];
//      t->edges[i].ray[2] = out.vedgelist[i].vnormal[2];
    }
  }
  int_ptr_unordered_map_free(nodes_for_edge);

  // Build the cell->face mapping and the face->cell mapping.
  int_ptr_unordered_map_t* faces_for_cell = int_ptr_unordered_map_new();
  int_ptr_unordered_map_t* cells_for_face = int_ptr_unordered_map_new();
  FORALLvertices
  {
    qh_order_vertexneighbors(vertex);

    bool infinity_seen = false;
    facetT *neighbor, **neighborp;
    int_slist_t* cell_faces = int_slist_new();
    FOREACHneighbor_(vertex)
    {
      if (!neighbor->seen)
      {
        int c1 = vertex->id, c2 = -1;
        if (!neighbor->upperdelaunay)
        {
//          c2 = FIXME;
        }
        int* face_cells = int_tuple_new(2);
        face_cells[0] = c1;
        face_cells[1] = c2;
        int_ptr_unordered_map_insert_with_v_dtor(cells_for_face, neighbor->id, face_cells, DTOR(int_tuple_free));
        int_slist_append(cell_faces, neighbor->id);
        neighbor->seen = true;
      }
    }
    int_ptr_unordered_map_insert_with_v_dtor(faces_for_cell, vertex->id, cell_faces, DTOR(int_slist_free));
  }

  // Face <-> edge/cell connectivity.
  t->num_faces = edges_for_face->size;
  t->faces = (voronoi_face_t*)malloc(t->num_faces*sizeof(voronoi_face_t));
  for (int f = 0; f < t->num_faces; ++f)
  {
    int* face_cells = *int_ptr_unordered_map_get(cells_for_face, f);
    t->faces[f].cell1 = face_cells[0];
    t->faces[f].cell2 = face_cells[1];
    int_slist_t* face_edges = *int_ptr_unordered_map_get(edges_for_face, f);
    int Ne = face_edges->size;
    t->faces[f].num_edges = Ne;
    t->faces[f].edges = (int*)malloc(sizeof(int)*Ne);
    int_slist_node_t* fe = face_edges->front;
    for (int e = 0; e < Ne; ++e)
    {
      t->faces[f].edges[e] = fe->value;
      fe = fe->next;
    }
  }
  int_ptr_unordered_map_free(edges_for_face);
  int_ptr_unordered_map_free(cells_for_face);

  // Cell -> face connectivity.
  t->num_cells = qh num_good; //qh num_vertices - qh_setsize(qh del_vertices);
  t->cells = (voronoi_cell_t*)malloc(t->num_cells*sizeof(voronoi_cell_t));
  for (int i = 0; i < t->num_cells; ++i)
  {
    int_slist_t* cell_faces = *int_ptr_unordered_map_get(faces_for_cell, i);
    t->cells[i].num_faces = cell_faces->size;
    int_slist_node_t* cf = cell_faces->front;
    for (int f = 0; f < cell_faces->size; ++f)
    {
      t->cells[i].faces[f] = cf->value;
      cf = cf->next;
    }
  }
  int_ptr_unordered_map_free(faces_for_cell);
#endif

  return t;
}

