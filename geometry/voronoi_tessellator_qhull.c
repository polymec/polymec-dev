// This implementation of the Voronoi tessellator uses QHull.

#include "libqhull/libqhull.h"
#include "libqhull/mem.h"
#include "libqhull/qset.h"

#include "core/polymec.h"
#include "core/unordered_map.h"
#include "core/tuple.h"
#include "core/slist.h"
#include "core/table.h"
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
  char *input = malloc(sizeof(char) * (num_points+1) * 256);
  int input_size = 0;
  {
    char first[128];
    snprintf(first, 128, "3\n%d\n", num_points);
    int len = strlen(first);
    memcpy(input, first, sizeof(char)*len);
    input_size += sizeof(char)*len;
  }
  for (int i = 0; i < num_points; ++i)
  {
    char next[256];
    snprintf(next, 256, "%g %g %g\n", points[3*i], points[3*i+1], points[3*i+2]);
    int len = strlen(next);
    memcpy(&input[input_size], next, sizeof(char)*len);
    input_size += sizeof(char)*len;
  }

  FILE* fin = fmemopen(input, input_size, "r");

  char* output;
  size_t output_size;
  FILE* fout = open_memstream(&output, &output_size);

  // Initialize QHull.
  int argc = 3;
  char* argv[] = {"qvoronoi", // The program itself
                  "p",        // Coordinates of nodes
//                  "Fn",       // List of neighbors for each node for building edges
                  "Fv"};      // List of ridges (edges) for faces / cell pairs
  qh_init_A(fin, fout, stderr, argc, argv);
  int status = setjmp(qh errexit); 

  // Run the thing.
  if (!status)
  {
    qh_option("voronoi  _bbound-last  _coplanar-keep", NULL, NULL);
    log_debug("voronoi_tessellator_qhull: Tessellating %d points\n", num_points);
    qh DELAUNAY = true;
    qh VORONOI = true;
    qh SCALElast = true;
    qh_checkflags(qh qhull_command, hidden_options);
    qh_initflags(qh qhull_command);
    int num_points, dim;
    boolT is_malloc;
    coordT* points = qh_readpoints(&num_points, &dim, &is_malloc);
    qh_init_B(points, num_points, dim, is_malloc);
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
  free(input);

  // At this point, the output buffer contains the tessellation. We simply 
  // need to parse it into our own data.

  // Count the number of lines in the output buffer.
  int num_lines = 0;
  for (size_t i = 0; i < output_size; ++i)
  {
    if (output[i] == '\n')
      ++num_lines;
  }
  int line_offsets[num_lines];
  line_offsets[0] = 0;
  {
    int i = 0, j = 0;
    while (i < output_size)
    {
      if (output[i] == '\n')
        line_offsets[++j] = i+1;
      ++i;
    }
  }

  // Parse the dimension in the file and make sure it's 3.
  int dim;
  sscanf(output, "%d\n", &dim);
  ASSERT(dim == 3);

  // Start assembling our tessellation.
  voronoi_tessellation_t* t = (voronoi_tessellation_t*)malloc(sizeof(voronoi_tessellation_t));

  // Nodes and node positions ("p" output).
  sscanf(&output[line_offsets[1]], "%d\n", &t->num_nodes);
  ASSERT(t->num_nodes > 0);
  t->nodes = malloc(3 * sizeof(double) * t->num_nodes);
  for (int i = 0; i < t->num_nodes; ++i)
  {
    sscanf(&output[line_offsets[2+i]], "%lg %lg %lg\n", 
           &(t->nodes[3*i]), &(t->nodes[3*i+1]), &(t->nodes[3*i+2]));
  }

  // Faces and edges ("Fv" output).
  // We do not construct faces that contain the point-at-infinity.
  sscanf(&output[line_offsets[2+t->num_nodes]], "%d\n", &t->num_faces);
  t->faces = malloc(sizeof(voronoi_face_t) * t->num_faces);
  memset(t->faces, 0, sizeof(voronoi_face_t) * t->num_faces);
  int_table_t* edge_for_nodes = int_table_new();
  int_ptr_unordered_map_t* faces_for_cell = int_ptr_unordered_map_new();
  t->num_edges = 0;
  int face_offset = 0;
  for (int f = 0; f < t->num_faces; ++f)
  {
    // Read the line and parse it into a tuple of face node indices.
    char face_nodes_str[8*100];
    int line_size = line_offsets[3+t->num_nodes+f+1] - line_offsets[3+t->num_nodes+f];
    memcpy(face_nodes_str, &output[line_offsets[3+t->num_nodes+f]], 
           sizeof(char) * line_size);
    char* str_p = &face_nodes_str[0];
    char* nn_str = strsep(&str_p, " ");
    int num_nodes = atoi(nn_str) - 2;
    char* g1_str = strsep(&str_p, " ");
    int g1 = atoi(g1_str);
    char* g2_str = strsep(&str_p, " ");
    int g2 = atoi(g2_str);

    int face_nodes[num_nodes];
    bool face_contains_point_at_infinity = false;
    for (int n = 0; n < num_nodes; ++n)
    {
      char* n_str;
      if (n < num_nodes - 1)
        n_str = strsep(&str_p, " ");
      else
        n_str = str_p;
      face_nodes[n] = atoi(n_str) - 1; // subtract 1 from output node index!
      if (face_nodes[n] < 0)
      {
        // This face contains the point-at-infinity, so it won't be created.
        face_contains_point_at_infinity = true;
        break;
      }
    }
    if (face_contains_point_at_infinity)
      continue;

    // Hook up the face to the (generator) cells.
    voronoi_face_t* face = &t->faces[face_offset];
    face->cell1 = g1;
    int_slist_t** cell1_faces = (int_slist_t**)int_ptr_unordered_map_get(faces_for_cell, g1);
    if (cell1_faces == NULL)
    {
      int_slist_t* faces = int_slist_new();
      int_ptr_unordered_map_insert_with_v_dtor(faces_for_cell, g1, faces, DTOR(int_slist_free));
      cell1_faces = (int_slist_t**)int_ptr_unordered_map_get(faces_for_cell, g1);
    }
    int_slist_append(*cell1_faces, face_offset);
//printf("cell %d sees face %d\n", g1, face_offset);
    face->cell2 = g2;
    int_slist_t** cell2_faces = (int_slist_t**)int_ptr_unordered_map_get(faces_for_cell, g2);
    if (cell2_faces == NULL)
    {
      int_slist_t* faces = int_slist_new();
      int_ptr_unordered_map_insert_with_v_dtor(faces_for_cell, g2, faces, DTOR(int_slist_free));
      cell2_faces = (int_slist_t**)int_ptr_unordered_map_get(faces_for_cell, g2);
    }
    int_slist_append(*cell2_faces, face_offset);
//printf("cell %d sees face %d\n", g2, face_offset);
    
    // Build edges out of each pair of nodes.
    face->edges = malloc(sizeof(int) * num_nodes);
    face->num_edges = 0;
    for (int n = 0; n < num_nodes; ++n)
    {
      int n1 = face_nodes[n];
      int n2 = face_nodes[(n+1)%num_nodes];

      // Insert this node pair as an edge if we don't already have it.
      int* edge_p = int_table_get(edge_for_nodes, MIN(n1, n2), MAX(n1, n2)), edge = -1;
      if (edge_p == NULL)
      {
        edge = t->num_edges;
        int_table_insert(edge_for_nodes, MIN(n1, n2), MAX(n1, n2), edge);
        ++(t->num_edges);
      }
      else 
        edge = *edge_p;

      // The face gets this edge if it doesn't already have it.
      bool has_edge = false;
      for (int nn = 0; nn < face->num_edges; ++nn)
      {
        if (face->edges[nn] == edge)
        {
          has_edge = true;
          break;
        }
      }
      if (!has_edge)
      {
        face->edges[face->num_edges] = edge;
        ++face->num_edges;
      }
    }
    ++face_offset;
  }
  t->num_faces = face_offset;

  // Construct the edges for the tessellation.
  t->edges = malloc(sizeof(voronoi_edge_t) * t->num_edges);
  int_table_cell_pos_t pos = int_table_start(edge_for_nodes);
  int n1, n2, e;
  while (int_table_next_cell(edge_for_nodes, &pos, &n1, &n2, &e))
  {
    if (n1 == -1)
    {
      t->edges[e].node1 = n2;
      t->edges[e].node2 = n1;
    }
    else
    {
      t->edges[e].node1 = n1;
      t->edges[e].node2 = n2;
    }
  }

  // Make sure each Voronoi cell has a face.
  for (int c = 0; c < num_points; ++c)
  {
    if (!int_ptr_unordered_map_contains(faces_for_cell, c))
      polymec_error("Voronoi tessellation failed: input site %d at (%g, %g, %g) has no faces!", c, points[3*c], points[3*c+1], points[3*c+2]);
  }

  // Construct the cells for the tessellation. We do not construct any cell 
  // that has fewer than 4 faces.
  int cell_offset = 0;
  t->num_cells = num_points; 
  t->cells = malloc(t->num_cells*sizeof(voronoi_cell_t));
  for (int c = 0; c < num_points; ++c)
  {
    int_slist_t* cell_faces = *int_ptr_unordered_map_get(faces_for_cell, c);
    if (cell_faces->size < 4)
      continue;
    t->cells[cell_offset].num_faces = cell_faces->size;
    t->cells[cell_offset].faces = malloc(sizeof(voronoi_face_t) * cell_faces->size);
    int_slist_node_t* f = cell_faces->front;
    int fi = 0;
    while (f != NULL)
    {
      t->cells[cell_offset].faces[fi++] = f->value;
      f = f->next;
    }
    ++cell_offset;
  }
  if (cell_offset == 0)
    polymec_error("Voronoi tessellation produced no cells!");

  t->num_cells = cell_offset;
  t->cells = realloc(t->cells, t->num_cells*sizeof(voronoi_cell_t));

  // Now calculate the outward rays using the bounding hyperplanes for 
  // the exterior cells ("Fo" output).
  // FIXME
 
  int_table_free(edge_for_nodes);
  int_ptr_unordered_map_free(faces_for_cell);
//  int_ptr_unordered_map_free(rays_for_edge);

  return t;
}

