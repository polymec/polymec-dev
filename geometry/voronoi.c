#include "geometry/voronoi.h"
#include "mpi.h"
#include "ctetgen.h"
#include "core/avl_tree.h"

#ifdef __cplusplus
extern "C" {
#endif

// This AVL node visitor appends the tree's data to a tag.
static void append_to_tag(avl_node_t* node, void* p)
{
  int* tag_p = (int*)p;
  *tag_p = (int)node->attribute;
  ++tag_p;
}

mesh_t* voronoi_tessellation(point_t* points, int num_points, 
                             point_t* ghost_points, int num_ghost_points)
{
  ASSERT(points != NULL);
  ASSERT(num_points >= 2);
  ASSERT(num_ghost_points >= 0);

  // Set Tetgen's options. We desire a Voronoi mesh.
  TetGenOpts opts;
  TetGenOptsInitialize(&opts);
  opts.voroout = 1;

  // Allocate an input PLC with all the points.
  PLC* in;
  PLCCreate(&in);
  in->numberofpoints = num_points + num_ghost_points;
  PetscMalloc(3*in->numberofpoints*sizeof(double), &in->pointlist);
  memcpy(in->pointlist, (double*)points, 3*num_points);
  memcpy(&in->pointlist[3*num_points], (double*)ghost_points, 3*num_ghost_points);

  // Allocate storage for the output PLC, which will store the 
  // tetrahedra and Voronoi polyhedra.
  PLC* out;
  PLCCreate(&out);

  // Tetrahedralize.
  TetGenTetrahedralize(&opts, in, out);
  ASSERT(out->numberofvcells == (num_points + num_ghost_points));

  // Construct the Voronoi graph.
  mesh_t* mesh = mesh_new(num_points - num_ghost_points, 
                          num_ghost_points,
                          out->numberofvfacets,
                          out->numberofvedges,
                          out->numberofvpoints);
  
  // Node coordinates.
  for (int i = 0; i < mesh->num_nodes; ++i)
  {
    mesh->nodes[i]->x = out->vpointlist[3*i];
    mesh->nodes[i]->y = out->vpointlist[3*i+1];
    mesh->nodes[i]->z = out->vpointlist[3*i+2];
  }

  // Edge <-> node connectivity.
  // NOTE: We keep track of "outer" edges and tag them accordingly.
  avl_tree_t* outer_edges = int_avl_tree_new();
  int num_outer_edges = 0;
  for (int i = 0; i < mesh->num_edges; ++i)
  {
    mesh->edges[i]->node1 = mesh->nodes[out->vedgelist[i].v1];
    int n2 = out->vedgelist[i].v2; // -1 if ghost
    if (n2 == -1)
    {
      avl_tree_insert(outer_edges, (void*)i);
      ++num_outer_edges;
      mesh->edges[i]->node2 = NULL;
    }
    else
    {
      mesh->edges[i]->node2 = mesh->nodes[n2];
    }
  }
  // Tag the outer edges as such.
  if (num_outer_edges > 0)
  {
    int* outer_edge_tag = mesh_create_tag(mesh->edge_tags, "outer", num_outer_edges);
    avl_node_t* root = avl_tree_root(outer_edges);
    int* tag_p = outer_edge_tag;
    avl_node_visit(root, &append_to_tag, (void*)tag_p);
  }

  // Face <-> edge/cell connectivity.
  for (int i = 0; i < mesh->num_faces; ++i)
  {
    mesh->faces[i]->cell1 = mesh->cells[out->vfacetlist[i].c1];
    mesh->faces[i]->cell2 = mesh->cells[out->vfacetlist[i].c2];
    int Ne = out->vfacetlist[i].elist[0];
    mesh->faces[i]->num_edges = Ne;
    mesh->faces[i]->edges = malloc(Ne*sizeof(edge_t*));
    for (int j = 0; j < Ne; ++j)
      mesh->faces[i]->edges[j] = mesh->edges[out->vfacetlist[i].elist[j+1]];
  }

  // Cell <-> face connectivity.
  for (int i = 0; i < mesh->num_cells; ++i)
  {
    int Nf = out->vcelllist[i][0];
    mesh->cells[i]->num_faces = Nf;
    mesh->cells[i]->faces = malloc(Nf*sizeof(face_t*));
    for (int j = 0; j < Nf; ++j)
      mesh->cells[i]->faces[j] = mesh->faces[out->vcelllist[i][j+1]];
  }

  // Clean up.
  PLCDestroy(&in);
  PLCDestroy(&out);

  return mesh;
}

void voronoi_intersect_with_boundary(mesh_t* mesh, sp_func_t* boundary)
{
}

#ifdef __cplusplus
}
#endif

