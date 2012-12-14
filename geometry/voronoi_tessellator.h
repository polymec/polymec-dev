#ifndef POLYMEC_VORONOI_TESSELLATOR_H
#define POLYMEC_VORONOI_TESSELLATOR_H

#ifdef __cplusplus
extern "C" {
#endif

// This represents a Voronoi cell in a tessellation.
typedef struct
{
  int num_faces;
  int* faces;
} voronoi_cell_t;

// This represents a Voronoi face separating two cells.
typedef struct
{
  int num_edges;
  int* edges;
  int cell1, cell2;
} voronoi_face_t;

// This represents a Voronoi edge connecting two nodes.
typedef struct
{
  int node1, node2;
  double ray[3];
} voronoi_edge_t;

// This type holds a Voronoi tessellation produced by a tessellator.
typedef struct
{
  int num_cells;
  voronoi_cell_t* cells;

  int num_faces;
  voronoi_face_t* faces;

  int num_edges;
  voronoi_edge_t* edges;

  int num_nodes;
  double* nodes;

} voronoi_tessellation_t;

// Destroys a tessellation that has been created by a tessellator.
void voronoi_tessellation_free(voronoi_tessellation_t* tessellation);

// The tessellator class creates Voronoi tessellations from sets of 
// generator points. Objects of this type are garbage-collected.
typedef struct voronoi_tessellator_t voronoi_tessellator_t;

// Constructs a new Voronoi tessellator.
voronoi_tessellator_t* voronoi_tessellator_new();

// Creates a tessellation with the given set of generators. The tessellation 
// is unbounded, and its "outermost" cells have infinite extent. The edges 
// on these cells have node2 == -1, and their "ray" field is a vector 
// pointing outward to infinity. points is a (3*num_points) array containing 
// the coordinates of the generator points in point-major order.
voronoi_tessellation_t* 
voronoi_tessellator_tessellate(voronoi_tessellator_t* tessellator,
                               double* points, int num_points);

#ifdef __cplusplus
}
#endif

#endif

