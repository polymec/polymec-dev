#ifndef POLYMEC_DELAUNAY_TRIANGULATOR_H
#define POLYMEC_DELAUNAY_TRIANGULATOR_H

#ifdef __cplusplus
extern "C" {
#endif

// This represents a Delaunay tetrahedron.
typedef struct
{
  int vertices[4];
} delaunay_tet_t;

// This type holds a Delaunay triangulation.
typedef struct
{
  int num_tets;
  delaunay_tet_t* tets;

  int num_vertices;
  double* vertices;

} delaunay_triangulation_t;

// Creates a new Delaunay triangulation with the given numbers of tetrahedra
// and vertices.
delaunay_triangulation_t* delaunay_triangulation_new(int num_tets, int num_vertices);

// Destroys a triangulation that has been created.
void delaunay_triangulation_free(delaunay_triangulation_t* triangulation);

// The tessellator class creates Delaunay triangulations from sets of 
// vertices. Objects of this type are garbage-collected.
typedef struct delaunay_triangulator_t delaunay_triangulator_t;

// Constructs a new Delaunay triangulator.
delaunay_triangulator_t* delaunay_triangulator_new();

// Creates a triangulation with the given set of vertices. 
// points is a (3*num_points) array containing the coordinates of the 
// vertices in point-major order.
delaunay_triangulation_t* 
delaunay_triangulator_triangulate(delaunay_triangulator_t* triangulator,
                                  double* points, int num_points);

#ifdef __cplusplus
}
#endif

#endif

