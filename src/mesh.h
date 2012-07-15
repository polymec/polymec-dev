#ifndef ARBI_MESH_H
#define ARBI_MESH_H

#include "arbi.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct face_t face_t;

// A node in 1, 2, or 3D space.
typedef struct
{
  double x, y, z;
} node_t;

// An edge, which belongs to a face and which is bounded by two nodes.
typedef struct 
{
  // The face for this edge.
  face_t* face;
  // Nodes that bound this edge. node1 is always the node with the lower index.
  node_t* node1;
  node_t* node2;
  // Length of the edge.
  double length;
} edge_t;

// A full-featured cell. This cell is bounded by a set of faces.
typedef struct
{
  // Faces bounding this cell. If this cell is a ghost cell, the faces 
  // do not bound the entire cell--only those that connect this cell 
  // with interior cells are included.
  face_t** faces;
  // Number of faces attached to this cell.
  int num_faces;
  // Cell volume.
  double volume;
} cell_t;

// A full-featured face. This face connects exactly two cells and is 
// bounded by edges.
struct face_t
{
  // The cells connected by this face. cell1 is always the face with 
  // the lower index in the mesh.
  cell_t* cell1;
  cell_t* cell2;
  // Edges that bound this face.
  edge_t** edges;
  // The number of edges bounding this face.
  int num_edges;
  // Face area.
  double area;
};

// This data type represents a "light" unstructured mesh, consisting of 
// stationary cells and the faces connecting them.
typedef struct 
{
  // Mesh cells, indexed from 0 to C-1.
  cell_t** cells;
  // Total number of (locally-owned) cells in the mesh.
  int num_cells;
  // Total number of ghost cells in the mesh.
  int num_ghost_cells;
  // Mesh faces, indexed from 0 to F-1.
  face_t** faces;
  // Total number of faces in the mesh.
  int num_faces;
  // Mesh edges, indexed from 0 to E-1.
  edge_t** edges;
  // Total number of edges in the mesh.
  int num_edges;
  // Mesh nodes, indexed from 0 to N-1.
  node_t** nodes;
  // Total number of nodes in the mesh.
  int num_nodes;
} mesh_t;

// Construct a new mesh with the given number of cells, ghost cells, 
// faces, edges, and nodes. This function does not provide any description
// of the mesh's topology and is only useful in the construction of mesh 
// generation algorithms.
mesh_t* mesh_new(int num_cells, int num_ghost_cells, int num_faces,
                 int num_edges, int num_nodes);

// Destroys the given mesh.
void mesh_free(mesh_t* mesh);

#ifdef __cplusplus
}
#endif

#endif

