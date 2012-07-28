#ifndef ARBI_LITE_MESH_H
#define ARBI_LITE_MESH_H

#include "core/arbi.h"

#ifdef __cplusplus
extern "C" {
#endif

// Mesh centerings.
typedef enum
{
  LITE_MESH_FACE,
  LITE_MESH_CELL
} lite_mesh_centering_t;

typedef struct lite_face_t lite_face_t;

// A "light" cell. This cell is stationary for its lifetime and is 
// bounded by a set of faces.
typedef struct
{
  // Faces bounding this cell. If this cell is a ghost cell, the faces 
  // do not bound the entire cell--only those that connect this cell 
  // with interior cells are included.
  lite_face_t** faces;
  // Number of faces attached to this cell.
  int num_faces;
  // Cell volume.
  double volume;
} lite_cell_t;

// A "light" face. This face is stationary for its lifetime and connects
// exactly two cells.
struct lite_face_t
{
  // The cells connected by this face. cell1 is always the face with 
  // the lower index in the mesh.
  lite_cell_t* cell1;
  lite_cell_t* cell2;
  // Face area.
  double area;
};

// This data type represents a "light" unstructured mesh, consisting of 
// stationary cells and the faces connecting them.
typedef struct 
{
  // Mesh cells, indexed from 0 to N-1.
  lite_cell_t** cells;
  // Total number of (locally-owned) cells in the mesh.
  int num_cells;
  // Total number of ghost cells in the mesh.
  int num_ghost_cells;
  // Mesh faces, indexed from 0 to M-1.
  lite_face_t** faces;
  // Total number of faces in the mesh.
  int num_faces;
} lite_mesh_t;

// Construct a new lite_mesh with the given number of cells, ghost cells, 
// and faces.
lite_mesh_t* lite_mesh_new(int num_cells, int num_ghost_cells, int num_faces);

// Destroys the given lite_mesh.
void lite_mesh_free(lite_mesh_t* mesh);

#ifdef __cplusplus
}
#endif

#endif

