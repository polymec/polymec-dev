// Copyright 2012-2013 Jeffrey Johnson.
// 
// This file is part of Polymec, and is licensed under the Apache License, 
// Version 2.0 (the "License"); you may not use this file except in 
// compliance with the License. You may may find the text of the license in 
// the LICENSE file at the top-level source directory, or obtain a copy of 
// it at
// 
// http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef POLYMEC_CUBIC_LATTICE_H
#define POLYMEC_CUBIC_LATTICE_H

#include "core/unordered_map.h"
#include "core/periodic_bc.h"

// This class defines an indexing scheme for a cubic lattice.
// Objects of this type are garbage-collected.
typedef struct 
{
  // Number of cells in x, y, and z.
  int nx, ny, nz;
} cubic_lattice_t;

// Constructs a representation of a cubic lattice with the given numbers of 
// cells in x, y, and z.
cubic_lattice_t* cubic_lattice_new(int nx, int ny, int nz);

// Returns the number of cells in the lattice.
static inline int cubic_lattice_num_cells(cubic_lattice_t* l)
{
  return l->nx * l->ny * l->nz;
}

// Returns the number of x-faces in the lattice.
static inline int cubic_lattice_num_x_faces(cubic_lattice_t* l)
{
  return (l->nx+1) * l->ny * l->nz;
}

// Returns the number of y-faces in the lattice.
static inline int cubic_lattice_num_y_faces(cubic_lattice_t* l)
{
  return l->nx * (l->ny+1) * l->nz;
}

// Returns the number of z-faces in the lattice.
static inline int cubic_lattice_num_z_faces(cubic_lattice_t* l)
{
  return l->nx * l->ny * (l->nz+1);
}

// Returns the total number of faces in the lattice.
static inline int cubic_lattice_num_faces(cubic_lattice_t* l)
{
  return cubic_lattice_num_x_faces(l) + cubic_lattice_num_y_faces(l) + 
         cubic_lattice_num_z_faces(l);
}

// Returns the number of x-edges in the lattice.
static inline int cubic_lattice_num_x_edges(cubic_lattice_t* l)
{
  return l->nx * (l->ny+1) * (l->nz+1);
}

// Returns the number of y-edges in the lattice.
static inline int cubic_lattice_num_y_edges(cubic_lattice_t* l)
{
  return (l->nx+1) * l->ny * (l->nz+1);
}

// Returns the number of z-edges in the lattice.
static inline int cubic_lattice_num_z_edges(cubic_lattice_t* l)
{
  return (l->nx+1) * (l->ny+1) * l->nz;
}

// Returns the total number of edges in the lattice.
static inline int cubic_lattice_num_edges(cubic_lattice_t* l)
{
  return cubic_lattice_num_x_edges(l) + cubic_lattice_num_y_edges(l) + 
         cubic_lattice_num_z_edges(l);
}

// Returns the number of nodes in the lattice.
static inline int cubic_lattice_num_nodes(cubic_lattice_t* l)
{
  return (l->nx+1) * (l->ny+1) * (l->nz+1);
}

// Returns the index of the cell corresponding to (i, j, k).
static inline int cubic_lattice_cell(cubic_lattice_t* l, int i, int j, int k)
{
  return l->nx*l->ny*k + l->nx*j + i;
}

// Returns the index of the x-face corresponding to (i-1/2, j, k).
static inline int cubic_lattice_x_face(cubic_lattice_t* l, int i, int j, int k)
{
  return (l->nx+1)*(l->ny)*k + (l->nx+1)*j + i;
}

// Returns the index of the y-face corresponding to (i, j-1/2, k).
static inline int cubic_lattice_y_face(cubic_lattice_t* l, int i, int j, int k)
{
  return cubic_lattice_num_x_faces(l) + 
         (l->nx)*(l->ny+1)*k + (l->nx)*j + i;
}

// Returns the index of the z-face corresponding to (i, j, k-1/2).
static inline int cubic_lattice_z_face(cubic_lattice_t* l, int i, int j, int k)
{
  return cubic_lattice_num_x_faces(l) + cubic_lattice_num_y_faces(l) + 
         (l->nx)*(l->ny)*k + (l->nx)*j + i;
}

// Returns the index of the x-edge corresponding to (i, j-1/2, k-1/2).
static inline int cubic_lattice_x_edge(cubic_lattice_t* l, int i, int j, int k)
{
  return (l->nx)*(l->ny+1)*k + (l->nx)*j + i;
}

// Returns the index of the y-edge corresponding to (i-1/2, j, k-1/2).
static inline int cubic_lattice_y_edge(cubic_lattice_t* l, int i, int j, int k)
{
  return cubic_lattice_num_x_edges(l) + 
         (l->nx+1)*(l->ny)*k + (l->nx+1)*j + i;
}

// Returns the index of the z-edge corresponding to (i-1/2, j-1/2, k).
static inline int cubic_lattice_z_edge(cubic_lattice_t* l, int i, int j, int k)
{
  return cubic_lattice_num_x_edges(l) + cubic_lattice_num_y_edges(l) + 
         (l->nx+1)*(l->ny+1)*k + (l->nx+1)*j + i;
}

// Returns the index of the node corresponding to (i-1/2, j-1/2, k-1/2).
static inline int cubic_lattice_node(cubic_lattice_t* l, int i, int j, int k)
{
  return (l->nx+1)*(l->ny+1)*k + (l->nx+1)*j + i;
}

// These functions can be used to generate mappings between faces on cubic 
// lattice meshes with periodic boundary conditions. 

// Creates a periodic boundary condition, identifying x-faces with the given tags.
periodic_bc_t* cubic_lattice_x_periodic_bc_new(const char* tag1, const char* tag2);

// Creates a periodic boundary condition, identifying y-faces with the given tags.
periodic_bc_t* cubic_lattice_y_periodic_bc_new(const char* tag1, const char* tag2);

// Creates a periodic boundary condition, identifying z-faces with the given tags.
periodic_bc_t* cubic_lattice_z_periodic_bc_new(const char* tag1, const char* tag2);

#endif

