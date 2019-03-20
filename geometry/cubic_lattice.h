// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_CUBIC_LATTICE_H
#define POLYMEC_CUBIC_LATTICE_H

#include "core/unordered_map.h"
#include "core/serializer.h"

/// \addtogroup geometry geometry
///@{

/// \class cubic_lattice
/// This class defines an indexing scheme for a cubic lattice.
/// \refcounted
typedef struct
{
  // Number of cells in x, y, and z.
  index_t nx, ny, nz;
} cubic_lattice_t;

/// Constructs a representation of a cubic lattice with the given numbers of
/// cells in x, y, and z.
/// \memberof cubic_lattice
cubic_lattice_t* cubic_lattice_new(index_t nx, index_t ny, index_t nz);

/// Returns the number of cells in the lattice.
/// \memberof cubic_lattice
static inline index_t cubic_lattice_num_cells(cubic_lattice_t* l)
{
  return l->nx * l->ny * l->nz;
}

/// Returns the number of x-faces in the lattice.
/// \memberof cubic_lattice
static inline index_t cubic_lattice_num_xfaces(cubic_lattice_t* l)
{
  return (l->nx+1) * l->ny * l->nz;
}

/// Returns the number of y-faces in the lattice.
/// \memberof cubic_lattice
static inline index_t cubic_lattice_num_yfaces(cubic_lattice_t* l)
{
  return l->nx * (l->ny+1) * l->nz;
}

/// Returns the number of z-faces in the lattice.
/// \memberof cubic_lattice
static inline index_t cubic_lattice_num_zfaces(cubic_lattice_t* l)
{
  return l->nx * l->ny * (l->nz+1);
}

/// Returns the total number of faces in the lattice.
/// \memberof cubic_lattice
static inline index_t cubic_lattice_num_faces(cubic_lattice_t* l)
{
  return cubic_lattice_num_xfaces(l) + cubic_lattice_num_yfaces(l) +
         cubic_lattice_num_zfaces(l);
}

/// Returns the number of x-edges in the lattice.
/// \memberof cubic_lattice
static inline index_t cubic_lattice_num_xedges(cubic_lattice_t* l)
{
  return l->nx * (l->ny+1) * (l->nz+1);
}

/// Returns the number of y-edges in the lattice.
/// \memberof cubic_lattice
static inline index_t cubic_lattice_num_yedges(cubic_lattice_t* l)
{
  return (l->nx+1) * l->ny * (l->nz+1);
}

/// Returns the number of z-edges in the lattice.
/// \memberof cubic_lattice
static inline index_t cubic_lattice_num_zedges(cubic_lattice_t* l)
{
  return (l->nx+1) * (l->ny+1) * l->nz;
}

/// Returns the total number of edges in the lattice.
/// \memberof cubic_lattice
static inline index_t cubic_lattice_num_edges(cubic_lattice_t* l)
{
  return cubic_lattice_num_xedges(l) + cubic_lattice_num_yedges(l) +
         cubic_lattice_num_zedges(l);
}

/// Returns the number of nodes in the lattice.
/// \memberof cubic_lattice
static inline index_t cubic_lattice_num_nodes(cubic_lattice_t* l)
{
  return (l->nx+1) * (l->ny+1) * (l->nz+1);
}

// Returns the index of the cell corresponding to (i, j, k).
/// \memberof cubic_lattice
static inline index_t cubic_lattice_cell(cubic_lattice_t* l, index_t i, index_t j, index_t k)
{
  return l->nx*l->ny*k + l->nx*j + i;
}

/// Computes the triple (i, j, k) corresponding to the cell with the given index.
/// \memberof cubic_lattice
static inline void cubic_lattice_get_cell_triple(cubic_lattice_t* l, index_t index, index_t* i, index_t* j, index_t* k)
{
  *k = index/(l->nx*l->ny);
  *j = (index - l->nx*l->ny*(*k))/l->nx;
  *i = index - l->nx*l->ny*(*k) - l->nx*(*j);
}

/// Returns the index of the x-face corresponding to (i-1/2, j, k).
/// \memberof cubic_lattice
static inline index_t cubic_lattice_xface(cubic_lattice_t* l, index_t i, index_t j, index_t k)
{
  return (l->nx+1)*(l->ny)*k + (l->nx+1)*j + i;
}

/// Returns the index of the y-face corresponding to (i, j-1/2, k).
/// \memberof cubic_lattice
static inline index_t cubic_lattice_yface(cubic_lattice_t* l, index_t i, index_t j, index_t k)
{
  return cubic_lattice_num_xfaces(l) +
         (l->nx)*(l->ny+1)*k + (l->nx)*j + i;
}

/// Returns the index of the z-face corresponding to (i, j, k-1/2).
/// \memberof cubic_lattice
static inline index_t cubic_lattice_zface(cubic_lattice_t* l, index_t i, index_t j, index_t k)
{
  return cubic_lattice_num_xfaces(l) + cubic_lattice_num_yfaces(l) +
         (l->nx)*(l->ny)*k + (l->nx)*j + i;
}

/// Returns the index of the x-edge corresponding to (i, j-1/2, k-1/2).
/// \memberof cubic_lattice
static inline index_t cubic_lattice_xedge(cubic_lattice_t* l, index_t i, index_t j, index_t k)
{
  return (l->nx)*(l->ny+1)*k + (l->nx)*j + i;
}

/// Returns the index of the y-edge corresponding to (i-1/2, j, k-1/2).
/// \memberof cubic_lattice
static inline index_t cubic_lattice_yedge(cubic_lattice_t* l, index_t i, index_t j, index_t k)
{
  return cubic_lattice_num_xedges(l) +
         (l->nx+1)*(l->ny)*k + (l->nx+1)*j + i;
}

/// Returns the index of the z-edge corresponding to (i-1/2, j-1/2, k).
/// \memberof cubic_lattice
static inline index_t cubic_lattice_zedge(cubic_lattice_t* l, index_t i, index_t j, index_t k)
{
  return cubic_lattice_num_xedges(l) + cubic_lattice_num_yedges(l) +
         (l->nx+1)*(l->ny+1)*k + (l->nx+1)*j + i;
}

/// Returns the index of the node corresponding to (i-1/2, j-1/2, k-1/2).
/// \memberof cubic_lattice
static inline index_t cubic_lattice_node(cubic_lattice_t* l, index_t i, index_t j, index_t k)
{
  return (l->nx+1)*(l->ny+1)*k + (l->nx+1)*j + i;
}

/// Returns a serializer for cubic lattice objects.
/// \memberof cubic_lattice
serializer_t* cubic_lattice_serializer(void);

///@}

#endif

