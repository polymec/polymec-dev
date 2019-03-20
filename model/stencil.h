// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_STENCIL_H
#define POLYMEC_STENCIL_H

#include "core/polymec.h"
#include "core/exchanger.h"
#include "core/serializer.h"
#include "core/unordered_set.h"
#include "core/adj_graph.h"
#include "geometry/point_cloud.h"
#include "solvers/matrix_sparsity.h"

/// \addtogroup model model
///@{

/// \class stencil
/// A stencil is a set of indices associated with a stencil for some spatial
/// discretization. Stencils can be constructed for any set of indices
/// representing a spatial domain, and the underlying indices are stored in
/// compressed arrays for efficiency.
typedef struct
{
  char* name;
  size_t num_indices;
  int* offsets;
  int* indices;
  size_t num_ghosts; // Number of remotely-owned indices.
  exchanger_t* ex; // Used to perform exchanges to fill all values for indices on a domain.
} stencil_t;

/// Creates a stencil object using compressed arrays for the offsets of
/// the different stencils for each index and the indices of each stencil
/// therein. Returns a newly-allocated stencil. These arrays must be allocated
/// using polymec_malloc, and are consumed by the stencil. Likewise, the
/// stencil's exchanger is consumed by the stencil and is used for its
/// exchanges.
/// \memberof stencil
stencil_t* stencil_new(const char* name, size_t num_indices,
                       int* offsets, int* indices,
                       size_t num_ghosts, exchanger_t* ex);

/// Destroys the given stencil object.
/// \memberof stencil
void stencil_free(stencil_t* stencil);

/// Creates a (deep) copy of the given stencil.
/// \memberof stencil
stencil_t* stencil_clone(stencil_t* stencil);

/// Given an array of sets of neighbors to remove from each index (of length
/// stencil->num_indices), removes those neighbors from their indices within
/// the stencil.
/// \param [in] neighbors_to_trim An array specifying sets of neighbors to trim.
///             The ith entry holds the set of neighbors to trim from index i
///             (or can be NULL if there are no such neighbors).
/// \memberof stencil
void stencil_trim(stencil_t* stencil, int_unordered_set_t** neighbors_to_trim);

/// Performs a synchronous exchange of the values for this stencil for the
/// given data. This method has the same signature as exchanger_exchange().
/// \memberof stencil
void stencil_exchange(stencil_t* stencil, void* data, int stride, int tag, MPI_Datatype type);

/// Begins an asynchronous exchange of the values for this stencil for the
/// given data. This method has the same signature as exchanger_start_exchange().
/// \memberof stencil
int stencil_start_exchange(stencil_t* stencil, void* data, int stride, int tag, MPI_Datatype type);

/// Concludes the asynchronous exchange corresponding to the given token for
/// this stencil. This method has the same signature as exchanger_finish_exchange().
/// \memberof stencil
void stencil_finish_exchange(stencil_t* stencil, int token);

/// Provides direct access to the stencil's \ref exchanger.
/// \memberof stencil
exchanger_t* stencil_exchanger(stencil_t* stencil);

/// Returns the number of indices for which the stencil has data.
/// \memberof stencil
static inline size_t stencil_num_indices(stencil_t* stencil)
{
  return stencil->num_indices;
}

/// Returns the number of indices in the stencil for the given index i.
/// \memberof stencil
static inline size_t stencil_size(stencil_t* stencil, int i)
{
  ASSERT(i < stencil->num_indices);
  return (size_t)(stencil->offsets[i+1] - stencil->offsets[i]);
}

/// Copies the indices of the neighbors for the index i into the given neighbors array.
/// \memberof stencil
static inline void stencil_get_neighbors(stencil_t* stencil, int i, int* neighbors)
{
  ASSERT(i < stencil->num_indices);
  int j1 = stencil->offsets[i], j2 = stencil->offsets[i+1];
  memcpy(neighbors, &stencil->indices[j1], sizeof(int) * (j2 - j1));
}

/// Traverses the stencil for a given index i, returning true if the traversal
/// has more indices remaining and false if it has completed.
/// \memberof stencil
/// \param [in] i The index for which the stencil indices are traversed.
/// \param [in,out] pos Controls the interation. Set to 0 to reset the iteration.
/// \param [out] j Stores the next stencil index.
static inline bool stencil_next(stencil_t* stencil, int i, int* pos, int* j)
{
  ASSERT(i < stencil->num_indices);
  if (*pos >= (stencil->offsets[i+1] - stencil->offsets[i]))
    return false;
  int k = stencil->offsets[i] + *pos;
  *j = stencil->indices[k];
  ASSERT(*j != -1);
  ++(*pos);
  return true;
}

/// Returns the number of remotely-owned indices in the stencil.
/// \memberof stencil
static inline size_t stencil_num_ghosts(stencil_t* stencil)
{
  return stencil->num_ghosts;
}

/// Returns an adjacency graph that represents the stencil. This graph "borrows"
/// the stencil's data and consumes minimal resources.
/// \memberof stencil
/// \collective Collective on the stencil's exchanger communicator.
adj_graph_t* stencil_as_graph(stencil_t* stencil);

/// Returns a serializer object that can read/write stencils from/to byte arrays.
/// \memberof stencil
serializer_t* stencil_serializer(void);

/// Returns a newly-created matrix sparsity constructed from this stencil.
/// \memberof stencil
matrix_sparsity_t* matrix_sparsity_from_stencil(stencil_t* stencil);

/// This pre-fab function creates a stencil for points in a cloud that have
/// neighbors within a radius given by R[i] for the ith point. num_ghost_points
/// will store the number of ghost points in the stencil.
/// \memberof stencil
stencil_t* distance_based_point_stencil_new(point_cloud_t* points,
                                            real_t* R);

///@}

#endif
