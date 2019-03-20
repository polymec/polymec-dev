// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_NEIGHBOR_PAIRING_H
#define POLYMEC_NEIGHBOR_PAIRING_H

#include "model/stencil.h"

/// \addtogroup model model
///@{

/// \class neighbor_pairing
/// A neighbor pairing is a set of index pairs associated with some spatial
/// discretization (usually a mesh-free method). Pairings can be constructed
/// for any set of indices representing a spatial domain, and the underlying
/// indices are stored in compressed arrays for efficiency.
typedef struct
{
  char* name;
  size_t num_pairs;
  int* pairs;
  exchanger_t* ex; // Used to perform exchanges to fill all values for indices on a domain.
} neighbor_pairing_t;

/// Creates a neighbor pairing object with the given pairs.
/// pairs is an array of length 2*num_pairs whose (2*i)th element is the first
/// index in the ith pair, and whose (2*i+1)th is the second such index.
/// Returns a newly-allocated neighbor pairing. The array must be allocated
/// using polymec_malloc, and is consumed by the pairing. Likewise,
/// the pairing's exchanger is consumed and is used for its exchanges. An
/// exchanger must be given even for serial configurations.
/// \memberof neighbor_pairing
neighbor_pairing_t* neighbor_pairing_new(const char* name, size_t num_pairs,
                                         int* pairs, exchanger_t* ex);

/// Destroys the given neighbor pairing object.
/// \memberof neighbor_pairing
void neighbor_pairing_free(neighbor_pairing_t* pairing);

/// Returns an internal pointer to the exchanger associated with this
/// neighbor pairing.
/// \memberof neighbor_pairing
exchanger_t* neighbor_pairing_exchanger(neighbor_pairing_t* pairing);

/// Performs a synchronous exchange of the values for this pairing for the
/// given data. This method has the same signature as exchanger_exchange().
/// \memberof neighbor_pairing
void neighbor_pairing_exchange(neighbor_pairing_t* pairing, void* data, int stride, int tag, MPI_Datatype type);

// Begins a synchronous exchange of the values for this pairing for the
// given data. This method has the same signature as exchanger_start_exchange().
/// \memberof neighbor_pairing
int neighbor_pairing_start_exchange(neighbor_pairing_t* pairing, void* data, int stride, int tag, MPI_Datatype type);

/// Concludes the asynchronous exchange corresponding to the given token for
/// this pairing. This method has the same signature as exchanger_finish_exchange().
/// \memberof neighbor_pairing
void neighbor_pairing_finish_exchange(neighbor_pairing_t* pairing, int token);

/// Returns the number of neighbor pairs contained within this pairing.
/// \memberof neighbor_pairing
static inline size_t neighbor_pairing_num_pairs(neighbor_pairing_t* pairing)
{
  return pairing->num_pairs;
}

/// Retrieves the index pair (i, j) of the pair with the given index within
/// the neighbor pairing.
/// \memberof neighbor_pairing
static inline void neighbor_pairing_get(neighbor_pairing_t* pairing,
                                        int pair_index,
                                        int* i, int* j)
{
  ASSERT(pair_index >= 0);
  ASSERT(pair_index < pairing->num_pairs);

  *i = pairing->pairs[2*pair_index];
  *j = pairing->pairs[2*pair_index+1];
}

/// Traverses the neighbor pairing, returning true if the traversal
/// has more pairs remaining and false if it has completed. The pos pointer
/// must be set to zero to reset the traversal. The i, j pointers
/// store the next (i, j) index pair, respectively.
/// \memberof neighbor_pairing
static inline bool neighbor_pairing_next(neighbor_pairing_t* pairing, int* pos,
                                         int* i, int* j)
{
  int k = *pos;
  if(k >= (int)pairing->num_pairs)
    return false;
  neighbor_pairing_get(pairing, k, i, j);
  ++(*pos);
  return true;
}

/// Returns a \ref serializer that can read/write neighbor pairings from/to byte arrays.
/// \memberof neighbor_pairing
serializer_t* neighbor_pairing_serializer(void);

/// Creates a neighbor pairing from the given stencil, rendering a symmetric
/// representation of the neighbor relations in that stencil.
/// \memberof neighbor_pairing
neighbor_pairing_t* neighbor_pairing_from_stencil(stencil_t* stencil);

/// Creates a stencil from the given point cloud and neighbor pairing.
/// \memberof neighbor_pairing
/// \collective Collective on point_cloud's communicator.
stencil_t* stencil_from_point_cloud_and_neighbors(point_cloud_t* points,
                                                  neighbor_pairing_t* neighbors);

/// This function creates an adjacency graph for the given point cloud with
/// the given neighbor pairing.
/// \memberof neighbor_pairing
adj_graph_t* graph_from_point_cloud_and_neighbors(point_cloud_t* points,
                                                  neighbor_pairing_t* neighbors);

/// This function creates a matrix sparsity pattern for the given point cloud
/// with the given neighbor pairing.
/// \memberof neighbor_pairing
matrix_sparsity_t* sparsity_from_point_cloud_and_neighbors(point_cloud_t* points,
                                                           neighbor_pairing_t* neighbors);

/// This is a shake-n-bake method for constructing pairs of neighboring
/// points using a neighbor relationship that relates points (xi, xj) that are
/// closer to each other than a specified distance Dij = max(R[i], R[j]), where
/// R[i] and R[j] are radii of influence associated with point i and j.
/// This method finds all neighbors on processes within the MPI communicator
/// of the given point cloud, and creates indices for ghost points corresponding
/// to points on remote processes. It also stores the number of created ghost
/// indices in *num_ghost_points.
/// \memberof neighbor_pairing
/// \returns A newly constructed neighbor_pairing.
neighbor_pairing_t* distance_based_neighbor_pairing_new(point_cloud_t* points,
                                                        real_t* R,
                                                        int* num_ghost_points);

///@}

#endif
