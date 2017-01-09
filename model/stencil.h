// Copyright (c) 2012-2017, Jeffrey N. Johnson
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
#include "core/point_cloud.h"
#include "core/adj_graph.h"
#include "core/matrix_sparsity.h"
#include "core/silo_file.h"

// A stencil is a set of indices and weights associated with a stencil for 
// some spatial discretization. Stencils can be constructed for any set of 
// indices representing a spatial domain, and the underlying indices and 
// weights are stored in compressed arrays for efficiency.
typedef struct 
{
  char* name;
  int num_indices;
  int* offsets;
  int* indices;
  real_t* weights;
  int num_ghosts; // Number of remotely-owned indices.
  exchanger_t* ex; // Used to perform exchanges to fill all values for indices on a domain.
} stencil_t;

// Creates a stencil object using compressed arrays for (1) the offsets of 
// the different stencils for each index, (2) the indices of each stencil 
// therein, and (3) the weights associated with the indices in each stencil.
// Returns a newly-allocated stencil. These arrays must be allocated using 
// polymec_malloc, and are consumed by the stencil. Likewise, the stencil's
// exchanger is consumed by the stencil and is used for its exchanges.
stencil_t* stencil_new(const char* name, int num_indices, 
                       int* offsets, int* indices, real_t* weights,
                       int num_ghosts, exchanger_t* ex);

// Creates a stencil object with weights all equal to one. See the other 
// constructor for details.
stencil_t* unweighted_stencil_new(const char* name, int num_indices, 
                                  int* offsets, int* indices, 
                                  int num_ghosts, exchanger_t* ex);

// Destroys the given stencil object.
void stencil_free(stencil_t* stencil);

// Creates a (deep) copy of the given stencil.
stencil_t* stencil_clone(stencil_t* stencil);

// This function sets the weights for the given stencil after its construction. 
// Any existing weights are deleted, and the weights array is consumed.
void stencil_set_weights(stencil_t* stencil, real_t* weights);

// Given an array of sets of neighbors to remove from each index (of length 
// stencil->num_indices), removes those neighbors from their indices within 
// the stencil. neighbors_to_trim[i] holds the set of neighbors to trim from 
// index i (or can be NULL if there are no such neighbors).
void stencil_trim(stencil_t* stencil, int_unordered_set_t** neighbors_to_trim);

// Performs a synchronous exchange of the values for this stencil for the 
// given data. This method has the same signature as exchanger_exchange().
void stencil_exchange(stencil_t* stencil, void* data, int stride, int tag, MPI_Datatype type);

// Begins an asynchronous exchange of the values for this stencil for the 
// given data. This method has the same signature as exchanger_start_exchange().
int stencil_start_exchange(stencil_t* stencil, void* data, int stride, int tag, MPI_Datatype type);

// Concludes the asynchronous exchange corresponding to the given token for 
// this stencil. This method has the same signature as exchanger_finish_exchange().
void stencil_finish_exchange(stencil_t* stencil, int token);

// Provides direct access to the stencil's exchanger.
exchanger_t* stencil_exchanger(stencil_t* stencil);

// Returns the number of indices for which the stencil has data.
static inline int stencil_num_indices(stencil_t* stencil)
{
  return stencil->num_indices;
}

// Returns the number of indices in the stencil for the given index i.
static inline int stencil_size(stencil_t* stencil, int i)
{
  ASSERT(i < stencil->num_indices);
  return stencil->offsets[i+1] - stencil->offsets[i];
}

// Copies the indices of the neighbors for the index i into the given neighbors array.
static inline void stencil_get_neighbors(stencil_t* stencil, int i, int* neighbors)
{
  ASSERT(i < stencil->num_indices);
  int j1 = stencil->offsets[i], j2 = stencil->offsets[i+1];
  memcpy(neighbors, &stencil->indices[j1], sizeof(int) * (j2 - j1));
}

// Traverses the stencil for a given index i, returning true if the traversal
// has more indices remaining and false if it has completed. The pos pointer 
// must be set to zero to reset the traversal. The j and weight pointers will 
// store the next stencil index and weight, respectively. If weight is NULL, 
// the weight obviously can't be returned.
static inline bool stencil_next(stencil_t* stencil, int i, int* pos,
                                int* j, real_t* weight)
{
  ASSERT(i < stencil->num_indices);
  if (*pos >= (stencil->offsets[i+1] - stencil->offsets[i]))
    return false;
  int k = stencil->offsets[i] + *pos;
  *j = stencil->indices[k];
  ASSERT(*j != -1);
  if (weight != NULL)
    *weight = (stencil->weights != NULL) ? stencil->weights[k] : 1.0;
  ++(*pos);
  return true;
}

// Returns the number of remotely-owned indices in the stencil.
static inline int stencil_num_ghosts(stencil_t* stencil)
{
  return stencil->num_ghosts;
}

// Returns an adjacency graph that represents the stencil. This graph "borrows" 
// the stencil's data and consumes minimal resources.
adj_graph_t* stencil_as_graph(stencil_t* stencil);

// Returns a serializer object that can read/write stencils from/to byte arrays.
serializer_t* stencil_serializer(void);

// Returns a newly-created matrix sparsity constructed from this stencil.
matrix_sparsity_t* sparsity_from_stencil(stencil_t* stencil);

// This function extends the silo_file type to allow it to write out a 
// stencil object to an entry with the given name.
void silo_file_write_stencil(silo_file_t* file,
                             const char* stencil_name,
                             stencil_t* stencil);

// This function extends the silo_file type to allow it to read in and 
// return a newly-allocated stencil object from the entry in the 
// file with the given name. The exchanger for the stencil is assigned
// to the given MPI communicator.
stencil_t* silo_file_read_stencil(silo_file_t* file,
                                  const char* stencil_name,
                                  MPI_Comm comm);

// This pre-fab function creates a stencil for points in a cloud that have 
// neighbors within a radius given by R[i] for the ith point. num_ghost_points
// will store the number of ghost points in the stencil. No weights are 
// assigned.
stencil_t* distance_based_point_stencil_new(point_cloud_t* points,
                                            real_t* R);

#endif
