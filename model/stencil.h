// Copyright (c) 2012-2014, Jeffrey N. Johnson
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without 
// modification, are permitted provided that the following conditions are met:
// 
// 1. Redistributions of source code must retain the above copyright notice, this 
// list of conditions and the following disclaimer.
// 
// 2. Redistributions in binary form must reproduce the above copyright notice, 
// this list of conditions and the following disclaimer in the documentation 
// and/or other materials provided with the distribution.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE 
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#ifndef POLYMEC_STENCIL_H
#define POLYMEC_STENCIL_H

#include "core/polymec.h"
#include "core/mesh.h"
#include "core/point_cloud.h"

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
                       exchanger_t* ex);

// Creates a stencil object with weights all equal to one. See the other 
// constructor for details.
stencil_t* unweighted_stencil_new(const char* name, int num_indices, 
                                  int* offsets, int* indices, exchanger_t* ex);

// Destroys the given stencil object.
void stencil_free(stencil_t* stencil);

// Performs a synchronous exchange of the values for this stencil for the 
// given data. This method has the same signature as exchanger_exchange().
void stencil_exchange(stencil_t* stencil, void* data, int stride, int tag, MPI_Datatype type);

// Begins a synchronous exchange of the values for this stencil for the 
// given data. This method has the same signature as exchanger_start_exchange().
int stencil_start_exchange(stencil_t* stencil, void* data, int stride, int tag, MPI_Datatype type);

// Concludes the asynchronous exchange corresponding to the given token for 
// this stencil. This method has the same signature as exchanger_finish_exchange().
void stencil_finish_exchange(stencil_t* stencil, int token);

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

// Traverses the stencil for a given index i, returning true if the traversal
// has more indices remaining and false if it has completed. The pos pointer 
// must be set to zero to reset the traversal. The j and weight pointers will 
// store the next stencil index and weight, respectively.
static inline bool stencil_next(stencil_t* stencil, int i, int* pos,
                                int* j, real_t* weight)
{
  ASSERT(i < stencil->num_indices);
  if (*pos >= (stencil->offsets[i+1] - stencil->offsets[i]))
    return false;
  int k = stencil->offsets[i] + *pos;
  *j = stencil->indices[k];
  ASSERT(*j != -1);
  *weight = (stencil->weights != NULL) ? stencil->weights[k] : 1.0;
  ++(*pos);
  return true;
}

// Creates a stencil for the cells that share at least one face with a given 
// cell. The stencil is constructed for every cell in the given mesh.
stencil_t* cell_face_stencil_new(mesh_t* mesh);

// Creates a stencil for the cells that share at least one edge with a given 
// cell. The stencil is constructed for every cell in the given mesh.
stencil_t* cell_edge_stencil_new(mesh_t* mesh);

// Creates a stencil for the cells that share at least one node with a given 
// cell. The stencil is constructed for every cell in the given mesh.
stencil_t* cell_node_stencil_new(mesh_t* mesh);

#endif
