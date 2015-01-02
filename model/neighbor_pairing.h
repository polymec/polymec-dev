// Copyright (c) 2012-2015, Jeffrey N. Johnson
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

#ifndef POLYMEC_NEIGHBOR_PAIRING_H
#define POLYMEC_NEIGHBOR_PAIRING_H

#include "core/polymec.h"
#include "core/exchanger.h"
#include "core/serializer.h"
#include "core/adj_graph.h"
#include "core/point_cloud.h"
#include "core/silo_file.h"

// A neighbor pairing is a set of index pairs associated with some spatial 
// discretization (usually a mesh-free method). Pairings can be constructed 
// for any set of indices representing a spatial domain, and the underlying 
// indices and weights are stored in compressed arrays for efficiency.
typedef struct 
{
  char* name;
  int num_pairs;
  int* pairs;
  real_t* weights;
  exchanger_t* ex; // Used to perform exchanges to fill all values for indices on a domain.
} neighbor_pairing_t;

// Creates a neighbor pairing object with the given pairs and weights.
// pairs is an array of length 2*num_pairs whose (2*i)th element is the first 
// index in the ith pair, and whose (2*i+1)th is the second such index.
// The weights array is of length num_pairs whose ith element is the weight 
// associated with the ith pair. Returns a newly-allocated neighbor pairing. 
// The arrays must be allocated using polymec_malloc, and are consumed by the 
// pairing. Likewise, the pairing's exchanger is consumed and is used for its 
// exchanges.
neighbor_pairing_t* neighbor_pairing_new(const char* name, int num_pairs,
                                         int* pairs, real_t* weights,
                                         exchanger_t* ex);

// Creates a neighbor pairing object with weights all equal to one. See the 
// other constructor for details.
neighbor_pairing_t* unweighted_neighbor_pairing_new(const char* name, int num_pairs, 
                                                    int* pairs, exchanger_t* ex);

// Destroys the given neighbor pairing object.
void neighbor_pairing_free(neighbor_pairing_t* pairing);

// Returns an internal pointer to the exchanger associated with this 
// neighbor pairing.
exchanger_t* neighbor_pairing_exchanger(neighbor_pairing_t* pairing);

// Performs a synchronous exchange of the values for this pairing for the 
// given data. This method has the same signature as exchanger_exchange().
void neighbor_pairing_exchange(neighbor_pairing_t* pairing, void* data, int stride, int tag, MPI_Datatype type);

// Begins a synchronous exchange of the values for this pairing for the 
// given data. This method has the same signature as exchanger_start_exchange().
int neighbor_pairing_start_exchange(neighbor_pairing_t* pairing, void* data, int stride, int tag, MPI_Datatype type);

// Concludes the asynchronous exchange corresponding to the given token for 
// this pairing. This method has the same signature as exchanger_finish_exchange().
void neighbor_pairing_finish_exchange(neighbor_pairing_t* pairing, int token);

// Returns the number of neighbor pairs contained within this pairing.
static inline int neighbor_pairing_num_pairs(neighbor_pairing_t* pairing)
{
  return pairing->num_pairs;
}

// Traverses the neighbor pairing, returning true if the traversal
// has more pairs remaining and false if it has completed. The pos pointer 
// must be set to zero to reset the traversal. The i, j and weight pointers 
// will store the next (i, j) index pair and its weight, respectively.
// If weight is NULL, the weight obviously can't be returned.
static inline bool neighbor_pairing_next(neighbor_pairing_t* pairing, int* pos,
                                         int* i, int* j, real_t* weight)
{
  int k = *pos;
  if(k >= pairing->num_pairs)
    return false;
  *i = pairing->pairs[2*k];
  *j = pairing->pairs[2*k+1];
  if (weight != NULL)
    *weight = (pairing->weights != NULL) ? pairing->weights[k] : 1.0;
  ++(*pos);
  return true;
}

// Returns a serializer object that can read/write neighbor pairings from/to byte arrays.
serializer_t* neighbor_pairing_serializer();

// This function creates an adjacency graph for the given point cloud with 
// the given neighbor pairing.
adj_graph_t* graph_from_point_cloud_and_neighbors(point_cloud_t* points, 
                                                  neighbor_pairing_t* neighbors);

// This function extends the silo_file type to allow it to write out a 
// neighbor_pairing object to an entry with the given name.
void silo_file_write_neighbor_pairing(silo_file_t* file,
                                      const char* neighbors_name,
                                      neighbor_pairing_t* neighbors);

// This function extends the silo_file type to allow it to read in and 
// return a newly-allocated neighbor_pairing object from the entry in the 
// file with the given name. The exchanger for the neighbor pairing is assigned
// to the given MPI communicator.
neighbor_pairing_t* silo_file_read_neighbor_pairing(silo_file_t* file,
                                                    const char* neighbors_name,
                                                    MPI_Comm comm);


// This is a shake-n-bake method for constructing pairs of neighboring 
// points using a neighbor relationship that relates points (xi, xj) that are 
// closer to each other than a specified distance Dij = max(R[i], R[j]), where 
// R[i] and R[j] are radii of influence associated with point i and j.
// This method finds all neighbors on processes within the MPI communicator 
// of the given point cloud, and creates indices for ghost points corresponding
// to points on remote processes. It also stores the number of created ghost 
// indices in *num_ghost_points.
neighbor_pairing_t* distance_based_neighbor_pairing_new(point_cloud_t* points,
                                                        real_t* R,
                                                        int* num_ghost_points);

#endif
