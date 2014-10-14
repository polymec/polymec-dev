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

#include "model/neighbor_pairing.h"

neighbor_pairing_t* neighbor_pairing_new(const char* name, int num_pairs, 
                                         int* pairs, real_t* weights,
                                         exchanger_t* ex)
{
  ASSERT(num_pairs > 0);
  ASSERT(pairs != NULL);
  neighbor_pairing_t* p = polymec_malloc(sizeof(neighbor_pairing_t));
  p->name = string_dup(name);
  p->num_pairs = num_pairs;
  p->pairs = pairs;
  p->weights = weights;
  p->ex = ex;
  return p;
}

neighbor_pairing_t* unweighted_neighbor_pairing_new(const char* name, int num_pairs, 
                                                    int* pairs, exchanger_t* ex)
{
  return neighbor_pairing_new(name, num_pairs, pairs, NULL, ex);
}

void neighbor_pairing_free(neighbor_pairing_t* pairing)
{
  polymec_free(pairing->name);
  polymec_free(pairing->pairs);
  if (pairing->weights != NULL)
    polymec_free(pairing->weights);
  exchanger_free(pairing->ex);
  polymec_free(pairing);
}

void neighbor_pairing_exchange(neighbor_pairing_t* pairing, void* data, int stride, int tag, MPI_Datatype type)
{
  exchanger_exchange(pairing->ex, data, stride, tag, type);
}

int neighbor_pairing_start_exchange(neighbor_pairing_t* pairing, void* data, int stride, int tag, MPI_Datatype type)
{
  return exchanger_start_exchange(pairing->ex, data, stride, tag, type);
}

void neighbor_pairing_finish_exchange(neighbor_pairing_t* pairing, int token)
{
  exchanger_finish_exchange(pairing->ex, token);
}

