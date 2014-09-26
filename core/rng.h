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

#ifndef POLYMEC_RNG_H
#define POLYMEC_RNG_H

#include "core/polymec.h"

// This type represents a random number generator with a specific min, max, 
// and underlying algorithm. Objects of this type are garbage-collected.
typedef struct rng_t rng_t;

// This virtual table must be implemented by any space-time function.
typedef struct 
{
  void (*set_seed)(void* context, uint32_t seed);
  uint32_t (*get)(void* context);
  real_t (*uniform)(void* context);
  real_t (*uniform_positive)(void* context);
  uint32_t (*uniform_int)(void* context, uint32_t n);
  void (*dtor)(void* context);
} rng_vtable;

// Constructs a random number generator from the given name, context, and 
// virtual table.
rng_t* rng_new(const char* name, void* context, 
               uint32_t min, uint32_t max, rng_vtable vtable);

// Returns the name of the random number generator.
const char* rng_name(rng_t* rng);

// Sets the seed used by the random number generator.
void rng_set_seed(rng_t* rng, uint32_t seed);

// Returns the minimum random number that can be generated by this generator.
uint32_t rng_min(rng_t* rng);

// Returns the maximum random number that can be generated by this generator.
uint32_t rng_max(rng_t* rng);

// Generates and returns a random integer within the range [min, max], where 
// these bounds are, respectively, the minimum and maximum numbers that can 
// be generated by the random number generator. All integers in this range 
// are equally likely.
uint32_t rng_get(rng_t* rng);

// Generates a random floating point number within the range [0, 1), with 
// uniform probability.
real_t rng_uniform(rng_t* rng);

// Generates a random floating point number within the range (0, 1), with 
// uniform probability.
real_t rng_uniform_positive(rng_t* rng);

// Generates an integer within the range [0, n-1], with 
// uniform probability.
uint32_t rng_uniform_int(rng_t* rng, uint32_t n);

// Creates the standard C random number generator implemented using rand().
// This is available on every platform.
rng_t* rand_rng_new();

// Creates the best vanilla random number generator available on this system.
rng_t* host_rng_new();

#endif

