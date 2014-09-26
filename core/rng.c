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

#include <gc/gc.h>
#include "core/rng.h"

struct rng_t 
{
  char* name;
  void* context;
  rng_vtable vtable;
  uint32_t min, max;
};

static void rng_free(void* ctx, void* dummy)
{
  rng_t* rng = ctx;
  string_free(rng->name);
  if ((rng->vtable.dtor != NULL) && (rng->context != NULL))
    rng->vtable.dtor(rng->context);
}

rng_t* rng_new(const char* name, void* context, 
               uint32_t min, uint32_t max, rng_vtable vtable)
{
  ASSERT(min < max);
  ASSERT(vtable.get != NULL);
  ASSERT(vtable.uniform != NULL);
  ASSERT(vtable.uniform_positive != NULL);
  ASSERT(vtable.uniform_int != NULL);
  rng_t* rng = GC_MALLOC(sizeof(rng_t));
  rng->name = string_dup(name);
  rng->context = context;
  rng->min = min;
  rng->max = max;
  rng->vtable = vtable;
  GC_register_finalizer(rng, rng_free, rng, NULL, NULL);
  return rng;
}

const char* rng_name(rng_t* rng)
{
  return (const char*)rng->name;
}

void rng_set_seed(rng_t* rng, uint32_t seed)
{
  if (rng->vtable.set_seed != NULL)
    rng->vtable.set_seed(rng->context, seed);
}

uint32_t rng_min(rng_t* rng)
{
  return rng->min;
}

uint32_t rng_max(rng_t* rng)
{
  return rng->max;
}

uint32_t rng_get(rng_t* rng)
{
  return rng->vtable.get(rng->context);
}

real_t rng_uniform(rng_t* rng)
{
  return rng->vtable.uniform(rng->context);
}

real_t rng_uniform_positive(rng_t* rng)
{
  return rng->vtable.uniform_positive(rng->context);
}

uint32_t rng_uniform_int(rng_t* rng, uint32_t n)
{
  return rng->vtable.uniform_int(rng->context, n);
}

// POSIX random() generator.

static const uint32_t posix_max = (uint32_t)(-1);

static void posix_set_seed(void* context, uint32_t seed)
{
  char* state = context;
  setstate((const char*)state);
  srandom(seed);
}

static uint32_t posix_get(void* context)
{
  return (uint32_t)random();
}

static real_t posix_uniform(void* context)
{
  return (real_t)(1.0 * random() / (posix_max + 1.0));
}

static real_t posix_uniform_positive(void* context)
{
  while (true)
  {
    real_t val = posix_uniform(context);
    if (val != 0.0)
      return val;
  }
}

static uint32_t posix_uniform_int(void* context, uint32_t n)
{
  uint32_t reduction_factor = posix_max / n;
  return posix_get(context) / reduction_factor;
}

static void posix_free(void* context)
{
  polymec_free(context);
}

rng_t* posix_rng_new()
{
  static const int num_bytes = 256;
  char* state = polymec_malloc(sizeof(char) * num_bytes);
  rng_vtable vtable = {.set_seed = posix_set_seed,
                       .get = posix_get,
                       .uniform = posix_uniform,
                       .uniform_positive = posix_uniform_positive,
                       .uniform_int = posix_uniform_int,
                       .dtor = posix_free};
  initstate(random(), state, num_bytes);
  return rng_new("posix RNG", state, 0, posix_max, vtable);
}

#if APPLE 

// ARC4 random generator.

static const uint32_t arc4_max = (uint32_t)(-1);

static uint32_t arc4_get(void* context)
{
  return arc4random();
}

static real_t arc4_uniform(void* context)
{
  return (real_t)(1.0 * arc4random() / (arc4_max + 1.0));
}

static real_t arc4_uniform_positive(void* context)
{
  while (true)
  {
    real_t val = arc4_uniform(context);
    if (val != 0.0)
      return val;
  }
}

static uint32_t arc4_uniform_int(void* context, uint32_t n)
{
  return arc4random_uniform(n);
}

static rng_t* arc4_rng_new()
{
  rng_vtable vtable = {.get = arc4_get,
                       .uniform = arc4_uniform,
                       .uniform_positive = arc4_uniform_positive,
                       .uniform_int = arc4_uniform_int};
  return rng_new("arc4 RNG", NULL, 0, arc4_max, vtable);
}
#endif

rng_t* host_rng_new()
{
#if APPLE 
  return arc4_rng_new();
#else
  return posix_rng_new();
#endif
}

