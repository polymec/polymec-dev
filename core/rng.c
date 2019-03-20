// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/rng.h"

struct rng_t
{
  char* name;
  void* context;
  rng_vtable vtable;
  uint32_t min, max;
  bool has_global_state, is_thread_safe;
};

static void rng_free(void* ctx)
{
  rng_t* rng = ctx;
  string_free(rng->name);
  if ((rng->vtable.dtor != NULL) && (rng->context != NULL))
    rng->vtable.dtor(rng->context);
}

rng_t* rng_new(const char* name, void* context,
               uint32_t min, uint32_t max, rng_vtable vtable,
               bool has_global_state, bool is_thread_safe)
{
  ASSERT(min < max);
  ASSERT(vtable.get != NULL);
  rng_t* rng = polymec_refcounted_malloc(sizeof(rng_t), rng_free);
  rng->name = string_dup(name);
  rng->context = context;
  rng->min = min;
  rng->max = max;
  rng->vtable = vtable;
  rng->has_global_state = has_global_state;
  rng->is_thread_safe = is_thread_safe;
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
  if (rng->vtable.uniform != NULL)
    return rng->vtable.uniform(rng->context);
  else
    return (real_t)(1.0 * rng_get(rng) / rng->max);
}

real_t rng_uniform_positive(rng_t* rng)
{
  if (rng->vtable.uniform_positive != NULL)
    return rng->vtable.uniform_positive(rng->context);
  else
  {
    while (true)
    {
      real_t val = rng_uniform(rng);
      if (val > 0.0)
        return val;
    }
  }
}

uint32_t rng_uniform_int(rng_t* rng, uint32_t n)
{
  if (rng->vtable.uniform_int != NULL)
    return rng->vtable.uniform_int(rng->context, n);
  else
  {
    uint32_t reduction_factor = rng->max / n;
    return rng_get(rng) / reduction_factor;
  }
}

#ifdef _BSD_SOURCE

// POSIX random() generator.

static const uint32_t posix_max = UINT_MAX;

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

static void posix_free(void* context)
{
  polymec_free(context);
}

rng_t* posix_rng_new(size_t state_size)
{
  static const int num_bytes = 256;
  char* state = polymec_malloc(sizeof(char) * num_bytes);
  rng_vtable vtable = {.set_seed = posix_set_seed,
                       .get = posix_get,
                       .dtor = posix_free};
  initstate(random(), state, num_bytes);
  return rng_new("posix RNG", state, 0, posix_max, vtable, true, false);
}
#endif

#if defined APPLE && APPLE

// ARC4 random generator.

static const uint32_t arc4_max = UINT_MAX;

static uint32_t arc4_get(void* context)
{
  return arc4random();
}

static uint32_t arc4_uniform_int(void* context, uint32_t n)
{
  return arc4random_uniform(n);
}

static rng_t* arc4_rng_new()
{
  rng_vtable vtable = {.get = arc4_get,
                       .uniform_int = arc4_uniform_int};
  return rng_new("arc4 RNG", NULL, 0, arc4_max, vtable, true, true);
}
#endif

// Standard C rand() generator -- always available.

static const uint32_t rand_max = RAND_MAX;

static void rand_set_seed(void* context, uint32_t seed)
{
  srand(seed);
}

static uint32_t rand_get(void* context)
{
  return (uint32_t)rand();
}

rng_t* rand_rng_new()
{
  rng_vtable vtable = {.set_seed = rand_set_seed,
                       .get = rand_get};
  return rng_new("rand (standard C) RNG", NULL, 0, rand_max, vtable,
                 true, false);
}

rng_t* host_rng_new()
{
#if defined APPLE && APPLE
  return arc4_rng_new();
#elif defined(_BSD_SOURCE)
  return posix_rng_new();
#else
  return rand_rng_new();
#endif
}

