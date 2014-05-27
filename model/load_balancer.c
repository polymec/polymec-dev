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

#include "core/array.h"
#include "model/load_balancer.h"

struct load_balancer_t 
{
  void* context;
  load_balancer_vtable vtable;
  real_array_t* imbalance_history;
  int workload_stride;
};

load_balancer_t* load_balancer_new(void* context, 
                                   load_balancer_vtable vtable, 
                                   int workload_stride)
{
  ASSERT(workload_stride > 0);

  load_balancer_t* balancer = polymec_malloc(sizeof(load_balancer_t));
  balancer->context = context;
  balancer->vtable = vtable;
  balancer->imbalance_history = real_array_new();
  balancer->workload_stride = workload_stride;
  return balancer;
}

void load_balancer_free(load_balancer_t* balancer)
{
  if ((balancer->vtable.dtor != NULL) && (balancer->context != NULL))
    balancer->vtable.dtor(balancer->context);
  real_array_free(balancer->imbalance_history);
  polymec_free(balancer);
}

void* load_balancer_context(load_balancer_t* balancer)
{
  return balancer->context;
}

real_t load_balancer_imbalance(load_balancer_t* balancer)
{
  // Compute the workload for each process.
  real_t local_work[balancer->workload_stride];
  balancer->vtable.compute_workload(balancer->context, local_work);

  // Compute the load imbalance.
  return balancer->vtable.compute_imbalance(balancer->context, local_work);
}

bool load_balancer_rebalance(load_balancer_t* balancer, real_t imbalance)
{
  // Append the imbalance to our imbalance history.
  real_array_append(balancer->imbalance_history, imbalance);

  // Do we need to rebalance?
  bool needs_balancing = balancer->vtable.needs_balancing(balancer->context, balancer->imbalance_history->data, balancer->imbalance_history->size);

  // If so, do it and clear the history.
  if (needs_balancing)
  {
    balancer->vtable.rebalance(balancer->context);
    real_array_clear(balancer->imbalance_history);
  }

  return needs_balancing;
}

