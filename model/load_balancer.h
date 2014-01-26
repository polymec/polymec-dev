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

#ifndef POLYMEC_LOAD_BALANCER_H
#define POLYMEC_LOAD_BALANCER_H

#include "core/polymec.h"

// A load balancer measures the workload on each parallel domain, decides 
// when (and how) to rebalance this workload, and performs the balancing.
typedef struct load_balancer_t load_balancer_t;

// A function for computing the workload on a parallel domain.
// The workload, which can be a single number or an array of numbers, is 
// stored in workload.
typedef void (*load_balancer_compute_workload_func)(void* context, real_t* local_workload);

// A function for computing the load imbalance on all processors, given 
// the workload on the local process.
typedef real_t (*load_balancer_compute_imbalance_func)(void* context, real_t* local_workload);

// A function for deciding whether the workload needs to be rebalanced 
// given the set of maximum load imbalances recorded since the last 
// rebalancing. The imbalance_history array contains all the load imbalances 
// recorded since the last rebalancing.
typedef bool (*load_balancer_needs_balancing_func)(void* context, real_t* imbalance_history, int num_imbalances);

// A function for rebalancing the workloads.
typedef void (*load_balancer_rebalance_func)(void* context);

// A function for destroying the data context for the balancer.
typedef void (*load_balancer_dtor)(void* context);

// This virtual table must be implemented by any load balancer.
typedef struct 
{
  load_balancer_compute_workload_func   compute_workload;
  load_balancer_compute_imbalance_func  compute_imbalance;
  load_balancer_needs_balancing_func    needs_balancing;
  load_balancer_rebalance_func          rebalance;
  load_balancer_dtor                    dtor;
} load_balancer_vtable;

// Creates an instance of a load balancer with the given context, virtual 
// table, and stride (number of numbers in) a workload.
load_balancer_t* load_balancer_new(void* context, 
                                   load_balancer_vtable vtable, 
                                   int workload_stride);

// Destroys the load balancer.
void load_balancer_free(load_balancer_t* balancer);

// Returns the context object associated with the load balancer.
void* load_balancer_context(load_balancer_t* balancer);

// Computes the imbalance on the balancer's working set.
real_t load_balancer_imbalance(load_balancer_t* balancer);

// Given the current maximum imbalance, rebalances the workload if necessary, 
// returning true if the rebalance was performed and false if not. If no 
// rebalancing was done, the imbalance is recorded in the imbalance history.
bool load_balancer_rebalance(load_balancer_t* balancer, real_t imbalance);

#endif

