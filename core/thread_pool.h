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

#ifndef POLYMEC_THREAD_POOL_H
#define POLYMEC_THREAD_POOL_H

// This type is intended to provide a pool of threads for performing a 
// given task.
typedef struct thread_pool_t thread_pool_t;

// Creates a thread pool with a number of threads equal to the number of 
// cores on the machine. Uses polymec_num_cores(). Upon return, the pool 
// is created and has started a number of threads.
thread_pool_t* thread_pool_new();

// Creates a thread pool with the given number of threads. Upon return, 
// the pool is created and has started a number of threads.
thread_pool_t* thread_pool_with_threads(int num_threads);

// Destroys a thread pool.
void thread_pool_free(thread_pool_t* pool);

// Returns the number of threads in the given thread pool.
int thread_pool_num_threads(thread_pool_t* pool);

// Schedules a task within the thread pool, providing a function to 
// perform the work and a context to be used by the function. Here, 
// do_work is a function that receives a context object that contains any 
// data to be used. The thread pool assumes no ownership over the context.
void thread_pool_schedule(thread_pool_t* pool,
                          void* context,
                          void (*do_work)(void*));

// Runs the thread pool, executing its queue of tasks, blocking until all
// are completed.
void thread_pool_execute(thread_pool_t* pool);

#endif
