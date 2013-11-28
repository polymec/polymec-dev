// Copyright 2012-2013 Jeffrey Johnson.
// 
// This file is part of Polymec, and is licensed under the Apache License, 
// Version 2.0 (the "License"); you may not use this file except in 
// compliance with the License. You may may find the text of the license in 
// the LICENSE file at the top-level source directory, or obtain a copy of 
// it at
// 
// http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

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
