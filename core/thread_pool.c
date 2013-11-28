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

#include <pthread.h>
#include "core/thread_pool.h"
#include "core/slist.h"

typedef struct
{
  int id;
  thread_pool_t* pool;
} thread_context_t;

typedef struct 
{
  void* context;
  void (*do_work)(void*);
} thread_work_t;

struct thread_pool_t 
{
  int num_threads;
  ptr_slist_t* queue;
  pthread_mutex_t queue_lock;
  pthread_cond_t queue_cond;
  pthread_t* threads;
  bool* thread_started;
  pthread_attr_t thread_attr;
  bool shutting_down;
};

static void* thread_life(void* context)
{
  thread_context_t* thread = context;
  thread_pool_t* pool = thread->pool;
  
  while (true)
  {
    // Wait on the condition variable for work to show up in the queue.
    pthread_mutex_lock(&pool->queue_lock);
    pool->thread_started[thread->id] = true;
    pthread_cond_wait(&pool->queue_cond, &pool->queue_lock);

    // Is the pool shutting down? If so, we should exit.
    if (pool->shutting_down)
    {
      pthread_mutex_unlock(&pool->queue_lock);
      free(thread);
      pthread_exit(NULL);
    }

    while (!ptr_slist_empty(pool->queue))
    {
      // Grab the next item in the queue and process it.
      void (*dtor)(void*);
      thread_work_t* work = ptr_slist_pop(pool->queue, &dtor);
      pthread_mutex_unlock(&pool->queue_lock);

      // Do the work.
      work->do_work(work->context);

      // Delete the work unit.
      dtor(work);
    }
  }
}

thread_pool_t* thread_pool_new(int num_threads)
{
  ASSERT(num_threads > 0);

  thread_pool_t* pool = malloc(sizeof(thread_pool_t));
  pool->num_threads = num_threads;
  pool->threads = malloc(sizeof(pthread_t) * num_threads);
  pool->thread_started = malloc(sizeof(bool) * num_threads);
  memset(pool->thread_started, 0, sizeof(bool) * num_threads);
  pool->queue = ptr_slist_new();
  pool->shutting_down = false;

  // Thread creation attribute(s).
  pthread_attr_init(&pool->thread_attr);
  pthread_attr_setdetachstate(&pool->thread_attr, PTHREAD_CREATE_JOINABLE);

  // Ready the lock for the queue.
  pthread_mutex_init(&pool->queue_lock, NULL);

  // Ready the condition variable to signal the threads that there is 
  // work available (or that it's time to leave).
  pthread_cond_init(&pool->queue_cond, NULL);

  // Gentlemen, start your engines!
  for (int i = 0; i < num_threads; ++i)
  {
    thread_context_t* thread = malloc(sizeof(thread_context_t));
    thread->id = i;
    thread->pool = pool;
    int stat = pthread_create(&pool->threads[i], &pool->thread_attr, thread_life, thread);
    if (stat != 0)
      polymec_error("thread_pool: could not create thread %d.", i);
  }

  // Wait for the threads to start up.
  int num_started;
  do 
  {
    num_started = 0;
    for (int i = 0; i < num_threads; ++i)
      num_started += (pool->thread_started[i]) ? 1 : 0;
  }
  while (num_started < num_threads);

  return pool;
}

void thread_pool_free(thread_pool_t* pool)
{
  // Wait for the threads to join.
  pool->shutting_down = true;
  pthread_mutex_lock(&pool->queue_lock);
  pthread_cond_broadcast(&pool->queue_cond);
  pthread_mutex_unlock(&pool->queue_lock);
  for (int i = 0; i < pool->num_threads; ++i)
    pthread_join(pool->threads[i], NULL);

  pthread_cond_destroy(&pool->queue_cond);
  pthread_mutex_destroy(&pool->queue_lock);
  pthread_attr_destroy(&pool->thread_attr);
  free(pool->threads);
  ptr_slist_free(pool->queue);
  free(pool);
}

void thread_pool_schedule(thread_pool_t* pool,
                          void* context,
                          void (*do_work)(void*))
{
  thread_work_t* work = malloc(sizeof(thread_work_t));
  work->context = context;
  work->do_work = do_work;
  pthread_mutex_lock(&pool->queue_lock);
  ptr_slist_append_with_dtor(pool->queue, work, free);
  pthread_mutex_unlock(&pool->queue_lock);
}

void thread_pool_execute(thread_pool_t* pool)
{
  // Signal the threads that there's work in the pool.
  pthread_mutex_lock(&pool->queue_lock);
  pthread_cond_broadcast(&pool->queue_cond);
  pthread_mutex_unlock(&pool->queue_lock);
}

