// Copyright (c) 2012-2013, Jeffrey N. Johnson
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

  // Indicates to threads whether the pool is in the process of shutting down.
  bool shutting_down;

  // Task queue and number of remaining items. 
  ptr_slist_t* queue;
  int work_remaining;

  // Mutexes for controlling access to the queue and the information 
  // about finished jobs.
  pthread_mutex_t queue_lock, finished_lock;

  // Condition variables that are used for signaling other threads 
  // regarding ready and finished work.
  pthread_cond_t queue_cond, finished_cond;

  // Attributes of created threads are stored in this thing.
  pthread_attr_t thread_attr;

  // Threads!
  pthread_t* threads;

  // An array used only to indicate when threads have started.
  bool* thread_started;
};

// It's a thread's life! This function represents the entire lifetime of 
// a worker thread.
static void* thread_life(void* context)
{
  thread_context_t* thread = context;
  thread_pool_t* pool = thread->pool;
  pool->thread_started[thread->id] = true;
  
  while (true)
  {
    // Acquire the queue lock.
    pthread_mutex_lock(&pool->queue_lock);

    // Wait on the condition variable for work to show up in the queue.
    while (!pool->shutting_down && (pool->work_remaining == 0))
      pthread_cond_wait(&pool->queue_cond, &pool->queue_lock);

    if (pool->shutting_down)
      break;

    if (pool->work_remaining == 0)
      pthread_mutex_unlock(&pool->queue_lock);
    else
    {
      // Grab the next item in the queue and process it.
      void (*dtor)(void*);
      thread_work_t* work = ptr_slist_pop(pool->queue, &dtor);

      // Decrement the work remaining while we have the lock.
      --pool->work_remaining;

      // Surrender the queue lock.
      pthread_mutex_unlock(&pool->queue_lock);

      // Do the work.
      work->do_work(work->context);

      // Delete the work unit.
      dtor(work);

      // Mark the job as completed.
      pthread_mutex_lock(&pool->finished_lock);

      // If there's no more work, tell the host thread that we're finished.
      if (pool->work_remaining == 0)
        pthread_cond_signal(&pool->finished_cond);
      pthread_mutex_unlock(&pool->finished_lock);
    }
  }

  pthread_mutex_unlock(&pool->queue_lock);
  free(thread); // Kill the thread's context.
  pthread_exit(NULL);
  return NULL;
}

thread_pool_t* thread_pool_new()
{
  return thread_pool_with_threads(polymec_num_cores());
}

thread_pool_t* thread_pool_with_threads(int num_threads)
{
  ASSERT(num_threads > 0);

  thread_pool_t* pool = malloc(sizeof(thread_pool_t));
  pool->num_threads = num_threads;
  pool->threads = malloc(sizeof(pthread_t) * num_threads);
  pool->queue = ptr_slist_new();
  pool->shutting_down = false;
  pool->work_remaining = 0;

  // Thread creation attribute(s).
  pthread_attr_init(&pool->thread_attr);
  pthread_attr_setdetachstate(&pool->thread_attr, PTHREAD_CREATE_JOINABLE);

  // Ready the locks.
  pthread_mutex_init(&pool->queue_lock, NULL);
  pthread_mutex_init(&pool->finished_lock, NULL);

  // Ready the condition variable to signal the threads that there is 
  // work available (or that it's time to leave).
  pthread_cond_init(&pool->queue_cond, NULL);

  // Now for the condition variable that signals when thread work is finished.
  pthread_cond_init(&pool->finished_cond, NULL);

  // Prevent newly-created threads from doing anything silly by acquiring
  // a lock on our queue.
  pthread_mutex_lock(&pool->queue_lock);

  // Gentlemen, start your engines!
  pool->thread_started = malloc(sizeof(bool) * num_threads);
  memset(pool->thread_started, 0, sizeof(bool) * num_threads);
  for (int i = 0; i < num_threads; ++i)
  {
    thread_context_t* thread = malloc(sizeof(thread_context_t));
    thread->id = i;
    thread->pool = pool;
    int stat = pthread_create(&pool->threads[i], &pool->thread_attr, thread_life, thread);
    if (stat != 0)
      polymec_error("thread_pool: could not create thread %d.", i);
  }

  // Unlock the queue and wait for the threads to start up.
  pthread_mutex_unlock(&pool->queue_lock);
  int num_started;
  do 
  {
    num_started = 0;
    for (int i = 0; i < num_threads; ++i)
      num_started += (pool->thread_started[i]) ? 1 : 0;
  }
  while (num_started < num_threads);
  free(pool->thread_started); // <-- no longer needed.

  return pool;
}

void thread_pool_free(thread_pool_t* pool)
{
  if (pool->work_remaining != 0)
    polymec_error("thread_pool_free() called within thread_pool_execute().");

  // Tell the threads we're shutting down.
  pthread_mutex_lock(&pool->queue_lock);
  pool->shutting_down = true;
  pthread_cond_broadcast(&pool->queue_cond);
  pthread_mutex_unlock(&pool->queue_lock);

  // Wait for all the threads to join.
  for (int i = 0; i < pool->num_threads; ++i)
    pthread_join(pool->threads[i], NULL);

  pthread_cond_destroy(&pool->finished_cond);
  pthread_cond_destroy(&pool->queue_cond);
  pthread_mutex_destroy(&pool->finished_lock);
  pthread_mutex_destroy(&pool->queue_lock);
  pthread_attr_destroy(&pool->thread_attr);
  free(pool->threads);
  ptr_slist_free(pool->queue);
  free(pool);
}

int thread_pool_num_threads(thread_pool_t* pool)
{
  return pool->num_threads;
}

void thread_pool_schedule(thread_pool_t* pool,
                          void* context,
                          void (*do_work)(void*))
{
  if (pool->work_remaining != 0)
    polymec_error("thread_pool_schedule() called within thread_pool_execute().");

  thread_work_t* work = malloc(sizeof(thread_work_t));
  work->context = context;
  work->do_work = do_work;
  pthread_mutex_lock(&pool->queue_lock);
  ptr_slist_append_with_dtor(pool->queue, work, free);
  pthread_mutex_unlock(&pool->queue_lock);
}

void thread_pool_execute(thread_pool_t* pool)
{
  pthread_mutex_lock(&pool->finished_lock);

  // Signal the threads that there's work in the pool.
  pthread_mutex_lock(&pool->queue_lock);
  pool->work_remaining = pool->queue->size;
  pthread_cond_broadcast(&pool->queue_cond);
  pthread_mutex_unlock(&pool->queue_lock);

  // Now wait for the threads to finish.
  pthread_cond_wait(&pool->finished_cond, &pool->finished_lock);
  pthread_mutex_unlock(&pool->finished_lock);
}

