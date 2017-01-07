// Copyright (c) 2012-2017, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// At the moment, we're between POSIX threads and C11 threads.
#if !defined(__STDC_NO_THREADS__) || __STDC_NO_THREADS__
#define USE_PTHREADS 1
#else
#define USE_PTHREADS 0
#endif

#if USE_PTHREADS
#include <pthread.h>
#else
#include <threads.h>
#endif

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
  int work_in_progress;

#if USE_PTHREADS
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
#else
  // Mutexes for controlling access to the queue and the information 
  // about finished jobs.
  mtx_t queue_lock, finished_lock;

  // Condition variables that are used for signaling other threads 
  // regarding ready and finished work.
  cnd_t queue_cond, finished_cond;

  // Threads!
  thrd_t* threads;
#endif

  // An array used only to indicate when threads have started.
  bool* thread_started;
};

static inline void lock_queue(thread_pool_t* pool)
{
#if USE_PTHREADS
    pthread_mutex_lock(&pool->queue_lock);
#else
    mtx_lock(&pool->queue_lock);
#endif
}

static inline void unlock_queue(thread_pool_t* pool)
{
#if USE_PTHREADS
    pthread_mutex_unlock(&pool->queue_lock);
#else
    mtx_unlock(&pool->queue_lock);
#endif
}

// It's a thread's life! This function represents the entire lifetime of 
// a worker thread.
#if USE_PTHREADS
static void* thread_life(void* context)
#else
static int thread_life(void* context)
#endif
{
  thread_context_t* thread = context;
  thread_pool_t* pool = thread->pool;
  pool->thread_started[thread->id] = true;
  
  while (true)
  {
    // Acquire the queue lock.
    lock_queue(pool);

    // Wait on the condition variable for work to show up in the queue.
#if USE_PTHREADS
    while (!pool->shutting_down && (pool->work_remaining == 0))
      pthread_cond_wait(&pool->queue_cond, &pool->queue_lock);
#else
    while (!pool->shutting_down && (pool->work_remaining == 0))
      cnd_wait(&pool->queue_cond, &pool->queue_lock);
#endif

    // If the thread pool was told to shut down, we leave quietly.
    if (pool->shutting_down)
      break;

    if (pool->work_remaining == 0)
      unlock_queue(pool);
    else
    {
      // Grab the next item in the queue.
      void (*dtor)(void*);
      thread_work_t* work = ptr_slist_pop(pool->queue, &dtor);

      // Decrement the work remaining in the queue, and increment the
      // work in progress.
      --pool->work_remaining;
      ++pool->work_in_progress;

      // Surrender the queue lock.
      unlock_queue(pool);

      // Do the work.
      work->do_work(work->context);

      // Delete the work unit.
      dtor(work);

      // Get the queue lock and decrement the remaining work.
      lock_queue(pool);
      --pool->work_in_progress;
      unlock_queue(pool);

      // Mark the job as completed.
#if USE_PTHREADS
      pthread_mutex_lock(&pool->finished_lock);
#else
      mtx_lock(&pool->finished_lock);
#endif

      // If there's no more work, tell the host thread that we're finished.
#if USE_PTHREADS
      if ((pool->work_remaining == 0) && (pool->work_in_progress == 0))
        pthread_cond_signal(&pool->finished_cond);
      pthread_mutex_unlock(&pool->finished_lock);
#else
      if ((pool->work_remaining == 0) && (pool->work_in_progress == 0))
        cnd_signal(&pool->finished_cond);
      mtx_unlock(&pool->finished_lock);
#endif
    }
  }

  unlock_queue(pool);
  polymec_free(thread); // Kill the thread's context.

#if USE_PTHREADS
  pthread_exit(NULL);
  return NULL;
#else
  thrd_exit(0);
  return 0;
#endif
}

thread_pool_t* thread_pool_new()
{
  return thread_pool_with_threads(polymec_num_cores());
}

thread_pool_t* thread_pool_with_threads(int num_threads)
{
  ASSERT(num_threads > 0);

  thread_pool_t* pool = polymec_malloc(sizeof(thread_pool_t));
  pool->num_threads = num_threads;
  pool->threads = polymec_malloc(sizeof(pthread_t) * num_threads);
  pool->queue = ptr_slist_new();
  pool->shutting_down = false;
  pool->work_remaining = 0;
  pool->work_in_progress = 0;

#if USE_PTHREADS
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
#else
  // Ready the locks.
  mtx_init(&pool->queue_lock, mtx_plain);
  mtx_init(&pool->finished_lock, mtx_plain);

  // Ready the condition variable to signal the threads that there is 
  // work available (or that it's time to leave).
  cnd_init(&pool->queue_cond);

  // Now for the condition variable that signals when thread work is finished.
  cnd_init(&pool->finished_cond);

  // Prevent newly-created threads from doing anything silly by acquiring
  // a lock on our queue.
  mtx_lock(&pool->queue_lock);
#endif

  // Gentlemen, start your engines!
  pool->thread_started = polymec_malloc(sizeof(bool) * num_threads);
  memset(pool->thread_started, 0, sizeof(bool) * num_threads);
  for (int i = 0; i < num_threads; ++i)
  {
    thread_context_t* thread = polymec_malloc(sizeof(thread_context_t));
    thread->id = i;
    thread->pool = pool;
#if USE_PTHREADS
    int stat = pthread_create(&pool->threads[i], &pool->thread_attr, thread_life, thread);
#else
    int stat = thrd_create(&pool->threads[i], thread_life, thread);
#endif
    if (stat != 0)
      polymec_error("thread_pool: could not create thread %d.", i);
  }

  // Unlock the queue and wait for the threads to start up.
#if USE_PTHREADS
  pthread_mutex_unlock(&pool->queue_lock);
#else
  mtx_unlock(&pool->queue_lock);
#endif
  int num_started;
  do 
  {
    num_started = 0;
    for (int i = 0; i < num_threads; ++i)
      num_started += (pool->thread_started[i]) ? 1 : 0;
  }
  while (num_started < num_threads);
  polymec_free(pool->thread_started); // <-- no longer needed.

  return pool;
}

void thread_pool_free(thread_pool_t* pool)
{
  if (pool->work_remaining != 0)
    polymec_error("thread_pool_free() called within thread_pool_execute().");

#if USE_PTHREADS
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
#else
  // Tell the threads we're shutting down.
  mtx_lock(&pool->queue_lock);
  pool->shutting_down = true;
  cnd_broadcast(&pool->queue_cond);
  mtx_unlock(&pool->queue_lock);

  // Wait for all the threads to join.
  for (int i = 0; i < pool->num_threads; ++i)
    thrd_join(pool->threads[i], 0);

  cnd_destroy(&pool->finished_cond);
  cnd_destroy(&pool->queue_cond);
  mtx_destroy(&pool->finished_lock);
  mtx_destroy(&pool->queue_lock);
#endif
  polymec_free(pool->threads);
  ptr_slist_free(pool->queue);
  polymec_free(pool);
}

int thread_pool_num_threads(thread_pool_t* pool)
{
  return pool->num_threads;
}

void thread_pool_schedule(thread_pool_t* pool,
                          void* context,
                          void (*do_work)(void*))
{
  if ((pool->work_remaining != 0) || (pool->work_in_progress != 0))
    polymec_error("thread_pool_schedule() called within thread_pool_execute().");

  thread_work_t* work = polymec_malloc(sizeof(thread_work_t));
  work->context = context;
  work->do_work = do_work;
#if USE_PTHREADS
  pthread_mutex_lock(&pool->queue_lock);
  ptr_slist_append_with_dtor(pool->queue, work, polymec_free);
  pthread_mutex_unlock(&pool->queue_lock);
#else
  mtx_lock(&pool->queue_lock);
  ptr_slist_append_with_dtor(pool->queue, work, polymec_free);
  mtx_unlock(&pool->queue_lock);
#endif
}

void thread_pool_execute(thread_pool_t* pool)
{
#if USE_PTHREADS
  pthread_mutex_lock(&pool->finished_lock);

  // Signal the threads that there's work in the pool.
  pthread_mutex_lock(&pool->queue_lock);
  pool->work_remaining = pool->queue->size;
  pool->work_in_progress = 0;
  pthread_cond_broadcast(&pool->queue_cond);
  pthread_mutex_unlock(&pool->queue_lock);

  // Now wait for the threads to finish.
  pthread_cond_wait(&pool->finished_cond, &pool->finished_lock);
  pthread_mutex_unlock(&pool->finished_lock);
#else
  mtx_lock(&pool->finished_lock);

  // Signal the threads that there's work in the pool.
  mtx_lock(&pool->queue_lock);
  pool->work_remaining = pool->queue->size;
  cnd_broadcast(&pool->queue_cond);
  mtx_unlock(&pool->queue_lock);

  // Now wait for the threads to finish.
  cnd_wait(&pool->finished_cond, &pool->finished_lock);
  mtx_unlock(&pool->finished_lock);
#endif
}

