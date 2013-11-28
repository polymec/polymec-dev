// Copyright 2012-2013 Jeffrey Johnson.
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
  int work_remaining;
  pthread_mutex_t queue_lock, finished_lock;
  pthread_cond_t queue_cond, finished_cond;
  pthread_t* threads;
  bool* thread_started;
  pthread_attr_t thread_attr;
};

static void* thread_life(void* context)
{
  thread_context_t* thread = context;
  thread_pool_t* pool = thread->pool;
  pthread_setcancelstate(PTHREAD_CANCEL_ENABLE, NULL);
  pool->thread_started[thread->id] = true;
  
  while (true)
  {
    // Wait on the condition variable for work to show up in the queue.
    pthread_mutex_lock(&pool->queue_lock);
    pthread_cond_wait(&pool->queue_cond, &pool->queue_lock);

    while (!ptr_slist_empty(pool->queue))
    {
      // Grab the next item in the queue and process it.
      void (*dtor)(void*);
      thread_work_t* work = ptr_slist_pop(pool->queue, &dtor);

      // Surrender the queue lock.
      pthread_mutex_unlock(&pool->queue_lock);

      // Do the work.
      work->do_work(work->context);

      // Delete the work unit.
      dtor(work);

      // Mark the job as completed.
      pthread_mutex_lock(&pool->finished_lock);
      --pool->work_remaining;

      // If there's no more work, tell the host thread that we're finished.
      if (pool->work_remaining == 0)
        pthread_cond_signal(&pool->finished_cond);
      pthread_mutex_unlock(&pool->finished_lock);
    }
  }
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
  pool->thread_started = malloc(sizeof(bool) * num_threads);
  memset(pool->thread_started, 0, sizeof(bool) * num_threads);
  pool->queue = ptr_slist_new();
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
  if (pool->work_remaining != 0)
    polymec_error("thread_pool_free() called within thread pool job.");

  // Cancel all the threads -- since thread_pool_execute is blocking, 
  // we are guaranteed that this will work unless something is REALLY 
  // wrong.
  for (int i = 0; i < pool->num_threads; ++i)
    pthread_cancel(pool->threads[i]);

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
    polymec_error("thread_pool_schedule() called within thread pool job.");

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
  pool->work_remaining = pool->queue->size;
  pthread_mutex_lock(&pool->queue_lock);
  pthread_cond_broadcast(&pool->queue_cond);
  pthread_mutex_unlock(&pool->queue_lock);

  // Now wait for the threads to finish.
  pthread_cond_wait(&pool->finished_cond, &pool->finished_lock);
  pthread_mutex_unlock(&pool->finished_lock);
}

