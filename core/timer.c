// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/polymec.h"
#include "core/options.h"
#include "core/timer.h"
#include "core/array.h"

struct polymec_timer_t 
{
  char* name;

  // Timing data.
  double accum_time, timestamp;

  // Hierarchy information.
  polymec_timer_t* parent;
  ptr_array_t* children;
};

// Globals.
static bool use_timers = false;
static ptr_array_t* all_timers = NULL;
static polymec_timer_t* current_timer = NULL;

static void polymec_timer_free(polymec_timer_t* timer)
{
  string_free(timer->name);
  ptr_array_free(timer->children);
  polymec_free(timer);
}

static polymec_timer_t* polymec_timer_new(const char* name)
{
  ASSERT(use_timers == true);
  polymec_timer_t* t = polymec_malloc(sizeof(polymec_timer_t));
  t->name = string_dup(name);
  t->accum_time = 0.0;
  t->timestamp = MPI_Wtime();
  t->parent = current_timer;
  t->children = ptr_array_new();

  // Register the timer with the global list of all timers.
  if (all_timers == NULL)
    all_timers = ptr_array_new();
  ptr_array_append_with_dtor(all_timers, t, DTOR(polymec_timer_free));
  return t;
}

polymec_timer_t* polymec_timer_get(const char* name)
{
  polymec_timer_t* t = NULL;
  if (current_timer == NULL)
  {
    // Do we need timers?
    options_t* options = options_argv();
    char* timers_on = options_value(options, "timers");
    if ((timers_on != NULL) && 
        ((strcmp(timers_on, "1") == 0) || 
         (strcasecmp(timers_on, "yes") == 0) ||
         (strcasecmp(timers_on, "on") == 0) ||
         (strcasecmp(timers_on, "true") == 0)))
    {
      use_timers = true;
      log_debug("polymec: Enabled timers.");
    }

    if (use_timers)
    {
      // We don't have any timers yet! So make one.
      current_timer = polymec_timer_new(name);
      t = current_timer;
    }
  }
  else
  {
    if (use_timers)
    {
      // Search the current timer for a child of this name.
      for (int i = 0; i < current_timer->children->size; ++i)
      {
        polymec_timer_t* child = current_timer->children->data[i];
        if (strcmp(name, child->name) == 0)
        {
          t = child;
          break;
        }
      }
      if (t == NULL)
        t = polymec_timer_new(name);
    }
  }
  return t;
}

void polymec_timer_start(polymec_timer_t* timer)
{
  if (use_timers)
  {
    // This timer becomes the "current" timer.
    current_timer = timer;
    timer->timestamp = MPI_Wtime();
  }
}

void polymec_timer_stop(polymec_timer_t* timer)
{
  if (use_timers)
  {
    if (current_timer != all_timers->data[0])
    {
      // This timer's parent becomes the "current" timer.
      current_timer = timer->parent;
    }
    double t = MPI_Wtime();
    timer->accum_time += t - timer->timestamp;
    timer->timestamp = MPI_Wtime();
  }
}

void polymec_timer_report()
{
  if (use_timers)
  {
    // This is currently just a stupid enumeration to test that reporting 
    // is properly triggered.
    for (int i = 0; i < all_timers->size; ++i)
    {
      polymec_timer_t* t = all_timers->data[i];
      printf("Timer %s: %g s\n", t->name, t->accum_time);
    }
  }
}

void polymec_timer_finalize()
{
  // Now we delete all the timers! Since they're all stored in an array, 
  // we can delete the array and be done with it.
  if (all_timers != NULL)
  {
    ptr_array_free(all_timers);
    all_timers = NULL;
    current_timer = NULL;
  }
}

