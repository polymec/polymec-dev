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
  unsigned long long count;

  // Hierarchy information.
  polymec_timer_t* parent;
  ptr_array_t* children;
};

// Globals.
int mpi_rank = -1;
static bool use_timers = false;
static ptr_array_t* all_timers = NULL;
static polymec_timer_t* current_timer = NULL;
static char timer_report_file[FILENAME_MAX];

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
  t->count = 0;
  t->timestamp = MPI_Wtime();
  t->parent = current_timer;
  t->children = ptr_array_new();

  // Make sure our parent records us.
  if (t->parent != NULL)
    ptr_array_append(t->parent->children, t);

  // Register the timer with the global list of all timers.
  if (all_timers == NULL)
    all_timers = ptr_array_new();
  ptr_array_append_with_dtor(all_timers, t, DTOR(polymec_timer_free));
  return t;
}

polymec_timer_t* polymec_timer_get(const char* name)
{
  polymec_timer_t* t = NULL;
  if (all_timers == NULL)
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

      // Were we given a specific report file?
      char* timer_file = options_value(options, "timer_file");
      if (timer_file != NULL)
        strncpy(timer_report_file, timer_file, FILENAME_MAX);
      else
        strcpy(timer_report_file, "timer_report.txt");

      // Record our MPI rank.
      MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

      // Get threading information.
    }

    if (use_timers)
    {
      // We don't have any timers yet! So make one.
      t = polymec_timer_new(name);
    }
  }
  else
  {
    if (use_timers)
    {
      if (current_timer != NULL) // A timer is running...
      {
        if (strcmp(name, current_timer->name) == 0)
          t = current_timer;
        else
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
      else
      {
        // polymec_init, which runs the 'polymec' timer, hasn't been called.
        polymec_error("polymec_timer_get: tried to create a timer without first calling\n"
                      "polymec_init. Please make sure to call polymec_init first.");
      }
    }
  }
  return t;
}

void polymec_timer_start(polymec_timer_t* timer)
{
  if (use_timers)
  {
    if (timer == current_timer)
      polymec_error("polymec_timer_start: Can't start timer %s, which has already been started.", timer->name);

    // This timer becomes the "current" timer.
    current_timer = timer;
    timer->timestamp = MPI_Wtime();
    ++(timer->count);
  }
}

void polymec_timer_stop(polymec_timer_t* timer)
{
  if (use_timers)
  {
    if (current_timer == NULL)
      polymec_error("polymec_timer_stop: Can't stop timer %s: no timers are running.", timer->name);

    if (timer != current_timer)
    {
      polymec_error("polymec_timer_stop: Can't stop timer %s, which isn't the currently running one.\n"
                    "(current timer is %s)", timer->name, current_timer->name);
    }

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

void polymec_timer_stop_all()
{
  if (use_timers)
  {
    while ((current_timer != NULL) && (current_timer != all_timers->data[0]))
      polymec_timer_stop(current_timer);
  }
}

static void report_timer(polymec_timer_t* t, int indentation, FILE* file)
{
  polymec_timer_t* root = all_timers->data[0];
  double percent = (double)(100.0 * t->accum_time / root->accum_time);
  char call_string[9];
  if (t->count > 1)
    strcpy(call_string, "calls");
  else
    strcpy(call_string, "call");
  int name_len = strlen(t->name);
  fprintf(file, "%*s%*s%10.4f s  %5.1f%%  %10lld %s\n", indentation + name_len, t->name, 
          45 - indentation - name_len, " ", t->accum_time, percent, t->count, call_string);
  int num_children = t->children->size;
  for (int i = 0; i < num_children; ++i)
  {
    polymec_timer_t* child = t->children->data[i];
    report_timer(child, indentation+1, file);
  }
}

void polymec_timer_report()
{
  if (use_timers)
  {
    log_debug("polymec: writing timer report file '%s'.", timer_report_file);
    if (mpi_rank == 0)
    {
      // This is currently just a stupid enumeration to test that reporting 
      // is properly triggered.
      FILE* report_file = fopen(timer_report_file, "w");
      if (report_file == NULL)
        polymec_error("Could not open file '%s' for writing!", timer_report_file);
      fprintf(report_file, "-----------------------------------------------------------------------------------\n");
      fprintf(report_file, "                                   Timer summary:\n");
      fprintf(report_file, "-----------------------------------------------------------------------------------\n");
      fprintf(report_file, "Invocation: %s\n", polymec_invocation());
      time_t invoc_time = polymec_invocation_time();
      fprintf(report_file, "At: %s", ctime(&invoc_time));
      fprintf(report_file, "-----------------------------------------------------------------------------------\n");

      // Print a header.
      fprintf(report_file, "%s%*s%s\n", "Name:", 49-5, " ", "Time:     Percent:     Count:");
      fprintf(report_file, "-----------------------------------------------------------------------------------\n");


      polymec_timer_t* t = all_timers->data[0];
      report_timer(t, 0, report_file);
      fclose(report_file);
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

