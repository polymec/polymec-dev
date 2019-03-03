// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/polymec.h"
#include "core/options.h"
#include "core/timer.h"
#include "core/array.h"
#include "core/serializer.h"

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
static int mpi_rank = -1;
static int mpi_nproc = -1;
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

void polymec_enable_timers()
{
  use_timers = true;
  log_debug("polymec: Enabled timers.");

  // Set a default timer file.
  polymec_set_timer_file("timer_report.txt");

  // Record our MPI rank and number of processes.
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_nproc);

  // Get threading information.
}

const char* polymec_timer_file()
{
  return (const char*)timer_report_file;
}

void polymec_set_timer_file(const char* timer_file)
{
  strncpy(timer_report_file, timer_file, FILENAME_MAX);
}

polymec_timer_t* polymec_timer_get(const char* name)
{
  static bool first_time = true;

  polymec_timer_t* t = NULL;
  if (first_time)
  {
    // Do we need timers?
    options_t* options = options_argv();
    char* timers_on = options_value(options, "timers");
    if ((timers_on != NULL) && string_as_boolean(timers_on))
    {
      polymec_enable_timers();

      // Were we given a specific report file?
      char* timer_file = options_value(options, "timer_file");
      if (timer_file != NULL)
        polymec_set_timer_file(timer_file);
    }

    first_time = false;
  }

  if (use_timers)
  {
    if (all_timers == NULL)
    {
      // We don't have any timers yet! So make one.
      t = polymec_timer_new(name);
    }
    else
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
    polymec_timer_stop(current_timer); // Last one!
  }
}

static void report_timer(polymec_timer_t* root,
                         polymec_timer_t* t,
                         int indentation,
                         FILE* file)
{
  double percent = (double)(100.0 * t->accum_time / root->accum_time);
  char call_string[9];
  if (t->count > 1)
    strcpy(call_string, "calls");
  else
    strcpy(call_string, "call");
  int name_len = (int)strlen(t->name);
  fprintf(file, "%*s%*s%10.4f s  %5.1f%%  %10lld %s\n", indentation + name_len, t->name,
          45 - indentation - name_len, " ", t->accum_time, percent, t->count, call_string);
  size_t num_children = t->children->size;
  for (size_t i = 0; i < num_children; ++i)
  {
    polymec_timer_t* child = t->children->data[i];
    report_timer(root, child, indentation+1, file);
  }
}

static void polymec_timer_finalize()
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

static size_t polymec_timer_byte_size(void* obj)
{
  polymec_timer_t* timer = obj;

  size_t byte_size = 2 * sizeof(size_t) +                  // metadata
                     sizeof(char) * strlen(timer->name) +  // name
                     sizeof(double) * 2 +                  // timings
                     sizeof(unsigned long long);           // counts

  // Find the byte size of the children.
  for (size_t i = 0; i < timer->children->size; ++i)
    byte_size += polymec_timer_byte_size(timer->children->data[i]);

  return byte_size;
}

static void* polymec_timer_byte_read(byte_array_t* bytes, size_t* offset)
{
  // Read in the timer's name length and number of children.
  size_t name_len;
  byte_array_read_size_ts(bytes, 1, &name_len, offset);
  size_t num_children;
  byte_array_read_size_ts(bytes, 1, &num_children, offset);

  // Read in the timer's name.
  char name[name_len+1];
  byte_array_read_chars(bytes, name_len, name, offset);
  name[name_len] = '\0';

  // Create the timer. Avoid using polymec_timer_new, as that constructor
  // registers stuff in the local process's list of resources.
  polymec_timer_t* timer = polymec_malloc(sizeof(polymec_timer_t));
  timer->name = string_dup(name);
  timer->parent = NULL;
  timer->children = ptr_array_new();

  // Read in the timing data.
  byte_array_read_doubles(bytes, 1, &timer->accum_time, offset);
  byte_array_read_doubles(bytes, 1, &timer->timestamp, offset);
  byte_array_read_unsigned_long_longs(bytes, 1, &timer->count, offset);

  // Read in all the children and set their parent.
  for (size_t i = 0; i < num_children; ++i)
  {
    polymec_timer_t* child = polymec_timer_byte_read(bytes, offset);
    child->parent = timer;
    ptr_array_append_with_dtor(timer->children, child, DTOR(polymec_timer_free));
  }

  return timer;
}

static void polymec_timer_byte_write(void* obj, byte_array_t* bytes, size_t* offset)
{
  polymec_timer_t* timer = obj;

  // Write the timer's name length and number of children.
  size_t name_len = strlen(timer->name);
  byte_array_write_size_ts(bytes, 1, &name_len, offset);
  byte_array_write_size_ts(bytes, 1, &timer->children->size, offset);

  // Write the timer's name.
  byte_array_write_chars(bytes, name_len, timer->name, offset);

  // Write the timing data.
  byte_array_write_doubles(bytes, 1, &timer->accum_time, offset);
  byte_array_write_doubles(bytes, 1, &timer->timestamp, offset);
  byte_array_write_unsigned_long_longs(bytes, 1, &timer->count, offset);

  // Write the children's data.
  for (size_t i = 0; i < timer->children->size; ++i)
    polymec_timer_byte_write(timer->children->data[i], bytes, offset);
}

static serializer_t* timer_serializer()
{
  return serializer_new("timer",
                        polymec_timer_byte_size,
                        polymec_timer_byte_read,
                        polymec_timer_byte_write,
                        NULL);
}

void polymec_timer_report()
{
  if (use_timers)
  {
    log_debug("polymec: writing timer report file '%s'.", timer_report_file);

    FILE* report_file = NULL;
    if (mpi_rank == 0)
    {
      report_file = fopen(timer_report_file, "w");
      if (report_file == NULL)
        polymec_error("Could not open file '%s' for writing!", timer_report_file);

      // Print a header for the timer report.
      fprintf(report_file, "-----------------------------------------------------------------------------------\n");
      fprintf(report_file, "                                   Timer summary:\n");
      fprintf(report_file, "-----------------------------------------------------------------------------------\n");
      fprintf(report_file, "Invocation: %s\n", polymec_invocation());
      time_t invoc_time = polymec_invocation_time();
      fprintf(report_file, "At: %s", ctime(&invoc_time));
      fprintf(report_file, "-----------------------------------------------------------------------------------\n");
      fprintf(report_file, "%s%*s%s\n", "Name:", 49-5, " ", "Time:     Percent:     Count:");
      fprintf(report_file, "-----------------------------------------------------------------------------------\n");
    }

    for (int p = 0; p < mpi_nproc; ++p)
    {
      if (mpi_rank == 0)
      {
        // Get timer information for this rank.
        polymec_timer_t* timer = NULL;
        if (p == mpi_rank)
          timer = all_timers->data[0]; // timer data is available locally.
        else
        {
          // We receive timer data from rank p.
          int recv_size;
          MPI_Recv(&recv_size, 1, MPI_INT, p, p, MPI_COMM_WORLD,
                   MPI_STATUS_IGNORE);
          byte_array_t* bytes = byte_array_new();
          byte_array_resize(bytes, (size_t)recv_size);
          MPI_Recv(bytes->data, recv_size, MPI_UINT8_T, p, p, MPI_COMM_WORLD,
                   MPI_STATUS_IGNORE);

          // Unpack it.
          serializer_t* s = timer_serializer();
          size_t offset = 0;
          timer = serializer_read(s, bytes, &offset);

          // Clean up.
          byte_array_free(bytes);
        }

        // Write out a textual representation of the timer.
        if (mpi_nproc > 1)
        {
          fprintf(report_file, "\n------------\n");
          fprintf(report_file, "Rank %6d:\n", p);
          fprintf(report_file, "------------\n");
        }
        report_timer(timer, timer, 0, report_file);

        // Destroy any off-process timer data.
        if (p != mpi_rank)
          polymec_timer_free(timer);
      }
      else if (mpi_rank == p)
      {
        // We pack up our local timer into a send buffer.
        serializer_t* s = timer_serializer();
        byte_array_t* bytes = byte_array_new();
        size_t offset = 0;
        serializer_write(s, all_timers->data[0], bytes, &offset);

        // Send the data to rank 0.
        int send_size = (int)bytes->size;
        MPI_Send(&send_size, 1, MPI_INT, 0, mpi_rank, MPI_COMM_WORLD);
        MPI_Send(bytes->data, (int)bytes->size, MPI_UINT8_T, 0, mpi_rank, MPI_COMM_WORLD);
        byte_array_free(bytes);
      }
    }

    if (mpi_rank == 0)
      fclose(report_file);

    polymec_timer_finalize();
  }
}

