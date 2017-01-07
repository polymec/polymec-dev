// Copyright (c) 2012-2017, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "simple_model.h"

// This is an extremely simple model implementation that we use to test 
// the model interface.
typedef struct
{
  int index; // Simulation index.
  real_t dt; // Time step size.
  MPI_Comm comm;
} simple_t;

static void simple_set_global_comm(void* context, MPI_Comm comm)
{
  simple_t* s = context;
  s->comm = comm;
}

static void simple_read_input(void* context,
                              interpreter_t* interp,
                              options_t* opts)
{
  simple_t* s = context;
  s->index = (int)(interpreter_get_number(interp, "index"));
  s->dt = interpreter_get_number(interp, "dt");
  log_debug("simple: read input (index = %d, dt = %g)", s->index, s->dt);
}

static void name_simple_file(simple_t* s, char* filename)
{
  // If we are the master rank for our global communicator, write out a file 
  // named simple.x.y.z, where x is the simulation index, y is the number of 
  // processes belonging to MPI_COMM_WORLD, and z is the number of processes 
  // belonging to our global communicator.
  int g_nproc, l_nproc;
  MPI_Comm_size(MPI_COMM_WORLD, &g_nproc);
  MPI_Comm_size(s->comm, &l_nproc);

  snprintf(filename, FILENAME_MAX-1, "simple.%d.%d.%d", s->index, g_nproc, l_nproc);
}

static void simple_init(void* context, real_t t)
{
  simple_t* s = context;

  int rank;
  MPI_Comm_rank(s->comm, &rank);
  log_debug("simple: rank %d reporting for simulation %d.", rank, s->index);

  // Delete the file corresponding to this simulation if it exists.
  char filename[FILENAME_MAX];
  name_simple_file(s, filename);

  if (file_exists(filename))
  {
    log_debug("simple: removing file %s", filename);
    remove(filename);
  }
}

static real_t simple_advance(void* context, real_t max_dt, real_t t)
{
  simple_t* s = context;
  return MIN(max_dt, s->dt);
}

static void simple_finalize(void* context, int step, real_t t)
{
  simple_t* s = context;

  MPI_Barrier(s->comm);

  // If we're the master rank, create our file.
  int rank;
  MPI_Comm_rank(s->comm, &rank);
  if (rank == 0)
  {
    char filename[FILENAME_MAX];
    name_simple_file(s, filename);

    if (file_exists(filename))
      polymec_error("simple_finalize: Found file %s -- shouldn't exist!", filename);
    else
    {
      log_debug("simple: writing file %s", filename);
      FILE* file = fopen(filename, "w");
      fprintf(file, "simple sim finished at t = %g in %d steps! (dt = %g) \n", t, step-1, s->dt);
      fclose(file);
    }
  }
}

static void simple_dtor(void* context)
{
  simple_t* s = context;
  polymec_free(s);
}

model_t* simple_new()
{
  simple_t* s = polymec_malloc(sizeof(simple_t));
  model_vtable vtable = {.set_global_comm = simple_set_global_comm,
                         .read_input = simple_read_input,
                         .init = simple_init,
                         .advance = simple_advance,
                         .finalize = simple_finalize,
                         .dtor = simple_dtor};
  return model_new("simple", s, vtable, 
                   docstring_from_string("simple model"), 
                   MODEL_MPI);
}

