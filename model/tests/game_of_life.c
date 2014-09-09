// Copyright (c) 2012-2014, Jeffrey N. Johnson
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

#include "core/polymec.h"
#include "core/silo_file.h"
#include "geometry/create_uniform_mesh.h"
#include "model/model.h"

typedef struct
{
  int nx, ny;
  mesh_t* grid;
  real_t* state;
} gol_t;

static const char gol_desc[] = "Game of Life model\n"
  "This model demonstrates a parallel version of the Game of Life. For more\n"
  "details on Life and its variations, see\n"
  "See http://en.wikipedia.org/wiki/Conway%27s_Game_of_Life#Notable_Life_programs\n";

static void gol_read_custom_input(void* context, const char* input, options_t* options)
{
}

static void gol_init(void* context, real_t t)
{
}

static real_t gol_max_dt(void* context, real_t t, char* reason)
{
  return 1.0;
}

static real_t gol_advance(void* context, real_t max_dt, real_t t)
{
  return 1.0;
}

static void gol_load(void* context, const char* file_prefix, const char* directory, real_t* t, int step)
{
  *t = 1.0 * step;
}

static void gol_save(void* context, const char* file_prefix, const char* directory, real_t t, int step)
{
  gol_t* gol = context;
  ASSERT(t == 1.0*step);

  int rank;
  MPI_Comm_rank(gol->grid->comm, &rank);

  silo_file_t* silo = silo_file_new(gol->grid->comm, file_prefix, directory, 1, 0, step, t);
  silo_file_write_mesh(silo, "grid", gol->grid);
  silo_file_write_scalar_cell_field(silo, "life", "grid", gol->state);
  real_t* rank_field = polymec_malloc(sizeof(real_t) * gol->grid->num_cells);
  for (int i = 0; i < gol->grid->num_cells; ++i)
    rank_field[i] = 1.0 * rank;
  silo_file_write_scalar_cell_field(silo, "rank", "grid", rank_field);
  polymec_free(rank_field);
  silo_file_close(silo);
}

static void gol_dtor(void* context)
{
  polymec_free(context);
}

static model_t* gol_ctor()
{
  gol_t* gol = polymec_malloc(sizeof(gol_t));
  model_vtable vtable = {.read_custom_input = gol_read_custom_input,
                         .init = gol_init,
                         .max_dt = gol_max_dt,
                         .advance = gol_advance,
                         .load = gol_load,
                         .save = gol_save,
                         .dtor = gol_dtor};
  docstring_t* gol_doc = docstring_from_string(gol_desc);
  return model_new("game_of_life", gol, vtable, gol_doc);
}

// Main program.
int main(int argc, char* argv[])
{
  return model_main("game_of_life", gol_ctor, argc, argv);
}
