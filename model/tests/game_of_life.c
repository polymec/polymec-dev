// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/polymec.h"
#include "core/kd_tree.h"
#include "core/unordered_set.h"
#include "core/silo_file.h"
#include "core/text_buffer.h"
#include "geometry/create_uniform_mesh.h"
#include "model/model.h"
#include "model/stencil.h"

typedef struct
{
  // State information.
  mesh_t* grid;
  stencil_t* stencil;
  real_t* state;

  // The Game of Life rule we use.
  int alive_nums[9]; // alive_nums[x] is 1 if x neighbors keeps a cell alive, 0 if not.
  int born_nums[9];  // born_nums[x] is 1 if x neighbors births a cell, 0 if not.

  // The information we read from our custom input file.
  int x_min, x_max, y_min, y_max;
  int_tuple_unordered_set_t* xy_pairs;

  // Expected difference in life pattern.
  int expected_diff;
} gol_t;

static const char gol_desc[] = "Game of Life model (game_of_life)\n"
  "This model demonstrates a parallel version of the Game of Life. For more\n"
  "details on Life and its variations, see\n"
  "See http://en.wikipedia.org/wiki/Conway%27s_Game_of_Life#Notable_Life_programs\n\n"
  "game_of_life takes Life 1.05 (*.LIF) files as input and takes time steps\n"
  "of size 1. Important command line options are:\n\n"
  "  t1=T          (optional) - floor(T) specifies the start time for the\n"
  "                             simulation.\n"
  "  t2=T          (required) - floor(T) specifies the end time for the\n"
  "                             simulation.\n"
  "  expect_diff=D (optional) - computes the difference in the life pattern and\n"
  "                             compares it to D.\n"
  "  save_every=N  (optional) - causes a save (.silo) file to be written every\n"
  "                             N steps.\n"
  "  pause=D       (optional) - Pauses for D seconds before execution.\n\n";

static void gol_read_custom_input(void* context, const char* input, options_t* options)
{
  // Because we parse custom input, we rely on command line options for 
  // some of our required inputs.
  if (options_value(options, "t2") == NULL)
    polymec_error("t2 must be given as an argument.");

  gol_t* gol = context;

  // We read .LIF/.LIFE files in Life 1.05 format, populating an initial state.
  // See http://psoup.math.wisc.edu/mcell/ca_files_formats.html#Life%201.05 for details.
  text_buffer_t* text = text_buffer_from_string(input);

  // First line should be a header with the file format.
  int pos = 0, len;
  char* raw_line;
  char line[120];
  bool got_one = text_buffer_next_nonempty(text, &pos, &raw_line, &len);
  if (got_one == false)
    polymec_error("Input file was empty.");
  string_copy_from_raw(raw_line, len+1, line);
  if (strstr((const char*)line, "#Life 1.05") == NULL)
    polymec_error("Invalid file format -- should be Life 1.05.");

  // Up to 22 description lines should follow.
  int num_desc_lines = 0;
  while (num_desc_lines < 22)
  {
    got_one = text_buffer_next_nonempty(text, &pos, &raw_line, &len);
    if (got_one == false)
      polymec_error("Life 1.05 file has no content!");
    string_copy_from_raw(raw_line, len+1, line);
    if ((line[0] != '#') || ((line[1] != 'D')))
      break;
    ++num_desc_lines;
  }

  // Optional rule specification. Default is 23/3.
  bool rule_specified = false;
  memset(gol->alive_nums, 0, sizeof(int) * 9);
  memset(gol->born_nums, 0, sizeof(int) * 9);
  gol->alive_nums[2] = gol->alive_nums[3] = 1;
  gol->born_nums[3] = 1;
  if ((line[0] == '#') && (line[1] == 'N'))
    rule_specified = true;
  else if ((line[0] == '#') && (line[1] == 'R'))
  {
    bool has_slash = (strstr(line, "/") != NULL);
    if (has_slash)
    {
      int spos = 0, slen;
      char* str;
      string_next_token(line, "/", &spos, &str, &slen);
      if (slen > 8)
        polymec_error("Invalid living string: %s", str);
      for (int i = 0; i < slen; ++i)
      {
        if (!string_is_number(&str[i]))
          polymec_error("Invalid living string: %s", str);
        int count = atoi(&str[i]);
        if ((count < 0) || (count > 8))
          polymec_error("Illegal value for living count: %d\n", count);
        gol->alive_nums[count] = 1;
      }
      string_next_token(line, "/", &spos, &str, &slen);
      if (slen > 8)
        polymec_error("Invalid birthing string: %s", str);
      for (int i = 0; i < slen; ++i)
      {
        if (!string_is_number(&str[i]))
          polymec_error("Invalid birthing string: %s", str);
        int count = atoi(&str[i]);
        if ((count < 0) || (count > 8))
          polymec_error("Illegal value for birthing count: %d\n", count);
        gol->born_nums[count] = 1;
      }
      rule_specified = true;
    }
    else
      polymec_error("Invalid rule specification: %s", line);
  }
  
  // If the last line was a rule specification, we need to read the next line.
  if (rule_specified)
    got_one = text_buffer_next_nonempty(text, &pos, &raw_line, &len);
  if (got_one == false)
    polymec_error("Life 1.05 file has no cell blocks!");
  string_copy_from_raw(raw_line, len+1, line);
  
  // Now we move on to cell blocks.
  gol->xy_pairs = int_tuple_unordered_set_new();
  gol->x_min = INT_MAX;
  gol->x_max = -INT_MAX;
  gol->y_min = INT_MAX;
  gol->y_max = -INT_MAX;
  char num_buf[16];
  int num_cell_blocks = 0;
  while (got_one)
  {
    // A cell block header follows the form
    // #P X Y
    if ((len < 6) || (line[0] != '#') || (line[1] != 'P'))
      polymec_error("Invalid cell block header: %s", line);
    int spos = 0, slen;
    char* str;
    bool found_coords = string_next_token(&line[3], " ", &spos, &str, &slen);
    if (!found_coords)
      polymec_error("Invalid cell block offset: %s", line);
    string_copy_from_raw(str, slen+1, num_buf);
    int t = string_trim(num_buf);
    if (!string_is_number(&num_buf[t]))
      polymec_error("Invalid x position of cell block: %s", &num_buf[t]);
    int block_x = atoi(num_buf);
    found_coords = string_next_token(&line[3], " ", &spos, &str, &slen);
    if (!found_coords)
      polymec_error("Invalid cell block offset: %s", line);
    string_copy_from_raw(str, slen+1, num_buf);
    t = string_trim(num_buf);
    if (!string_is_number(&num_buf[t]))
      polymec_error("Invalid y position of cell block: %s", &num_buf[t]);
    int block_y = atoi(num_buf);
    
    if (block_x < gol->x_min)
      gol->x_min = block_x;
    if (block_y < gol->y_min)
      gol->y_min = block_y;

    // Now process the lines in the cell block.
    int x, y = block_y;
    do
    {
      got_one = text_buffer_next_nonempty(text, &pos, &raw_line, &len);
      if (!got_one) break;
      string_copy_from_raw(raw_line, len+1, line);
      for (int i = 0; i < len; ++i)
      {
        char c = line[i];
        x = block_x + i;
        if (c == '*') 
        {
          if (x > gol->x_max)
            gol->x_max = x;
          if (y > gol->y_max)
            gol->y_max = y;
          int* xy = int_tuple_new(2);
          xy[0] = x;
          xy[1] = y;
          int_tuple_unordered_set_insert_with_dtor(gol->xy_pairs, xy, int_tuple_free);
        }
        else if (c != '.')
          break;
      }
      ++y;
    } while (got_one && ((len < 6) || ((line[0] != '#') && (line[1] != 'P'))) );
    ++num_cell_blocks;
  }
  text_buffer_free(text);

  log_detail("game_of_life: Read %d cell blocks from input.", num_cell_blocks);

  char* exp_diff = options_value(options, "expect_diff");
  if (exp_diff != NULL)
  {
    if (!string_is_number(exp_diff))
      polymec_error("Invalid expected diff: %s\n", exp_diff);
    gol->expected_diff = atoi(exp_diff);
    if (gol->expected_diff < 0)
      polymec_error("Invalid expected diff: %s\n", exp_diff);
  }
  else
    gol->expected_diff = -1;
}

static void gol_init(void* context, real_t t)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  int nprocs;
  MPI_Comm_size(comm, &nprocs);

  gol_t* gol = context;
  ASSERT(gol->xy_pairs->size > 0);
  ASSERT(gol->x_min < gol->x_max);
  ASSERT(gol->y_min < gol->y_max);

  // Create the grid and the state vector. Make sure we're not going to run 
  // out of memory!
  int nx = gol->x_max - gol->x_min + 1;
  int ny = gol->y_max - gol->y_min + 1;
  int total_num_cells = nx * ny;
  if ((total_num_cells / nprocs) > 10000000)
    polymec_error("Excessively huge grid (%d cells on %d processes)!", total_num_cells, nprocs);

  log_detail("game_of_life: Creating %d x %d grid...", nx, ny);
  bbox_t bbox = {.x1 = 1.0*gol->x_min - 0.5, .x2 = 1.0*gol->x_max + 0.5,
                 .y1 = 1.0*gol->y_min - 0.5, .y2 = 1.0*gol->y_max + 0.5,
                 .z1 = -0.5, .z2 = 0.5};
  gol->grid = create_uniform_mesh(comm, nx, ny, 1, &bbox);
  gol->stencil = cell_edge_stencil_new(gol->grid);
  gol->state = polymec_malloc(sizeof(real_t) * gol->grid->num_cells);
  memset(gol->state, 0, sizeof(real_t) * gol->grid->num_cells);

  // Place the living cells into the grid. This is a little silly, but this is what 
  // we get for using a completely unstructured grid!
  int num_living_cells = gol->xy_pairs->size;
  log_detail("game_of_life: Placing %d living cells...", num_living_cells);
  kd_tree_t* tree = kd_tree_new(gol->grid->cell_centers, gol->grid->num_cells);
  int pos = 0, *xy;
  while (int_tuple_unordered_set_next(gol->xy_pairs, &pos, &xy))
  {
    point_t x = {.x = 1.0*xy[0], .y = 1.0*xy[1]};
    int cell_index = kd_tree_nearest(tree, &x);
    ASSERT(cell_index >= 0);
    ASSERT(cell_index < gol->grid->num_cells);
    gol->state[cell_index] = 1.0;
  }
  kd_tree_free(tree);
}

static real_t gol_max_dt(void* context, real_t t, char* reason)
{
  return 1.0;
}

static real_t gol_advance(void* context, real_t max_dt, real_t t)
{
  gol_t* gol = context;

  // Make sure that the state field is exchanged all parallel-like.
  stencil_exchange(gol->stencil, gol->state, 1, 0, MPI_REAL_T);

  // Now go over each cell and update according to living neighbor counts.
  for (int cell = 0; cell < gol->grid->num_cells; ++cell)
  {
    int pos = 0, neighbor;
    real_t count = 0.0, weight;
    while (stencil_next(gol->stencil, cell, &pos, &neighbor, &weight))
    {
      if (neighbor >= 0)
        count += gol->state[neighbor];
    }
    int icount = (int)count;
    int living = (int)gol->state[cell];
    if (living)
    {
      if (!gol->alive_nums[icount])
        gol->state[cell] = 0.0;
    }
    else if (gol->born_nums[icount])
      gol->state[cell] = 1.0;
  }
  return 1.0;
}

static void gol_load(void* context, const char* file_prefix, const char* directory, real_t* t, int step)
{
  gol_t* gol = context;
  silo_file_t* silo = silo_file_open(gol->grid->comm, file_prefix, directory, 0, step, t);
  *t = 1.0 * step;

  ASSERT(gol->grid == NULL);
  ASSERT(gol->state == NULL);
  gol->grid = silo_file_read_mesh(silo, "grid");
  gol->stencil = cell_edge_stencil_new(gol->grid);
  gol->state = silo_file_read_scalar_cell_field(silo, "state", "grid");
  silo_file_close(silo);
}

static void gol_save(void* context, const char* file_prefix, const char* directory, real_t t, int step)
{
  gol_t* gol = context;

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

static void gol_finalize(void* context, int step, real_t t)
{
  gol_t* gol = context;

  // Check for any expected difference.
  if (gol->expected_diff >= 0)
  {
    int local_diff = 0;
    int* xy = int_tuple_new(2);
    for (int c = 0; c < gol->grid->num_cells; ++c)
    {
      int val = (int)gol->state[c];
      point_t* x = &gol->grid->cell_centers[c];
      xy[0] = (int)x->x, xy[1] = (int)x->y;
      if (int_tuple_unordered_set_contains(gol->xy_pairs, xy))
        local_diff += abs(val - 1);
      else
        local_diff += val;
    }
    int_tuple_free(xy);
    int global_diff;
    MPI_Allreduce(&local_diff, &global_diff, 1, MPI_INT, MPI_SUM, gol->grid->comm);

    log_urgent("game_of_life: expected diff: %d", gol->expected_diff);
    log_urgent("game_of_life: computed diff: %d", global_diff);
    if (gol->expected_diff == global_diff)
      log_urgent("game_of_life: expected == computed");
    else
      log_urgent("game_of_life: expected != computed");
  }

  stencil_free(gol->stencil);
  mesh_free(gol->grid);
  polymec_free(gol->state);
}

static void gol_dtor(void* context)
{
  gol_t* gol = context;
  if (gol->xy_pairs != NULL)
    int_tuple_unordered_set_free(gol->xy_pairs);
  polymec_free(gol);
}

static model_t* gol_ctor()
{
  gol_t* gol = polymec_malloc(sizeof(gol_t));
  gol->grid = NULL;
  gol->stencil = NULL;
  gol->state = NULL;
  gol->xy_pairs = NULL;
  memset(gol->alive_nums, 0, sizeof(int) * 9);
  memset(gol->born_nums, 0, sizeof(int) * 9);
  model_vtable vtable = {.read_custom_input = gol_read_custom_input,
                         .init = gol_init,
                         .max_dt = gol_max_dt,
                         .advance = gol_advance,
                         .load = gol_load,
                         .save = gol_save,
                         .finalize = gol_finalize,
                         .dtor = gol_dtor};
  docstring_t* gol_doc = docstring_from_string(gol_desc);
  return model_new("game_of_life", gol, vtable, gol_doc);
}

// Main program.
int main(int argc, char* argv[])
{
  return model_main("game_of_life", gol_ctor, argc, argv);
}
