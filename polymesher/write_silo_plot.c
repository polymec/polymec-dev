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

// This implements polymesher's capability for writing Silo plots of meshes.

#include <string.h>
#include "core/polymec.h"
#include "core/interpreter.h"

// Lua stuff.
#include "lua.h"
#include "lualib.h"
#include "lauxlib.h"

#include "core/write_silo.h"

// write_silo_plot(args) -- This function writes a given mesh to a file 
// on disk. Arguments (passed in a table according to Chapter 5.3 of the 
// Lua reference manual) are:
//
// mesh -> mesh object 
// filename -> name of the file to write (1 file only)
int write_silo_plot(lua_State* lua)
{
  // Check the arguments.
  int num_args = lua_gettop(lua);
  bool is_mesh = false;
  if (lua_ismesh(lua, 1))
  {
    is_mesh = true;
    if (((num_args == 2) && (!lua_ismesh(lua, 1) || !lua_isstring(lua, 2))) || 
        ((num_args == 3) && (!lua_ismesh(lua, 1) || !lua_istable(lua, 2) || !lua_isstring(lua, 3))) || 
        ((num_args != 2) && (num_args != 3)))
    {
      return luaL_error(lua, "write_silo_plot: invalid arguments. Usage:\n"
                        "write_silo_plot(mesh, filename) OR\n"
                        "write_silo_plot(mesh, fields, filename)");
    }
  }
  else
  {
    if (((num_args == 2) && (!lua_ispointlist(lua, 1) || !lua_isstring(lua, 2))) || 
        ((num_args == 3) && (!lua_ispointlist(lua, 1) || !lua_istable(lua, 2) || !lua_isstring(lua, 3))) || 
        ((num_args != 2) && (num_args != 3)))
    {
      return luaL_error(lua, "write_silo_plot: invalid arguments. Usage:\n"
                        "write_silo_plot(points, filename) OR\n"
                        "write_silo_plot(points, fields, filename)");
    }
  }

  // Get the argument(s).
  mesh_t* mesh = NULL;
  point_t* points = NULL;
  int N;
  if (is_mesh)
  {
    mesh = lua_tomesh(lua, 1);
    N = mesh->num_cells;
    ASSERT(mesh != NULL);
  }
  else
    points = lua_topointlist(lua, 1, &N);
  bool has_fields = (num_args == 3) ? lua_istable(lua, 2) : false;
  char* filename = (num_args == 3) ? string_dup(lua_tostring(lua, 3)) : string_dup(lua_tostring(lua, 2));

  // Check the table of fields if it's there.
  if (has_fields)
  {
    lua_pushnil(lua);
    while (lua_next(lua, 2))
    {
      static const int key_index = -2;
      static const int val_index = -1;
      bool key_is_string = lua_isstring(lua, key_index);
      bool val_is_field = lua_issequence(lua, val_index) || lua_isvectorlist(lua, val_index);
      if (!key_is_string || !val_is_field)
      {
        lua_pop(lua, 2);
        polymec_error("write_silo_plot: argument 2 must be a table mapping field names to values.");
      }
      if (lua_issequence(lua, val_index))
      {
        int num_vals;
        double* vals = lua_tosequence(lua, val_index, &num_vals);
        if (num_vals != N)
        {
          lua_pop(lua, 2);
          polymec_error("write_silo_plot: a scalar field has %d values (should have %d).", num_vals, N);
        }
        vals = NULL;
      }
      else
      {
        int num_vals;
        vector_t* vals = lua_tovectorlist(lua, val_index, &num_vals);
        if (num_vals != N)
        {
          lua_pop(lua, 2);
          polymec_error("write_silo_plot: a vector field has %d values (should have %d).", num_vals, N);
        }
        vals = NULL;
      }
      lua_pop(lua, 1);
    }
  }

  // Construct a set of fields.
  string_ptr_unordered_map_t* fields = string_ptr_unordered_map_new();
  if (is_mesh)
  {
    double* volume = malloc(sizeof(double) * N);
    for (int c = 0; c < N; ++c)
      volume[c] = mesh->cell_volumes[c];
    string_ptr_unordered_map_insert_with_v_dtor(fields, "volume", volume, DTOR(free));
  }

  // Stick in any other fields.
  if (has_fields)
  {
    lua_pushnil(lua);
    while (lua_next(lua, 2))
    {
      static const int key_index = -2;
      static const int val_index = -1;
      const char* field_name = (const char*)lua_tostring(lua, key_index);
      if (lua_issequence(lua, val_index))
      {
        int num_vals;
        double* field_data = lua_tosequence(lua, val_index, &num_vals);
        string_ptr_unordered_map_insert(fields, (char*)field_name, field_data);
      }
      else
      {
        ASSERT(lua_isvectorlist(lua, val_index));
        int num_vals;
        vector_t* vector_data = lua_tovectorlist(lua, val_index, &num_vals);
        double* Fx = malloc(sizeof(double) * num_vals);
        double* Fy = malloc(sizeof(double) * num_vals);
        double* Fz = malloc(sizeof(double) * num_vals);
        for (int i = 0; i < num_vals; ++i)
        {
          Fx[i] = vector_data[i].x;
          Fy[i] = vector_data[i].y;
          Fz[i] = vector_data[i].z;
        }
        char Fx_name[128], Fy_name[128], Fz_name[128];
        snprintf(Fx_name, 128, "%s_x", field_name);
        snprintf(Fy_name, 128, "%s_y", field_name);
        snprintf(Fz_name, 128, "%s_z", field_name);
        string_ptr_unordered_map_insert(fields, Fx_name, Fx);
        string_ptr_unordered_map_insert(fields, Fy_name, Fy);
        string_ptr_unordered_map_insert(fields, Fz_name, Fz);
        free(Fx);
        free(Fy);
        free(Fz);
      }
      lua_pop(lua, 1);
    }
  }

  log_info("Writing SILO plot with prefix '%s'...", filename);

  // Write the thing to a plot file.
  if (is_mesh)
  {
    write_silo_mesh(mesh, NULL, NULL, NULL, fields, filename, ".", 0, 0.0, 
                    MPI_COMM_SELF, 1, 0);
  }
  else
  {
    write_silo_points(points, N, fields, filename, ".", 0, 0.0, 
                      MPI_COMM_SELF, 1, 0);
  }

  // Clean up.
  string_ptr_unordered_map_free(fields);

  return 1;
}

