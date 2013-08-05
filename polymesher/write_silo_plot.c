// This implements polymesher's capability for writing Silo plots of meshes.

#include <string.h>
#include "core/polymec.h"
#include "core/interpreter.h"
#include "io/silo_io.h"

// Lua stuff.
#include "lua.h"
#include "lualib.h"
#include "lauxlib.h"

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
  if (((num_args == 2) && (!lua_ismesh(lua, 1) || !lua_isstring(lua, 2))) || 
      ((num_args == 3) && (!lua_ismesh(lua, 1) || !lua_istable(lua, 2) || !lua_isstring(lua, 3))) || 
      ((num_args != 2) && (num_args != 3)))
  {
    lua_pushstring(lua, "write_silo_plot: invalid arguments. Usage:\n"
                        "write_silo_plot(mesh, filename) OR\n"
                        "write_silo_plot(mesh, fields, filename)");
    lua_error(lua);
    return LUA_ERRRUN;
  }

  // Get the argument(s).
  mesh_t* mesh = lua_tomesh(lua, 1);
  bool has_fields = (num_args == 3) ? lua_istable(lua, 2) : false;
  char* filename = (num_args == 3) ? strdup(lua_tostring(lua, 3)) : strdup(lua_tostring(lua, 2));

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
        if (num_vals != mesh->num_cells)
        {
          lua_pop(lua, 2);
          polymec_error("write_silo_plot: a scalar field has %d values (should have %d).", num_vals, mesh->num_cells);
        }
      }
      else
      {
        int num_vals;
        vector_t* vals = lua_tovectorlist(lua, val_index, &num_vals);
        if (num_vals != mesh->num_cells)
        {
          lua_pop(lua, 2);
          polymec_error("write_silo_plot: a vector field has %d values (should have %d).", num_vals, mesh->num_cells);
        }
      }
      lua_pop(lua, 1);
    }
  }

  log_info("Writing SILO plot with prefix '%s'...", filename);

  // Write the mesh to a plot file.
  io_interface_t* plot = silo_plot_io_new(MPI_COMM_SELF, 0, false);
  io_open(plot, filename, ".", IO_WRITE);
  io_dataset_t* dataset = io_dataset_new("default");
  io_dataset_put_mesh(dataset, mesh);

  // Stick in a couple of diagnostic fields.
  double volume[mesh->num_cells];
  for (int c = 0; c < mesh->num_cells; ++c)
    volume[c] = mesh->cells[c].volume;
  io_dataset_put_field(dataset, "volume", volume, 1, MESH_CELL, true);

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
        io_dataset_put_field(dataset, field_name, field_data, 1, MESH_CELL, true);
      }
      else
      {
        ASSERT(lua_isvectorlist(lua, val_index));
        int num_vals;
        vector_t* vector_data = lua_tovectorlist(lua, val_index, &num_vals);
        double* field_data = malloc(sizeof(double) * 3 * num_vals);
        for (int i = 0; i < num_vals; ++i)
        {
          field_data[3*i  ] = vector_data[i].x;
          field_data[3*i+1] = vector_data[i].y;
          field_data[3*i+2] = vector_data[i].z;
        }
        io_dataset_put_field(dataset, field_name, field_data, 3, MESH_CELL, false);
      }
      lua_pop(lua, 1);
    }
  }

  io_append_dataset(plot, dataset);
  io_close(plot);

  // Clean up.
  io_free(plot);

  return 1;
}

