// This implements polymesher's capability for writing VTK plots of meshes.

#include <string.h>
#include "core/polymec.h"
#include "core/interpreter.h"
#include "io/vtk_plot_io.h"

// Lua stuff.
#include "lua.h"
#include "lualib.h"
#include "lauxlib.h"

// write_tough_mesh(args) -- This function writes a given mesh to a file 
// on disk. Arguments (passed in a table according to Chapter 5.3 of the 
// Lua reference manual) are:
//
// filename -> name of the file to write (1 file only)
// format -> 'T2', 'T+' (mesh format to use)
// mesh -> mesh object 
// inactive_tag -> the tag within the mesh object that denotes inactive elements.
int write_vtk_plot(lua_State* lua)
{
  // Check the arguments.
  int num_args = lua_gettop(lua);
  if ((num_args != 2) || !lua_ismesh(lua, 1) || !lua_isstring(lua, 2))
  {
    lua_pushstring(lua, "write_tough_mesh: invalid arguments. Usage:\n"
                        "write_tough_mesh(mesh, filename)");
    lua_error(lua);
    return LUA_ERRRUN;
  }

  // Get the argument(s).
  mesh_t* mesh = lua_tomesh(lua, 1);
  char* filename = strdup(lua_tostring(lua, 2));

  log_info("Writing VTK plot with prefix '%s'...", filename);

  // Write the mesh to a plot file.
  io_interface_t* plot = vtk_plot_io_new(MPI_COMM_SELF, 0, false);
  io_open(plot, filename, ".", IO_WRITE);
  io_dataset_t* dataset = io_dataset_new("default");
  io_dataset_put_mesh(dataset, mesh);

  // Stick in a couple of diagnostic fields.
  double volume[mesh->num_cells];
  for (int c = 0; c < mesh->num_cells; ++c)
    volume[c] = mesh->cells[c].volume;
  io_dataset_put_field(dataset, "volume", volume, 1, MESH_CELL, true);

  io_append_dataset(plot, dataset);
  io_close(plot);

  // Clean up.
  io_free(plot);

  return 1;
}

