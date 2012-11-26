#include "geometry/interpreter_register_geometry_functions.h"
#include "geometry/create_cubic_lattice_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

// Lua stuff.
#include "lua.h"
#include "lualib.h"
#include "lauxlib.h"

static int cubic_lattice_mesh(lua_State* lua)
{
  // Check the arguments.
  int num_args = lua_gettop(lua);
  if ((num_args != 3) && (num_args != 4))
  {
    lua_pushstring(lua, "Invalid arguments. Usage:\nmesh = cubic_lattice_mesh(nx, ny, nz) OR\nmesh = cubic_lattice_mesh(nx, ny, nz, bounds)");
    lua_error(lua);
    return LUA_ERRRUN;
  }
  // Get the arguments.
  int nx = (int)lua_tonumber(lua, 1);
  int ny = (int)lua_tonumber(lua, 2);
  int nz = (int)lua_tonumber(lua, 3);
  if ((nx <= 0) || (ny <= 0) || (nz <= 0))
  {
    lua_pushstring(lua, "nx, ny, and nz must all be positive.");
    lua_error(lua);
    return LUA_ERRRUN;
  }

  // Bounding box? 
  bbox_t bbox = {.x1 = 0.0, .x2 = 1.0, .y1 = 0.0, .y2 = 1.0, .z1 = 0.0, .z2 = 1.0};
  if (num_args == 4)
  {
    if (!lua_istable(lua, 4))
    {
      lua_pushstring(lua, "bounds must be a table containing x1, x2, y1, y2, z1, z2.");
      lua_error(lua);
      return LUA_ERRRUN;
    }

    // Look for x1, x2, y1, y2, z1, z2 in the table.
    const char* bounds_names[] = {"x1", "x2", "y1", "y2", "z1", "z2"};
    double bounds_values[6];
    for (int i = 0; i < 6; ++i)
    {
      lua_pushstring(lua, bounds_names[i]); // moves bounds -> 5.
      lua_gettable(lua, 5); // Reads name from top, replaces with bounds[name].
      if (!lua_isnumber(lua, 1))
      {
        lua_pushstring(lua, "x1, x2, y1, y2, z1, z2, must all be numbers.");
        lua_error(lua);
        return LUA_ERRRUN;
      }
      bounds_values[i] = lua_tonumber(lua, 1);
      lua_pop(lua, 1); // Moves bounds -> 4.
    }
    bbox.x1 = bounds_values[0];
    bbox.x2 = bounds_values[1];
    bbox.y1 = bounds_values[2];
    bbox.y2 = bounds_values[3];
    bbox.z1 = bounds_values[4];
    bbox.z2 = bounds_values[5];

    // Now check the bounds.
    if (bbox.x1 >= bbox.x2)
    {
      lua_pushstring(lua, "x1 must be less than x2.");
      lua_error(lua);
      return LUA_ERRRUN;
    }
    if (bbox.y1 >= bbox.y2)
    {
      lua_pushstring(lua, "y1 must be less than y2.");
      lua_error(lua);
      return LUA_ERRRUN;
    }
    if (bbox.z1 >= bbox.z2)
    {
      lua_pushstring(lua, "z1 must be less than z2.");
      lua_error(lua);
      return LUA_ERRRUN;
    }
  }

  // Create the mesh and push it onto the stack.
  mesh_t* mesh = create_cubic_lattice_mesh_with_bbox(nx, ny, nz, &bbox);
  lua_pushmesh(lua, mesh);
  return 1;
}

void interpreter_register_geometry_functions(interpreter_t* interp)
{
  interpreter_register_function(interp, "cubic_lattice_mesh", cubic_lattice_mesh);
}

#ifdef __cplusplus
}
#endif

