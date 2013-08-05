#include "core/polymec.h"
#include "core/point.h"
#include "core/interpreter.h"

// Lua stuff.
#include "lua.h"
#include "lualib.h"
#include "lauxlib.h"

int point_factory_cubic_lattice(lua_State* lua)
{
  // Check the arguments.
  int num_args = lua_gettop(lua);
  if ((num_args != 4) && (num_args != 5))
  {
    lua_pushstring(lua, "Invalid arguments. Usage:\n"
                  "points = point_factory.cubic_lattice(nx, ny, nz, ng) OR\n"
                  "mesh = cubic_lattice_mesh(nx, ny, nz, bounding_box, ng)");
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

  // Number of ghost points.
  int ng = 0;
  if (num_args == 4)
  {
    if (!lua_isnumber(lua, 4))
    {
      lua_pushstring(lua, "ng must be an integer.");
      lua_error(lua);
      return LUA_ERRRUN;
    }
    ng = (int)lua_tonumber(lua, 4);
  }
  else
    ng = (int)lua_tonumber(lua, 5);
  if (ng < 0)
  {
    lua_pushstring(lua, "ng must be non-negative.");
    lua_error(lua);
    return LUA_ERRRUN;
  }

  // Bounding box? 
  bbox_t bbox = {.x1 = 0.0, .x2 = 1.0, .y1 = 0.0, .y2 = 1.0, .z1 = 0.0, .z2 = 1.0};
  if (num_args == 5)
  {
    if (!lua_isboundingbox(lua, 4))
    {
      lua_pushstring(lua, "bounding_box must be a bounding box.");
      lua_error(lua);
      return LUA_ERRRUN;
    }
    bbox = *lua_toboundingbox(lua, 4);

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

  // Pop all the previous arguments off the stack.
  lua_pop(lua, lua_gettop(lua));

  // Create the lattice of points.
  int num_points = (nx + 2*ng) * (ny + 2*ng) * (nz + 2*ng), offset = 0;
  point_t* points = malloc(sizeof(point_t) * num_points);
  double dx = (bbox.x2 - bbox.x1) / nx;
  double dy = (bbox.y2 - bbox.y1) / ny;
  double dz = (bbox.z2 - bbox.z1) / nz;
  for (int i = -ng; i < nx+ng; ++i)
  {
    double xi = (i+0.5) * dx;
    for (int j = -ng; j < ny+ng; ++j)
    {
      double yj = (j+0.5) * dy;
      for (int k = -ng; k < nz+ng; ++k, ++offset)
      {
        double zk = (k+0.5) * dz;
        points[offset].x = xi;
        points[offset].y = yj;
        points[offset].z = zk;
      }
    }
  }

  // Push the points onto the stack.
  lua_pushpointlist(lua, points, num_points);
  return 1;
}

