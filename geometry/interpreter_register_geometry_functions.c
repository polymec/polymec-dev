#include "core/boundary_cell_map.h"
#include "core/constant_st_func.h"
#include "geometry/interpreter_register_geometry_functions.h"
#include "geometry/cubic_lattice.h"
#include "geometry/create_cubic_lattice_mesh.h"
#include "geometry/generate_random_points.h"

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
      lua_pushstring(lua, bounds_names[i]);
      lua_gettable(lua, 4); // Reads name from top, replaces with bounds[name].
      if (!lua_isnumber(lua, -1))
      {
        lua_pushstring(lua, "x1, x2, y1, y2, z1, z2, must all be numbers.");
        lua_error(lua);
        return LUA_ERRRUN;
      }
      bounds_values[i] = lua_tonumber(lua, -1);
      lua_pop(lua, 1); 
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

  // Create the mesh.
  mesh_t* mesh = create_cubic_lattice_mesh_with_bbox(nx, ny, nz, &bbox);

  // Tag its faces.
  tag_cubic_lattice_mesh_faces(mesh, nx, ny, nz, "x1", "x2", "y1", "y2", "z1", "z2");

  // Push the mesh onto the stack.
  lua_pushmesh(lua, mesh);
  return 1;
}

static int cubic_lattice_periodic_bc(lua_State* lua)
{
  // Check the argument.
  int num_args = lua_gettop(lua);
  if (num_args != 2)
  {
    lua_pushstring(lua, "Arguments must be 2 boundary mesh (face) tags.");
    lua_error(lua);
    return LUA_ERRRUN;
  }
  for (int i = 1; i <= 2; ++i)
  {
    if (!lua_isstring(lua, i))
    {
      lua_pushfstring(lua, "Argument %d must be a face tag.", i);
      lua_error(lua);
      return LUA_ERRRUN;
    }
  }

  const char* tag1 = lua_tostring(lua, 1);
  const char* tag2 = lua_tostring(lua, 2);

  // Based on the names of the tags, we can decide which periodic BC 
  // mapping function to use.
  periodic_bc_t* bc;
  if (!strcmp(tag1, "x1"))
  {
    if (!strcmp(tag2, "x2"))
    {
      lua_pushfstring(lua, "Periodic boundary maps from 'x1' to '%s' (must be 'x2').", tag2);
      lua_error(lua);
      return LUA_ERRRUN;
    }
    bc = cubic_lattice_x_periodic_bc_new(tag1, tag2);
  }
  else if (!strcmp(tag1, "x2"))
  {
    if (!strcmp(tag2, "x1"))
    {
      lua_pushfstring(lua, "Periodic boundary maps from 'x2' to '%s' (must be 'x1').", tag2);
      lua_error(lua);
      return LUA_ERRRUN;
    }
    bc = cubic_lattice_x_periodic_bc_new(tag1, tag2);
  }
  else if (!strcmp(tag1, "y1"))
  {
    if (!strcmp(tag2, "y2"))
    {
      lua_pushfstring(lua, "Periodic boundary maps from 'y1' to '%s' (must be 'y2').", tag2);
      lua_error(lua);
      return LUA_ERRRUN;
    }
    bc = cubic_lattice_y_periodic_bc_new(tag1, tag2);
  }
  else if (!strcmp(tag1, "y2"))
  {
    if (!strcmp(tag2, "y1"))
    {
      lua_pushfstring(lua, "Periodic boundary maps from 'y2' to '%s' (must be 'y1').", tag2);
      lua_error(lua);
      return LUA_ERRRUN;
    }
    bc = cubic_lattice_y_periodic_bc_new(tag1, tag2);
  }
  else if (!strcmp(tag1, "z1"))
  {
    if (!strcmp(tag2, "z2"))
    {
      lua_pushfstring(lua, "Periodic boundary maps from 'z1' to '%s' (must be 'z2').", tag2);
      lua_error(lua);
      return LUA_ERRRUN;
    }
    bc = cubic_lattice_z_periodic_bc_new(tag1, tag2);
  }
  else if (!strcmp(tag1, "z2"))
  {
    if (!strcmp(tag2, "z1"))
    {
      lua_pushfstring(lua, "Periodic boundary maps from 'z2' to '%s' (must be 'z1').", tag2);
      lua_error(lua);
      return LUA_ERRRUN;
    }
    bc = cubic_lattice_z_periodic_bc_new(tag1, tag2);
  }
  lua_pushuserdefined(lua, bc, NULL);
  return 1;
}

static int random_points(lua_State* lua)
{
  // Check the arguments.
  int num_args = lua_gettop(lua);
  if (num_args != 2)
  {
    lua_pushstring(lua, "Invalid arguments. Usage:\npoints = random_points(N, bounding_box) OR\npoints = random_points(N, density, bounding_box)");
    lua_error(lua);
    return LUA_ERRRUN;
  }
  // Get the arguments.
  int N = (int)lua_tonumber(lua, 1);
  sp_func_t* density = NULL;
  bbox_t* bbox = NULL;
  if (num_args == 2)
  {
    if (!lua_isboundingbox(lua, 2))
    {
      lua_pushstring(lua, "Second argument must be a bounding box.");
      lua_error(lua);
      return LUA_ERRRUN;
    }
    bbox = lua_toboundingbox(lua, 2);
    double one = 1.0;
    density = constant_sp_func_new(1, &one);
  }
  else
  {
    if (!lua_isscalarfunction(lua, 2))
    {
      lua_pushstring(lua, "Second argument must be a scalar function.");
      lua_error(lua);
      return LUA_ERRRUN;
    }
    st_func_t* density_t = lua_toscalarfunction(lua, 2);
    density = st_func_freeze(density_t, 0.0);
    if (!lua_isboundingbox(lua, 3))
    {
      lua_pushstring(lua, "Third argument must be a bounding box.");
      lua_error(lua);
      return LUA_ERRRUN;
    }
    bbox = lua_toboundingbox(lua, 3);
  }

  point_t* points = malloc(sizeof(point_t) * N);
  generate_random_points(random, density, bbox, N, points);

  // Return the point list.
  lua_pushpointlist(lua, points, N);
  return 1;
}

void interpreter_register_geometry_functions(interpreter_t* interp)
{
  interpreter_register_function(interp, "cubic_lattice_mesh", cubic_lattice_mesh);
  interpreter_register_function(interp, "cubic_lattice_periodic_bc", cubic_lattice_periodic_bc);
  interpreter_register_function(interp, "random_points", random_points);
}

