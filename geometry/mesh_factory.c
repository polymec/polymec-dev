// mesh_factory.c - Implementations of interpreter functions for generating
// meshes.

#include "core/polymec.h"
#include "core/mesh.h"
#include "core/interpreter.h"
#include "geometry/cubic_lattice.h"
#include "geometry/create_cubic_lattice_mesh.h"
#include "geometry/create_voronoi_mesh.h"

// Lua stuff.
#include "lua.h"
#include "lualib.h"
#include "lauxlib.h"

int mesh_factory_cubic_lattice(lua_State* lua)
{
  // Check the arguments.
  int num_args = lua_gettop(lua);
  if ((num_args != 3) && (num_args != 4))
  {
    lua_pushstring(lua, "Invalid arguments. Usage:\n"
                  "mesh = cubic_lattice_mesh(nx, ny, nz) OR\n"
                  "mesh = cubic_lattice_mesh(nx, ny, nz, bounds)");
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

  // Pop all the previous arguments off the stack.
  lua_pop(lua, lua_gettop(lua));

  // Create the mesh.
  mesh_t* mesh = create_cubic_lattice_mesh_with_bbox(nx, ny, nz, &bbox);

  // Tag its faces.
  tag_cubic_lattice_mesh_faces(mesh, nx, ny, nz, "x1", "x2", "y1", "y2", "z1", "z2");

  // Push the mesh onto the stack.
  lua_pushmesh(lua, mesh);
  return 1;
}

int mesh_factory_cubic_lattice_periodic_bc(lua_State* lua)
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
  periodic_bc_t* bc = NULL;
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

int mesh_factory_voronoi(lua_State* lua)
{
  // Check the arguments.
  int num_args = lua_gettop(lua);
  if ((num_args != 1) || !lua_ispointlist(lua, 1))
  {
    lua_pushstring(lua, "Invalid argument(s). Usage:\n"
                        "mesh = mesh_factory.voronoi(generators)");
    lua_error(lua);
    return LUA_ERRRUN;
  }

  // Get the generators.
  int num_generators;
  point_t* generators = lua_topointlist(lua, 1, &num_generators);

  // Create the mesh.
  mesh_t* mesh = create_voronoi_mesh(generators, num_generators,
                                               NULL, 0);

  // Push the mesh onto the stack.
  lua_pushmesh(lua, mesh);
  return 1;
}

#if 0
int mesh_factory_cvt(lua_State* lua)
{
  // Check the arguments.
  int num_args = lua_gettop(lua);
  if ((num_args != 2) || !lua_istable(lua, 1) || !lua_istable(lua, 2))
  {
    lua_pushstring(lua, "Invalid argument(s). Usage:\n"
                        "mesh = mesh_factory.cvt(surface, options) OR\n"
                        "mesh = mesh_factory.cvt(surfaces, options)");
    lua_error(lua);
    return LUA_ERRRUN;
  }

  // Get the options for the Voronoi tessellation.
  lua_getfield(lua, 2, "num_interior_points");
  if (!lua_isnumber(lua, -1))
  {
    lua_pushstring(lua, "Option 'num_interior_points' must be a number.");
    lua_error(lua);
    return LUA_ERRRUN;
  }
  int num_interior_points = (int)lua_tonumber(lua, -1);
  lua_pop(lua, 1);
  lua_getfield(lua, 2, "iteration_method");
  if (!lua_isnil(lua, -1) && !lua_isstring(lua, -1))
  {
    lua_pushstring(lua, "Option 'iteration_method' must be a string.");
    lua_error(lua);
    return LUA_ERRRUN;
  }
  char iteration_method[1024];
  strcpy(iteration_method, "lloyd");
  if (!lua_isnil(lua, -1))
    strcpy(iteration_method, lua_tostring(lua, -1));
  lua_pop(lua, 1);
  if (!strcasecmp(iteration_method, "lloyd"))
  {
    lua_pushstring(lua, "Supported iteration methods are 'lloyd'.");
    lua_error(lua);
    return LUA_ERRRUN;
  }
  lua_getfield(lua, 2, "num_iterations");
  if (!lua_isnumber(lua, -1))
  {
    lua_pushstring(lua, "Option 'num_iterations' must be a number.");
    lua_error(lua);
    return LUA_ERRRUN;
  }
  int num_iterations = lua_tonumber(lua, -1);
  if (num_iterations < 1)
  {
    lua_pushstring(lua, "Option 'num_iterations' must be positive.");
    lua_error(lua);
    return LUA_ERRRUN;
  }
  lua_pop(lua, 1);

  // Determine whether the first argument is a surface or a table containing
  // several surfaces. We do this by traversing the table and looking for 
  // certain fields. In the case of a list of surfaces, we make a list of 
  // surface names.
  bool is_surface = false, found_points = false, found_normals = false;
  bool is_surface_list = false, is_invalid = false;

  str_array_t* surface_names = str_array_new();
  ptr_array_t* surface_points = ptr_array_new();
  ptr_array_t* surface_normals = ptr_array_new();
  ptr_array_t* surface_tags = ptr_array_new();
  lua_pushnil(lua);
  while (lua_next(lua, 1))
  {
    static const int key_index = -2;
    static const int val_index = -1;
    const char* key = lua_tostring(lua, key_index);
    if (!is_surface_list && !strcmp("points"))
    {
      if (!lua_ispointlist(lua, val_index))
      {
        is_invalid = true;
        break;
      }
      found_points = true;
    }
    else if (!is_surface_list && !strcmp("normals"))
    {
      if (!lua_isvectorlist(lua, val_index))
      {
        is_invalid = true;
        break;
      }
      found_normals = true;
    }
    else
    {
      // Something else is in there. It must be a table.
      if (!lua_istable(lua, key_index))
      {
        is_invalid = true;
        break;
      }
      lua_pushnil(lua);
      bool found_points = false, found_normals = false;
      while (lua_next(lua, -2))
      {
        const char* key = lua_tostring(lua, key_index);
        if (!is_surface_list && !strcmp("points"))
        {
          if (!lua_ispointlist(lua, val_index))
          {
            is_invalid = true;
            break;
          }
          found_points = true;
        }
        else if (!is_surface_list && !strcmp("normals"))
        {
          if (!lua_isvectorlist(lua, val_index))
          {
            is_invalid = true;
            break;
          }
          found_points = true;
        }
      }
      if (is_invalid) break;
      if (found_points && found_normals)
      {

        str_ptr_unordered_map_insert_with_kv_dtor(surfaces, strdup(key), lua_tofree_name_and_surface);
        is_surface_list = true;
      }
    }
    if (!is_surface_list && (found_points && found_normals)) 
      is_surface = true;
    lua_pop(lua, 1);
  }

  if (is_invalid || (!is_surface && !is_surface_list)) 
  {
    string_slist_free(surface_names);
    lua_pushstring(lua, "Argument 1 must be a surface (a table with points and normals fields) or a list of surfaces.");
    lua_error(lua);
    return LUA_ERRRUN;
  }

  // In the case of a single surface, we name it "all".
  if (is_surface)
    string_slist_append(surface_names, "all");

  // Process the surface or surfaces into a list of points with normals and 
  // tags.
  int num_points;
  point_t* surface_points;
  ptr_array_t* normal_lists;
  string_array_t* tags;
  string_slist_node_t* s_iter = surface_names->front;
  while (iter != NULL)
  {
    char* surface_name = s_iter->value;
    lua_getfield(lua, 1, surface_name);

    s_iter = s_iter->next;
  }
  lua_pop(lua, 1);

  // Clean up a bit.
  string_slist_free(surface_names);

  // For now, we do Lloyd iteration, fixing the boundary points and moving 
  // the interior points.
  mesh_t* mesh = NULL;
  int iteration = 0;
  do
  {
    // Construct a new tessellation.

    // Go through all the interior points and compute their centroids, 
    // moving each interior point to its centroid.

    // Delete the existing mesh so that we can retessellate.
//    mesh_free(mesh);
  }
  while (iteration < num_iterations);

  // Push the mesh onto the stack.
  lua_pushmesh(lua, mesh);
  return 1;
}
#endif

