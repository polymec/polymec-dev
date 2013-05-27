#include "core/boundary_cell_map.h"
#include "core/constant_st_func.h"
#include "geometry/interpreter_register_geometry_functions.h"
#include "geometry/cubic_lattice.h"
#include "geometry/create_cubic_lattice_mesh.h"
#include "geometry/generate_random_points.h"
#include "geometry/create_unbounded_voronoi_mesh.h"
#include "geometry/create_bounded_voronoi_mesh.h"
#include "geometry/rect_prism.h"
#include "geometry/prune_voronoi_mesh.h"
#include "geometry/bound_voronoi_mesh.h"

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
  if (N <= 0)
  {
    lua_pushstring(lua, "Invalid (nonpositive) number of points.");
    lua_error(lua);
    return LUA_ERRRUN;
  }
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
    ASSERT(bbox != NULL);
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

static int ccp_points(lua_State* lua)
{
  // Check the arguments.
  int num_args = lua_gettop(lua);
  if (num_args != 4)
  {
    lua_pushstring(lua, "Invalid arguments. Usage:\n"
                        "points = ccp_points(Nx, Ny, Nz, bounding_box)");
    lua_error(lua);
    return LUA_ERRRUN;
  }

  // Get the arguments.
  int Nx = (int)lua_tonumber(lua, 1);
  int Ny = (int)lua_tonumber(lua, 2);
  int Nz = (int)lua_tonumber(lua, 3);
  if ((Nx <= 0) || (Ny <= 0) || (Nz <= 0))
  {
    lua_pushstring(lua, "Nx, Ny, and Nz must all be positive.");
    lua_error(lua);
    return LUA_ERRRUN;
  }
  if (!lua_isboundingbox(lua, 4))
  {
    lua_pushstring(lua, "Fourth argument must be a bounding box.");
    lua_error(lua);
    return LUA_ERRRUN;
  }
  bbox_t* bbox = lua_toboundingbox(lua, 4);

  // Create the point list.
  ptr_slist_t* point_list = ptr_slist_new();
  double dx = (bbox->x2 - bbox->x1) / Nx,
         dy = (bbox->y2 - bbox->y1) / Ny,
         dz = (bbox->z2 - bbox->z1) / Nz;
  for (int i = 0; i < Nx; ++i)
  {
    bool x1face = (i > 0);
    bool x2face = (i < Nx-1);
    double x1 = bbox->x1 + i * dx;
    double x2 = x1 + dx;
    for (int j = 0; j < Ny; ++j)
    {
      bool y1face = (j > 0);
      bool y2face = (j < Ny-1);
      double y1 = bbox->x1 + j * dy;
      double y2 = y1 + dy;
      for (int k = 0; k < Nz; ++k)
      {
        bool z1face = (k > 0);
        bool z2face = (k < Ny-1);
        double z1 = bbox->x1 + k * dz;
        double z2 = z1 + dz;

        if (x1face)
        {
          // yz face center
          ptr_slist_append(point_list, point_new(x1, y1+0.5*dy, z1+0.5*dz));

          // nodes on the x1 face.
          if (y1face)
          {
            if (z1face)
              ptr_slist_append(point_list, point_new(x1, y1, z1));
            if (z2face)
              ptr_slist_append(point_list, point_new(x1, y1, z2));
          }
          if (y2face)
          {
            if (z1face)
              ptr_slist_append(point_list, point_new(x1, y2, z1));
            if (z2face)
              ptr_slist_append(point_list, point_new(x1, y2, z2));
          }
        }

        if (x2face)
        {
          // yz face center
          ptr_slist_append(point_list, point_new(x2, y1+0.5*dy, z1+0.5*dz));

          // nodes on the x2 face.
          if (y1face)
          {
            if (z1face)
              ptr_slist_append(point_list, point_new(x2, y1, z1));
            if (z2face)
              ptr_slist_append(point_list, point_new(x2, y1, z2));
          }
          if (y2face)
          {
            if (z1face)
              ptr_slist_append(point_list, point_new(x2, y2, z1));
            if (z2face)
              ptr_slist_append(point_list, point_new(x2, y2, z2));
          }
        }

        if (y1face)
        {
          // xz face center.
          ptr_slist_append(point_list, point_new(x1+0.5*dx, y1, z1+0.5*dz));
        }
        if (y2face)
        {
          // xz face center.
          ptr_slist_append(point_list, point_new(x1+0.5*dx, y2, z1+0.5*dz));
        }

        if (z1face)
        {
          // xy face center.
          ptr_slist_append(point_list, point_new(x1+0.5*dx, y1+0.5*dy, z1));
        }
        if (z2face)
        {
          // xy face center.
          ptr_slist_append(point_list, point_new(x1+0.5*dx, y2+0.5*dy, z2));
        }
      }
    }
  }

  // Pack it up and return it.
  int num_points = point_list->size;
  point_t* points = malloc(sizeof(point_t) * num_points);
  for (int i = 0; i < num_points; ++i)
    point_copy(&points[i], ptr_slist_pop(point_list, NULL));
  ptr_slist_free(point_list);
  lua_pushpointlist(lua, points, num_points);
  return 1;
}

static int unbounded_voronoi_mesh(lua_State* lua)
{
  // Check the arguments.
  int num_args = lua_gettop(lua);
  if ((num_args != 1) || !lua_ispointlist(lua, 1))
  {
    lua_pushstring(lua, "Invalid argument(s). Usage:\n"
                        "mesh = unbounded_voronoi_mesh(generators)");
    lua_error(lua);
    return LUA_ERRRUN;
  }

  // Get the generators.
  int num_generators;
  point_t* generators = lua_topointlist(lua, 1, &num_generators);

  // Create the mesh.
  mesh_t* mesh = create_unbounded_voronoi_mesh(generators, num_generators,
                                               NULL, 0);

  // Push the mesh onto the stack.
  lua_pushmesh(lua, mesh);
  return 1;
}

static int bounded_voronoi_mesh(lua_State* lua)
{
  // Check the arguments.
  int num_args = lua_gettop(lua);
  if ((num_args != 2) || !lua_ispointlist(lua, 1) || !lua_ispointlist(lua, 2))
  {
    lua_pushstring(lua, "Invalid argument(s). Usage:\n"
                        "mesh = bounded_voronoi_mesh(generators, boundary_generators).");
    lua_error(lua);
    return LUA_ERRRUN;
  }

  // Get the generators.
  int num_generators, num_boundary_generators;
  point_t* generators = lua_topointlist(lua, 1, &num_generators);
  point_t* boundary_generators = lua_topointlist(lua, 2, &num_boundary_generators);

  // Create an bounded mesh.
  mesh_t* mesh = create_bounded_voronoi_mesh(generators, num_generators,
                                             boundary_generators, num_boundary_generators,
                                             NULL, 0);

  // Log some information.
  log_detail("Generated bounded Voronoi mesh:");
  log_detail("  %d interior cells, %d boundary cells", num_generators, num_boundary_generators);
  log_detail("  %d faces, %d edges, %d nodes", mesh->num_faces, mesh->num_edges, mesh->num_nodes);

  // Push the mesh onto the stack.
  lua_pushmesh(lua, mesh);
  return 1;
}

#if 0
static int bounded_voronoi_mesh(lua_State* lua)
{
  // Check the arguments.
  int num_args = lua_gettop(lua);
  if ((num_args != 2) || 
      !lua_ispointlist(lua, 1) || 
      (!lua_isscalarfunction(lua, 2) && !lua_isboundingbox(lua, 2)))
  {
    lua_pushstring(lua, "Invalid argument(s). Usage:\n"
                        "mesh = bounded_voronoi_mesh(generators, boundary)\n"
                        "where boundary is a bounding box OR an implicit function\n"
                        "whose 0 level set indicates the boundary.");
    lua_error(lua);
    return LUA_ERRRUN;
  }

  // Get the generators.
  int num_generators;
  point_t* generators = lua_topointlist(lua, 1, &num_generators);

  // Get the implicit function representing the boundary.
  sp_func_t* boundary;
  if (lua_isscalarfunction(lua, 2))
  {
    st_func_t* boundary_t = lua_toscalarfunction(lua, 2);
    boundary = st_func_freeze(boundary_t, 0.0); // Freeze at t = 0.
  }
  else
  {
    bbox_t* bbox = lua_toboundingbox(lua, 2);
    boundary = rect_prism_new_from_bbox(bbox);
  }

  // Kill all generators outside of the boundary.
  int num_culled_generators = num_generators;
  point_t* culled_generators = malloc(sizeof(point_t) * num_generators);
  memcpy(culled_generators, generators, sizeof(point_t) * num_generators);
  for (int i = 0; i < num_generators; ++i)
  {
    double bval;
    sp_func_eval(boundary, &culled_generators[i], &bval);
    if (bval >= 0.0)
    {
      point_t temp = culled_generators[num_generators-1];
      culled_generators[num_generators-1] = culled_generators[i];
      culled_generators[i] = temp;
      --num_culled_generators;
    }
  }

  // Create an unbounded mesh.
  mesh_t* mesh = create_unbounded_voronoi_mesh(culled_generators, num_culled_generators,
                                               NULL, 0);

  // Reflect the points nearest to the boundary across the boundary.
  int num_outer_cells;
  int* outer_cells = mesh_tag(mesh->cell_tags, "outer_cells", &num_outer_cells);
  if (outer_cells == NULL)
    polymec_error("bounded_voronoi_mesh: could not find outer cells. An error occurred.");
  point_t* reflected_generators = malloc(sizeof(point_t) * num_outer_cells);
  for (int c = 0; c < num_outer_cells; ++c)
  {
    point_t* xi = &culled_generators[outer_cells[c]];
    double D;
    sp_func_eval(boundary, xi, &D);
    double grad_f[3];
    sp_func_eval_deriv(boundary, 1, xi, grad_f);
    double grad_mag = sqrt(grad_f[0]*grad_f[0] + grad_f[1]*grad_f[1] + grad_f[2]*grad_f[2]);
    vector_t n = {.x = -grad_f[0] / grad_mag, 
                  .y = -grad_f[1] / grad_mag, 
                  .z = -grad_f[2] / grad_mag};

    point_t xb = {.x = xi->x + 2.0 * D * n.x,
                  .y = xi->y + 2.0 * D * n.y,
                  .z = xi->z + 2.0 * D * n.z};
    reflected_generators[c] = xb;
  }

  // Toss out the mesh and construct a new one with the additional 
  // generators.
  mesh_free(mesh);
  int num_all_generators = num_culled_generators + num_outer_cells;
  point_t* all_generators = malloc(sizeof(point_t) * num_all_generators);
  memcpy(all_generators, culled_generators, num_culled_generators);
  memcpy(all_generators + num_generators, reflected_generators, num_outer_cells);
  mesh = create_unbounded_voronoi_mesh(all_generators, num_all_generators,
                                       NULL, 0);

  // Clean up.
  free(all_generators);
  free(reflected_generators);
  free(culled_generators);

  // Prune the mesh.
  prune_voronoi_mesh(mesh);

  // Push it onto the stack.
  lua_pushmesh(lua, mesh);
  return 1;
}
#endif

static int bound_voronoi_mesh_(lua_State* lua)
{
  // Check the arguments.
  int num_args = lua_gettop(lua);
  if ((num_args != 2) || !lua_ismesh(lua, 1) || 
      (!lua_isscalarfunction(lua, 2) && !lua_isboundingbox(lua, 2)))
  {
    lua_pushstring(lua, "Invalid argument(s). Usage:\n"
                        "bound_voronoi_mesh(mesh, boundary)\n"
                        "where mesh is an unbounded Voronoi mesh and \n"
                        "boundary is an implicit function.");
    lua_error(lua);
    return LUA_ERRRUN;
  }

  // Get the mesh.
  mesh_t* mesh = lua_tomesh(lua, 1);

  // Get the implicit function representing the boundary.
  sp_func_t* boundary;
  if (lua_isscalarfunction(lua, 2))
  {
    st_func_t* boundary_t = lua_toscalarfunction(lua, 2);
    boundary = st_func_freeze(boundary_t, 0.0); // Freeze at t = 0.
  }
  else
  {
    bbox_t* bbox = lua_toboundingbox(lua, 2);
    boundary = rect_prism_new_from_bbox(bbox);
  }

  // Make sure the mesh is unbounded.
  if (!mesh_has_tag(mesh->cell_tags, "outer_cells") || 
      !mesh_has_tag(mesh->edge_tags, "outer_edges"))
  {
    lua_pushstring(lua, "Given mesh is not unbounded (no outer cells/edges found).");
    lua_error(lua);
    return LUA_ERRRUN;
  }

  // Bound the mesh using the boundary function.
  mesh_diff_t* diff = bound_voronoi_mesh(mesh, boundary);
  mesh_diff_apply(diff, mesh);

  // No results.
  return 0;
}

static int prune_voronoi_mesh_(lua_State* lua)
{
  // Check the arguments.
  int num_args = lua_gettop(lua);
  if ((num_args != 1) || !lua_ismesh(lua, 1))
  {
    lua_pushstring(lua, "Invalid argument(s). Usage:\n"
                        "prune_voronoi_mesh(mesh)\n"
                        "where mesh is an unbounded Voronoi mesh.");
    lua_error(lua);
    return LUA_ERRRUN;
  }

  // Get the mesh.
  mesh_t* mesh = lua_tomesh(lua, 1);

  // Make sure the mesh is unbounded.
  if (!mesh_has_tag(mesh->cell_tags, "outer_cells") || 
      !mesh_has_tag(mesh->edge_tags, "outer_edges"))
  {
    lua_pushstring(lua, "Given mesh is not unbounded (no outer cells/edges found).");
    lua_error(lua);
    return LUA_ERRRUN;
  }

  // Prune the mesh.
  prune_voronoi_mesh(mesh);

  // No results.
  return 0;
}

extern void interpreter_register_spfuncs(interpreter_t* interp);

static int jostle_points(lua_State* lua)
{
  // Check the arguments.
  int num_args = lua_gettop(lua);
  if ((num_args != 3) || !lua_ispointlist(lua, 1) || 
      !lua_isnumber(lua, 2) || !lua_isnumber(lua, 3))
  {
    lua_pushstring(lua, "Invalid argument(s). Usage:\n"
                        "jostle_points(points, radius, factor)\n"
                        "-> jostles points within a given radius using a\n"
                        "given randomness factor.");
    lua_error(lua);
    return LUA_ERRRUN;
  }

  // Get the points.
  int num_points;
  point_t* points = lua_topointlist(lua, 1, &num_points);

  // Get the jostling radius.
  double radius = lua_tonumber(lua, 2);
  if (radius < 0.0)
  {
    lua_pushstring(lua, "Jostling radius must be non-negative.");
    lua_error(lua);
    return LUA_ERRRUN;
  }

  // Get the randomness factor.
  double randomness = lua_tonumber(lua, 5);
  if ((randomness < 0.0) || (randomness > 1.0))
  {
    lua_pushstring(lua, "Fifth argument must be a randomness factor between 0 and 1.");
    lua_error(lua);
    return LUA_ERRRUN;
  }
  
  for (int i = 0; i < num_points; ++i)
  {
    point_t* x = &points[i];
    x->x += (randomness*random()/RAND_MAX - 0.5)*radius;
    x->y += (randomness*random()/RAND_MAX - 0.5)*radius;
    x->z += (randomness*random()/RAND_MAX - 0.5)*radius;
  }

  // No results.
  return 0;
}

void interpreter_register_geometry_functions(interpreter_t* interp)
{
  interpreter_register_function(interp, "cubic_lattice_mesh", cubic_lattice_mesh);
  interpreter_register_function(interp, "cubic_lattice_periodic_bc", cubic_lattice_periodic_bc);
  interpreter_register_function(interp, "random_points", random_points);
  interpreter_register_function(interp, "ccp_points", ccp_points);
  interpreter_register_function(interp, "jostle_points", jostle_points);
  interpreter_register_function(interp, "unbounded_voronoi_mesh", unbounded_voronoi_mesh);
  interpreter_register_function(interp, "bounded_voronoi_mesh", bounded_voronoi_mesh);
  interpreter_register_function(interp, "prune_voronoi_mesh", prune_voronoi_mesh_);
  interpreter_register_function(interp, "bound_voronoi_mesh", bound_voronoi_mesh_);
  interpreter_register_spfuncs(interp);
}

