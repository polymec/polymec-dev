// Copyright (c) 2012-2017, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/kd_tree.h"
#include "core/slist.h"
#include "core/unordered_set.h"
#include "geometry/rect_prism_sp_func.h"
#include "model/interpreter.h"

// Lua stuff.
#include "lua.h"
#include "lualib.h"
#include "lauxlib.h"

// Various function libraries.
extern void interpreter_register_sd_functions(interpreter_t* interp);
extern void interpreter_register_scalar_functions(interpreter_t* interp);
extern void interpreter_register_vector_functions(interpreter_t* interp);
extern void interpreter_register_indicators(interpreter_t* interp);

// Functions for the point factory, which manufactures sets of points.
extern int point_factory_random_points(lua_State* lua);
extern docstring_t* point_factory_random_points_doc();
extern int point_factory_cubic_lattice(lua_State* lua);
extern docstring_t* point_factory_cubic_lattice_doc();
extern int point_factory_cylinder(lua_State* lua);
extern docstring_t* point_factory_cylinder_doc();
extern int point_factory_import_from_cad(lua_State* lua);
extern docstring_t* point_factory_import_from_cad_doc();

extern int mesh_factory_uniform(lua_State* lua);
extern docstring_t* mesh_factory_uniform_doc();
extern int mesh_factory_rectilinear(lua_State* lua);
extern docstring_t* mesh_factory_rectilinear_doc();

extern int cartesian_1d_uniform(lua_State* lua);
extern docstring_t* cartesian_1d_uniform_doc();
extern int cartesian_1d_logarithmic(lua_State* lua);
extern docstring_t* cartesian_1d_logarithmic_doc();
extern int cartesian_1d_nonuniform(lua_State* lua);
extern docstring_t* cartesian_1d_nonuniform_doc();

extern int cylindrical_1d_uniform(lua_State* lua);
extern docstring_t* cylindrical_1d_uniform_doc();
extern int cylindrical_1d_logarithmic(lua_State* lua);
extern docstring_t* cylindrical_1d_logarithmic_doc();
extern int cylindrical_1d_nonuniform(lua_State* lua);
extern docstring_t* cylindrical_1d_nonuniform_doc();

extern int spherical_1d_uniform(lua_State* lua);
extern docstring_t* spherical_1d_uniform_doc();
extern int spherical_1d_logarithmic(lua_State* lua);
extern docstring_t* spherical_1d_logarithmic_doc();
extern int spherical_1d_nonuniform(lua_State* lua);
extern docstring_t* spherical_1d_nonuniform_doc();

static int sample_bbox(lua_State* lua)
{
  // Check the arguments.
  int num_args = lua_gettop(lua);
  if ((num_args != 4) || !lua_isboundingbox(lua, 1) || 
      !lua_isnumber(lua, 2) || !lua_isnumber(lua, 3) ||
      !lua_isnumber(lua, 4))
  {
    return luaL_error(lua, "Invalid argument(s). Usage:\n"
                      "points = sample_bounding_box(bbox, nx, ny, nz)");
  }

  bbox_t* bbox = lua_toboundingbox(lua, 1);
  int nx = (int)(lua_tonumber(lua, 2));
  if (nx <= 0)
    return luaL_error(lua, "nx must be a positive number of x points.");

  int ny = (int)(lua_tonumber(lua, 3));
  if (ny <= 0)
    return luaL_error(lua, "ny must be a positive number of y points.");

  int nz = (int)(lua_tonumber(lua, 4));
  if (nz <= 0)
    return luaL_error(lua, "nz must be a positive number of z points.");

  int num_points = 2*(nx*ny + ny*nz + nz*nx);
  real_t dx = (bbox->x2 - bbox->x1) / nx;
  real_t dy = (bbox->y2 - bbox->y1) / ny;
  real_t dz = (bbox->z2 - bbox->z1) / nz;
  point_t* points = polymec_malloc(sizeof(point_t) * num_points);
  int offset = 0;

  // -x face.
  for (int i = 0; i < ny; ++i)
  {
    for (int j = 0; j < nz; ++j, ++offset)
    {
      points[offset].x = bbox->x1;
      points[offset].y = bbox->y1 + (i+0.5) * dy;
      points[offset].z = bbox->z1 + (j+0.5) * dz;
    }
  }

  // +x face.
  for (int i = 0; i < ny; ++i)
  {
    for (int j = 0; j < nz; ++j, ++offset)
    {
      points[offset].x = bbox->x2;
      points[offset].y = bbox->y1 + (i+0.5) * dy;
      points[offset].z = bbox->z1 + (j+0.5) * dz;
    }
  }

  // -y face.
  for (int i = 0; i < nx; ++i)
  {
    for (int j = 0; j < nz; ++j, ++offset)
    {
      points[offset].x = bbox->x1 + (i+0.5) * dx;
      points[offset].y = bbox->y1;
      points[offset].z = bbox->z1 + (j+0.5) * dz;
    }
  }

  // +y face.
  for (int i = 0; i < nx; ++i)
  {
    for (int j = 0; j < nz; ++j, ++offset)
    {
      points[offset].x = bbox->x1 + (i+0.5) * dx;
      points[offset].y = bbox->y2;
      points[offset].z = bbox->z1 + (j+0.5) * dz;
    }
  }

  // -z face.
  for (int i = 0; i < nx; ++i)
  {
    for (int j = 0; j < ny; ++j, ++offset)
    {
      points[offset].x = bbox->x1 + (i+0.5) * dx;
      points[offset].y = bbox->y1 + (j+0.5) * dy;
      points[offset].z = bbox->z1;
    }
  }

  // +z face.
  for (int i = 0; i < nx; ++i)
  {
    for (int j = 0; j < ny; ++j, ++offset)
    {
      points[offset].x = bbox->x1 + (i+0.5) * dx;
      points[offset].y = bbox->y1 + (j+0.5) * dy;
      points[offset].z = bbox->z2;
    }
  }

  lua_pushpointlist(lua, points, num_points);
  return 1;
}

static docstring_t* sample_bbox_doc()
{
  return docstring_from_string("sample_bounding_box(bbox, nx, ny, nz) -\n"
                               "  Returns a set points on a lattice that covers a bounding box.");
}

static int scaled_bounding_box(lua_State* lua)
{
  // Check the arguments.
  int num_args = lua_gettop(lua);
  if ((num_args != 2) || !lua_isboundingbox(lua, 1) || 
      !lua_isnumber(lua, 2))
  {
    return luaL_error(lua, "Invalid argument(s). Usage:\n"
                      "bbox2 = scaled_bounding_box(bbox, factor)");
  }

  bbox_t* bbox = lua_toboundingbox(lua, 1);
  real_t f = (real_t)lua_tonumber(lua, 2);
  if (f <= 0.0)
    return luaL_error(lua, "factor must be positive.");

  real_t xc = 0.5 * (bbox->x1 + bbox->x2);
  real_t yc = 0.5 * (bbox->y1 + bbox->y2);
  real_t zc = 0.5 * (bbox->z1 + bbox->z2);
  bbox_t* scaled_bbox = bbox_new((bbox->x1 - xc) * f + xc, 
                                 (bbox->x2 - xc) * f + xc,
                                 (bbox->y1 - yc) * f + yc,
                                 (bbox->y2 - yc) * f + yc,
                                 (bbox->z1 - zc) * f + zc,
                                 (bbox->z2 - zc) * f + zc);
  lua_pushboundingbox(lua, scaled_bbox);
  return 1;
}

static docstring_t* scaled_bounding_box_doc()
{
  return docstring_from_string("scaled_bounding_box(bbox, factor) -\n"
                               "  Returns a bounding box scaled by the given factor.");
}

static int sample_cyl_shell(lua_State* lua)
{
  // Check the arguments.
  int num_args = lua_gettop(lua);
  if ((num_args != 7) || 
      !lua_isnumber(lua, 1) || !lua_isnumber(lua, 2) || 
      !lua_isnumber(lua, 3) || !lua_isnumber(lua, 4) || 
      !lua_isnumber(lua, 5) || !lua_isnumber(lua, 6) || 
      !lua_isnumber(lua, 7))
  {
    return luaL_error(lua, "Invalid argument(s). Usage:\n"
                      "points = sample_cyl_shell(r1, r2, z1, z2, nr, nphi, nz)");
  }

  real_t r1 = (real_t)lua_tonumber(lua, 1);
  if (r1 < 0.0)
    return luaL_error(lua, "r1 must be a non-negative inner radius.");

  real_t r2 = (real_t)lua_tonumber(lua, 2);
  if (r2 <= r1)
    return luaL_error(lua, "r2 must be greater than r1.");

  real_t z1 = (real_t)lua_tonumber(lua, 3);
  real_t z2 = (real_t)lua_tonumber(lua, 4);
  if (z2 <= z1)
    return luaL_error(lua, "z2 must be greater than z1.");

  int nr = (int)(lua_tonumber(lua, 5));
  if (nr <= 0)
    return luaL_error(lua, "nr must be a positive number of radial points.");

  int nphi = (int)(lua_tonumber(lua, 6));
  if (nphi <= 0)
    return luaL_error(lua, "nphi must be a positive number of azimuthal points.");

  int nz = (int)(lua_tonumber(lua, 7));
  if (nz <= 0)
    return luaL_error(lua, "nz must be a positive number of axial points.");

  int num_points = nr * nphi * nz;
  real_t dr = (r2 - r1) / nr;
  real_t dphi = 2.0 * M_PI / nphi;
  real_t dz = (z2 - z1) / nz;
  point_t* points = polymec_malloc(sizeof(point_t) * num_points);
  int offset = 0;

  for (int i = 0; i < nr; ++i)
  {
    real_t r = r1 + (i+0.5) * dr;
    for (int j = 0; j < nphi; ++j)
    {
      real_t phi = (j+0.5) * dphi;
      for (int k = 0; k < nz; ++k, ++offset)
      {
        real_t z = z1 + (k+0.5) * dz;
        points[offset].x = r*cos(phi);
        points[offset].y = r*sin(phi);
        points[offset].z = z;
      }
    }
  }
  lua_pushpointlist(lua, points, num_points);
  return 1;
}

static docstring_t* sample_cyl_shell_doc()
{
  return docstring_from_string("sample_cyl_shell(r1, r2, z1, z2, nr, nphi, nz) - \n"
                               "  Returns a set points on a lattice that covers a cylindrical shell.");
}

static int translate_points(lua_State* lua)
{
  // Check the arguments.
  int num_args = lua_gettop(lua);
  if (((num_args != 2) && (num_args != 3)) ||
      ((num_args == 2) && (!lua_ispointlist(lua, 1) || (!lua_isvector(lua, 2) && !lua_isvectorlist(lua, 2)))) || 
      ((num_args == 3) && (!lua_ispointlist(lua, 1) || (!lua_isvector(lua, 2) && !lua_isvectorlist(lua, 2)) || (!lua_isnumber(lua, 3) && !lua_issequence(lua, 3)))))
  {
    return luaL_error(lua, "Invalid argument(s). Usage:\n"
                      "translate_points(points, vector) OR\n"
                      "translate_points(points, vector, factor) OR\n"
                      "translate_points(points, vectors) OR\n"
                      "translate_points(points, vectors, factor) OR\n"
                      "translate_points(points, vectors, factors)");
  }

  int num_points;
  point_t* points = lua_topointlist(lua, 1, &num_points);

  real_t factor = 1.0;
  real_t *factors = NULL;
  int num_factors = -1;

  if (num_args == 3)
  {
    if (lua_isnumber(lua, 3))
      factor = (real_t)lua_tonumber(lua, 3);
    else
    {
      factors = lua_tosequence(lua, 3, &num_factors);
      if (num_factors != num_points)
        return luaL_error(lua, "Number of scale factors must equal number of points.");
    }
  }

  if (lua_isvector(lua, 2))
  {
    vector_t* vector = lua_tovector(lua, 2);
    for (int i = 0; i < num_points; ++i)
    {
      real_t f = (factors != NULL) ? factors[i] : factor;
      points[i].x += f * vector->x;
      points[i].y += f * vector->y;
      points[i].z += f * vector->z;
    }
  }
  else
  {
    int num_vectors;
    vector_t* vectors = lua_tovectorlist(lua, 2, &num_vectors);
    if (num_vectors != num_points)
      return luaL_error(lua, "Number of vectors must equal number of points.");
    for (int i = 0; i < num_points; ++i)
    {
      real_t f = (factors != NULL) ? factors[i] : factor;
      points[i].x += f * vectors[i].x;
      points[i].y += f * vectors[i].y;
      points[i].z += f * vectors[i].z;
    }
  }

  // Modified in place.
  return 0;
}

static docstring_t* translate_points_doc()
{
  return docstring_from_string("translate_points(points, vector) OR\n"
                               "translate_points(points, vector, factor) OR\n"
                               "translate_points(points, vectors) OR\n"
                               "translate_points(points, vectors, factor) OR\n"
                               "translate_points(points, vectors, factors) -\n"
                               "  Translates a set of points by the given constant vector or corresponding vectors.");
}

static int rotate_points(lua_State* lua)
{
  // Check the arguments.
  int num_args = lua_gettop(lua);
  if ((num_args < 2) ||
      !lua_ispointlist(lua, 1) || 
      !lua_isnumber(lua, 2) || 
      ((num_args >= 3) && !lua_isvector(lua, 3)) ||
      ((num_args == 4) && !lua_ispoint(lua, 4)))
  {
    return luaL_error(lua, "Invalid argument(s). Usage:\n"
                      "rotate_points(points, angle, axis = {0, 0, 1}, origin = {0, 0, 0})");
  }

  int num_points;
  point_t* points = lua_topointlist(lua, 1, &num_points);
  real_t angle = (real_t)lua_tonumber(lua, 2);
  vector_t* axis = (num_args >= 3) ? lua_tovector(lua, 3) : vector_new(0.0, 0.0, 1.0);
  point_t* origin = (num_args == 4) ? lua_topoint(lua, 4) : point_new(0.0, 0.0, 0.0);

  for (int i = 0; i < num_points; ++i)
  {
    // Relative coordinates: y = x - origin.
    vector_t y;
    point_displacement(origin, &points[i], &y);

    // Set up an orthonormal basis for the given axis.
    vector_t e1, e2;
    compute_orthonormal_basis(axis, &e1, &e2);

    // Project y onto the e1 x e2 plane.
    real_t xi1 = vector_dot(&y, &e1);
    real_t xi2 = vector_dot(&y, &e2);
    real_t y3 = vector_dot(&y, axis);

    // Rotate it about the axis by the angle.
    real_t Rxi1 = xi1*cos(angle) - xi2*sin(angle);
    real_t Rxi2 = xi1*sin(angle) + xi2*cos(angle);

    // Back to original coordinate frame.
    points[i].x = origin->x + Rxi1*e1.x + Rxi2*e2.x + y3*axis->x;
    points[i].y = origin->y + Rxi1*e1.y + Rxi2*e2.y + y3*axis->y;
    points[i].z = origin->z + Rxi1*e1.z + Rxi2*e2.z + y3*axis->z;
  }

  // Modified in place.
  return 0;
}

static docstring_t* rotate_points_doc()
{
  return docstring_from_string("rotate_points(points, angle, axis = {0, 0, 1}, origin = {0, 0, 0}) -\n"
                               "  Rotates a set points about the axis by the given angle.");
}

static int copy_points(lua_State* lua)
{
  // Check the arguments.
  int num_args = lua_gettop(lua);
  if ((num_args != 1) || !lua_ispointlist(lua, 1))
  {
    return luaL_error(lua, "Invalid argument(s). Usage:\n"
                      "new_points = copy_points(points)");
  }

  int num_points;
  point_t* old_points = lua_topointlist(lua, 1, &num_points);
  point_t* new_points = polymec_malloc(sizeof(point_t) * num_points);
  memcpy(new_points, old_points, sizeof(point_t) * num_points);
  lua_pushpointlist(lua, new_points, num_points);

  return 1;
}

static docstring_t* copy_points_doc()
{
  return docstring_from_string("copy_points(points) -\n"
                               "  Creates a new copy of a list of points.");
}

static int select_points(lua_State* lua)
{
  // Check the arguments.
  int num_args = lua_gettop(lua);
  if ((num_args != 2) || !lua_ispointlist(lua, 1) || !lua_istable(lua, 2))
  {
    return luaL_error(lua, "Invalid argument(s). Usage:\n"
                      "selected_points = select_points(points, options)");
  }

  // Extract arguments.
  int num_points = 0;
  point_t* points = lua_topointlist(lua, 1, &num_points);
  int num_near_points = 0;
  point_t* near_points = NULL;
  real_t within_distance = 0.0;
  st_func_t* within_surface = NULL;
  bbox_t* within_bbox = NULL;
  real_t at_time = 0.0;
  const char* entries[] = {"near_points", "within_distance", "within_surface", "at_time"};
  for (int i = 0; i < 4; ++i)
  {
    lua_pushstring(lua, entries[i]);
    lua_gettable(lua, 2);
    if (i == 0) // near_points
    {
      if (lua_isnil(lua, -1))
      {
        lua_pop(lua, 1);
        lua_pushinteger(lua, 1);
        lua_gettable(lua, 2);
      }
      if (!lua_ispointlist(lua, -1))
        return luaL_error(lua, "near_points should be a list of points.");

      near_points = lua_topointlist(lua, -1, &num_near_points);
    }
    else if (i == 1) // within_distance
    {
      if (!lua_isnumber(lua, -1))
        return luaL_error(lua, "within_distance should be a number.");
      within_distance = (real_t)lua_tonumber(lua, -1);
      if (within_distance <= 0.0)
        return luaL_error(lua, "within_distance should be positive.");
    }
    else if (i == 2) // within_surface
    {
      if (lua_isboundingbox(lua, -1))
      {
        within_bbox = lua_toboundingbox(lua, -1);
        sp_func_t* bbox_surface = rect_prism_sp_func_from_bbox(within_bbox);
        within_surface = st_func_from_sp_func(bbox_surface);
        bbox_surface = NULL;
      }
      else if (lua_isscalarfunction(lua, -1))
        within_surface = lua_toscalarfunction(lua, -1);
      else if (!lua_isnil(lua, -1))
        return luaL_error(lua, "within_surface should be a scalar function.");
    }
    else // at_time
    {
      if (lua_isnil(lua, -1)) break;
      if (!lua_isnumber(lua, -1))
        return luaL_error(lua, "at_time should be a number.");

      at_time = (real_t)lua_tonumber(lua, -1);
    }
  }

  // Take stock of our arguments.
  if (near_points != NULL) 
  {
    if (reals_equal(within_distance, 0.0))
      return luaL_error(lua, "within_distance should be given with near_points.");
    if (within_surface != NULL)
      return luaL_error(lua, "within_surface cannot be specified with near_points.");
  }
  if ((near_points == NULL) && (within_surface == NULL))
    return luaL_error(lua, "no criteria for selecting points!");

  // Now go get the points.
  int_slist_t* selected_points = int_slist_new();
  if (near_points != NULL)
  {
    // Stick our points in a kd-tree.
    kd_tree_t* tree = kd_tree_new(near_points, num_near_points);
    for (int i = 0; i < num_points; ++i)
    {
      int j = kd_tree_nearest(tree, &points[i]);
      if (point_distance(&points[i], &near_points[j]) < within_distance)
        int_slist_append(selected_points, i);
    }

    // Clean up.
    kd_tree_free(tree);
  }
  else 
  {
    ASSERT(within_surface != NULL);
    for (int i = 0; i < num_points; ++i)
    {
      real_t F;
      st_func_eval(within_surface, &points[i], at_time, &F);
      if (F >= 0.0)
        int_slist_append(selected_points, i);
    }
  }
  ASSERT(selected_points->size > 0);

  // Construct a sequence consisting of the point indices.
  // NOTE: we use real_t because that's how lua represents numbers.
  real_t* s_points = polymec_malloc(sizeof(real_t) * selected_points->size);

  int i = 0;
  for (int_slist_node_t* node = selected_points->front; node != NULL; node = node->next, ++i)
    s_points[i] = (real_t)node->value;
  lua_pushsequence(lua, s_points, selected_points->size);

  // Clean up.
  within_surface = NULL;
  int_slist_free(selected_points);

  return 1;
}

static docstring_t* select_points_doc()
{
  return docstring_from_string("select_points(points, options) -\n"
                               "  Returns a list of indices of points from the given list that satisfy certain criteria.\n"
                               "  Criteria are specified by the options table, which contains one or more of the following\n"
                               "  entries:\n"
                               "  - near_points = [list of points] <-- selects points near those in this list\n" 
                               "  - within_distance = D  <-- should be specified with near_points\n" 
                               "  - within_surface = F   <-- selects points inside of the surface defined by the function F\n"
                               "  - at_time = t <-- should be given with within_surface if F accepts a time");
}

static int remove_points(lua_State* lua)
{
  // Check the arguments.
  int num_args = lua_gettop(lua);
  if ((num_args != 2) || !lua_ispointlist(lua, 1) || !lua_issequence(lua, 2))
  {
    return luaL_error(lua, "Invalid argument(s). Usage:\n"
                      "trimmed_points = remove_points(points, selected_indices)");
  }

  // Extract arguments.
  int num_points = 0, num_selected_points = 0;
  point_t* points = lua_topointlist(lua, 1, &num_points);
  real_t* selected_points = lua_tosequence(lua, 2, &num_selected_points);

  if (num_selected_points > num_points)
    return luaL_error(lua, "remove_points: Cannot remove %d points from a set of %d.", num_selected_points, num_points);

  // Now remove the points.
  int_unordered_set_t* removed_points = int_unordered_set_new();
  for (int i = 0; i < num_selected_points; ++i)
    int_unordered_set_insert(removed_points, (int)selected_points[i]);

  // Construct a copy of the list of points, with the removed points absent.
  point_t* trimmed_points = polymec_malloc(sizeof(point_t) * num_points - num_selected_points);
  int j = 0;
  for (int i = 0; i < num_points; ++i)
  {
    if (!int_unordered_set_contains(removed_points, i))
      trimmed_points[j++] = points[i];
  }
  lua_pushpointlist(lua, trimmed_points, num_points - num_selected_points);

  // Clean up.
  int_unordered_set_free(removed_points);

  return 1;
}

static docstring_t* remove_points_doc()
{
  return docstring_from_string("remove_points(points, selected_indices) -\n"
                               "  Returns a copy of points with the given selected points (identified by\n"
                               "  a sequence of indices) removed.");
}

void interpreter_register_geometry_functions(interpreter_t* interp);
void interpreter_register_geometry_functions(interpreter_t* interp)
{
  // Set up a point factory for generating points in 3D.
  interpreter_register_global_table(interp, "point_factory", docstring_from_string("point_factory - Functions for generating 3D point distributions.")); 
  interpreter_register_global_method(interp, "point_factory", "random_points", point_factory_random_points, point_factory_random_points_doc());
  interpreter_register_global_method(interp, "point_factory", "cubic_lattice", point_factory_cubic_lattice, point_factory_cubic_lattice_doc());
  interpreter_register_global_method(interp, "point_factory", "cylinder", point_factory_cylinder, point_factory_cylinder_doc());
  interpreter_register_global_method(interp, "point_factory", "import_from_cad", point_factory_import_from_cad, point_factory_import_from_cad_doc());

  // Set up a mesh factory for generating 3D meshes.
  interpreter_register_global_table(interp, "mesh_factory", docstring_from_string("mesh_factory - Functions for generating 3D polyhedral meshes.")); 
  interpreter_register_global_method(interp, "mesh_factory", "uniform", mesh_factory_uniform, mesh_factory_uniform_doc());
  interpreter_register_global_method(interp, "mesh_factory", "rectilinear", mesh_factory_rectilinear, mesh_factory_rectilinear_doc());

  // Set up a factory for generating 1D Cartesian meshes
  interpreter_register_global_table(interp, "cartesian_1d", docstring_from_string("cartesian_1d - Functions for creating 1D point distributions."));
  interpreter_register_global_method(interp, "cartesian_1d", "uniform", cartesian_1d_uniform, cartesian_1d_uniform_doc());
  interpreter_register_global_method(interp, "cartesian_1d", "logarithmic", cartesian_1d_logarithmic, cartesian_1d_logarithmic_doc());
  interpreter_register_global_method(interp, "cartesian_1d", "nonuniform", cartesian_1d_nonuniform, cartesian_1d_nonuniform_doc());

  // Set up a factory for generating 1D cylindrical meshes
  interpreter_register_global_table(interp, "cylindrical_1d", docstring_from_string("cylindrical_1d - Functions for creating 1D cylindrically-symmetric points."));
  interpreter_register_global_method(interp, "cylindrical_1d", "uniform", cylindrical_1d_uniform, cylindrical_1d_uniform_doc());
  interpreter_register_global_method(interp, "cylindrical_1d", "logarithmic", cylindrical_1d_logarithmic, cylindrical_1d_logarithmic_doc());
  interpreter_register_global_method(interp, "cylindrical_1d", "nonuniform", cylindrical_1d_nonuniform, cylindrical_1d_nonuniform_doc());

  // Set up a factory for generating 1D spherical meshes
  interpreter_register_global_table(interp, "spherical_1d", docstring_from_string("spherical_1d - Functions for creating 1D spherically-symmetric points."));
  interpreter_register_global_method(interp, "spherical_1d", "uniform", spherical_1d_uniform, spherical_1d_uniform_doc());
  interpreter_register_global_method(interp, "spherical_1d", "logarithmic", spherical_1d_logarithmic, spherical_1d_logarithmic_doc());
  interpreter_register_global_method(interp, "spherical_1d", "nonuniform", spherical_1d_nonuniform, spherical_1d_nonuniform_doc());

  interpreter_register_function(interp, "scaled_bounding_box", scaled_bounding_box, scaled_bounding_box_doc());
  interpreter_register_function(interp, "sample_bounding_box", sample_bbox, sample_bbox_doc());
  interpreter_register_function(interp, "sample_cyl_shell", sample_cyl_shell, sample_cyl_shell_doc());
  interpreter_register_function(interp, "translate_points", translate_points, translate_points_doc());
  interpreter_register_function(interp, "rotate_points", rotate_points, rotate_points_doc());
  interpreter_register_function(interp, "copy_points", copy_points, copy_points_doc());
  interpreter_register_function(interp, "select_points", select_points, select_points_doc());
  interpreter_register_function(interp, "remove_points", remove_points, remove_points_doc());

  interpreter_register_sd_functions(interp);
  interpreter_register_scalar_functions(interp);
  interpreter_register_vector_functions(interp);
  interpreter_register_indicators(interp);
}

