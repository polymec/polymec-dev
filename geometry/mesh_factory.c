// Copyright (c) 2012-2014, Jeffrey N. Johnson
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without 
// modification, are permitted provided that the following conditions are met:
// 
// 1. Redistributions of source code must retain the above copyright notice, this 
// list of conditions and the following disclaimer.
// 
// 2. Redistributions in binary form must reproduce the above copyright notice, 
// this list of conditions and the following disclaimer in the documentation 
// and/or other materials provided with the distribution.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE 
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

// mesh_factory.c - Implementations of interpreter functions for generating
// meshes.

#include <strings.h>
#include "core/polymec.h"
#include "core/array.h"
#include "core/mesh.h"
#include "core/kd_tree.h"
#include "core/interpreter.h"
#include "core/unordered_map.h"
#include "core/tuple.h"
#include "geometry/create_uniform_mesh.h"
#include "geometry/create_rectilinear_mesh.h"
#include "geometry/create_pebi_mesh.h"
#include "geometry/create_boundary_generators.h"
#include "geometry/rect_prism.h"

#if POLYMEC_HAVE_TETGEN
#include "geometry/create_voronoi_mesh.h"
//#include "geometry/create_cvt_with_lloyd_iteration.h"
#endif

// Lua stuff.
#include "lua.h"
#include "lualib.h"
#include "lauxlib.h"

int mesh_factory_uniform(lua_State* lua)
{
  // Check the arguments.
  int num_args = lua_gettop(lua);
  if (num_args != 4)
  {
    return luaL_error(lua, "Invalid arguments. Usage:\n"
                      "mesh = mesh_factory.uniform(nx, ny, nz, bounds)");
  }

  // Get the arguments.
  int nx = (int)lua_tonumber(lua, 1);
  int ny = (int)lua_tonumber(lua, 2);
  int nz = (int)lua_tonumber(lua, 3);
  if ((nx <= 0) || (ny <= 0) || (nz <= 0))
    return luaL_error(lua, "nx, ny, and nz must all be positive.");
  if (!lua_isboundingbox(lua, 4))
    return luaL_error(lua, "bounds must be a bounding box.");

  // Bounding box.
  bbox_t* bbox = lua_toboundingbox(lua, 4);

  // Pop all the previous arguments off the stack.
  lua_pop(lua, lua_gettop(lua));

  // Create the mesh.
  mesh_t* mesh = create_uniform_mesh(MPI_COMM_WORLD, nx, ny, nz, bbox);

  // Tag its faces.
  tag_rectilinear_mesh_faces(mesh, nx, ny, nz, "x1", "x2", "y1", "y2", "z1", "z2");

  // Push the mesh onto the stack.
  lua_pushmesh(lua, mesh);
  return 1;
}

int mesh_factory_rectilinear(lua_State* lua)
{
  // Check the arguments.
  int num_args = lua_gettop(lua);
  if ((num_args != 3) || !lua_issequence(lua, 1) || !lua_issequence(lua, 2) || !lua_issequence(lua, 3))
  {
    return luaL_error(lua, "Invalid arguments. Usage:\n"
                      "mesh = mesh_factory.rectilinear(xs, ys, zs)");
  }

  // Get the arguments.
  int nxs, nys, nzs;
  real_t* xs = lua_tosequence(lua, 1, &nxs);
  real_t* ys = lua_tosequence(lua, 2, &nys);
  real_t* zs = lua_tosequence(lua, 3, &nzs);
  if ((nxs < 2) || (nys < 2) || (nzs < 2))
  {
    return luaL_error(lua, "xs, ys, and zs must all contain at least 2 values.");
  }

  // Pop all the previous arguments off the stack.
  lua_pop(lua, lua_gettop(lua));

  // Create the mesh.
  mesh_t* mesh = create_rectilinear_mesh(MPI_COMM_WORLD, xs, nxs, ys, nys, zs, nzs);

  // Tag its faces.
  tag_rectilinear_mesh_faces(mesh, nxs-1, nys-1, nzs-1, "x1", "x2", "y1", "y2", "z1", "z2");

  // Push the mesh onto the stack.
  lua_pushmesh(lua, mesh);
  return 1;
}

#if 0
int mesh_factory_cubic_lattice_periodic_bc(lua_State* lua)
{
  // Check the argument.
  int num_args = lua_gettop(lua);
  if (num_args != 2)
    return luaL_error(lua, "Arguments must be 2 boundary mesh (face) tags.");

  for (int i = 1; i <= 2; ++i)
  {
    if (!lua_isstring(lua, i))
      return luaL_error(lua, "Argument %d must be a face tag.", i);
  }

  const char* tag1 = lua_tostring(lua, 1);
  const char* tag2 = lua_tostring(lua, 2);

  // Based on the names of the tags, we can decide which periodic BC 
  // mapping function to use.
  periodic_bc_t* bc = NULL;
  if (!strcmp(tag1, "x1"))
  {
    if (!strcmp(tag2, "x2"))
      return luaL_error(lua, "Periodic boundary maps from 'x1' to '%s' (must be 'x2').", tag2);

    bc = cubic_lattice_x_periodic_bc_new(tag1, tag2);
  }
  else if (!strcmp(tag1, "x2"))
  {
    if (!strcmp(tag2, "x1"))
      return luaL_error(lua, "Periodic boundary maps from 'x2' to '%s' (must be 'x1').", tag2);

    bc = cubic_lattice_x_periodic_bc_new(tag1, tag2);
  }
  else if (!strcmp(tag1, "y1"))
  {
    if (!strcmp(tag2, "y2"))
      return luaL_error(lua, "Periodic boundary maps from 'y1' to '%s' (must be 'y2').", tag2);

    bc = cubic_lattice_y_periodic_bc_new(tag1, tag2);
  }
  else if (!strcmp(tag1, "y2"))
  {
    if (!strcmp(tag2, "y1"))
      return luaL_error(lua, "Periodic boundary maps from 'y2' to '%s' (must be 'y1').", tag2);

    bc = cubic_lattice_y_periodic_bc_new(tag1, tag2);
  }
  else if (!strcmp(tag1, "z1"))
  {
    if (!strcmp(tag2, "z2"))
      return luaL_error(lua, "Periodic boundary maps from 'z1' to '%s' (must be 'z2').", tag2);

    bc = cubic_lattice_z_periodic_bc_new(tag1, tag2);
  }
  else if (!strcmp(tag1, "z2"))
  {
    if (!strcmp(tag2, "z1"))
      return luaL_error(lua, "Periodic boundary maps from 'z2' to '%s' (must be 'z1').", tag2);

    bc = cubic_lattice_z_periodic_bc_new(tag1, tag2);
  }
  lua_pushuserdefined(lua, bc, NULL);
  return 1;
}
#endif

static void free_string(char* str)
{
  free(str);
}

#if POLYMEC_HAVE_TETGEN
int mesh_factory_voronoi(lua_State* lua)
{
  // Check the arguments.
  int num_args = lua_gettop(lua);
  if (((num_args == 1) && !lua_ispointlist(lua, 1) && !lua_istable(lua, 1)) || 
      ((num_args == 2) && (!lua_istable(lua, 1) && !lua_isboundingbox(lua, 2))))
  {
    return luaL_error(lua, "Invalid argument(s). Usage:\n"
                      "mesh = mesh_factory.voronoi(generators) OR\n"
                      "mesh = mesh_factory.voronoi(generators, bounding_box)\n\n"
                      "Above, generators is a list of points or a table mapping\n"
                      "tag names to points.");
  }

  // Get the generators.
  int num_generators = 0;
  point_t* generators = NULL;
  string_ptr_unordered_map_t* tag_map = NULL;
  if (lua_ispointlist(lua, 1))
    generators = lua_topointlist(lua, 1, &num_generators);
  else
  {
    if (!lua_istable(lua, 1))
      return luaL_error(lua, "Argument 1 must be a list of points or a table mapping tag names to points.");

    lua_pushnil(lua);
    tag_map = string_ptr_unordered_map_new();
    while (lua_next(lua, 1))
    {
      static const int key_index = -2;
      static const int val_index = -1;
      const char* key = lua_tostring(lua, key_index);
      if (!lua_ispointlist(lua, val_index))
        return luaL_error(lua, "Table values must be lists of points.");

      // Extract the points associated with this tag.
      int Np;
      point_t* p = lua_topointlist(lua, val_index, &Np);
      num_generators += Np;
      generators = realloc(generators, sizeof(point_t) * num_generators);
      memcpy(&generators[num_generators-Np], p, sizeof(point_t) * Np);

      // Write down the generator offset and size associated with this tag.
      int* tuple = int_tuple_new(2);
      tuple[0] = num_generators - Np;
      tuple[1] = Np;
      string_ptr_unordered_map_insert_with_v_dtor(tag_map, (char*)key, tuple, DTOR(int_tuple_free));

      lua_pop(lua, 1);
    }
  }
  bbox_t* bbox = NULL;
  if (num_args == 2)
    bbox = lua_toboundingbox(lua, 2);

  // Create the mesh.
  mesh_t* mesh = NULL;
  if (bbox != NULL)
  {
    // Make sure all the points are within the bounding box.
    sp_func_t* boxy = rect_prism_new_from_bbox(bbox);
    for (int i = 0; i < num_generators; ++i)
    {
      real_t F = 0.0;
      sp_func_eval(boxy, &generators[i], &F);
      if (F > 0.0)
        return luaL_error(lua, "mesh_factory.voronoi: Generator %d is outside the given bounding box.", i);
    }
    boxy = NULL;
    mesh = create_voronoi_mesh_in_box(MPI_COMM_WORLD, generators, num_generators, NULL, 0, bbox);
  }
  else
  {
    mesh = create_voronoi_mesh(MPI_COMM_WORLD, generators, num_generators, NULL, 0);
  }
  ASSERT(mesh != NULL);

  // Tag the generators if necessary.
  if (tag_map != NULL)
  {
    int pos = 0, *tuple;
    char* tag_name;
    while (string_ptr_unordered_map_next(tag_map, &pos, &tag_name, (void*)(&tuple)))
    {
      int tag_offset = tuple[0];
      int tag_size = tuple[1];
      int* tag = mesh_create_tag(mesh->cell_tags, tag_name, tag_size);
      for (int i = 0; i < tag_size; ++i)
        tag[i] = tag_offset + i;
    }

    string_ptr_unordered_map_free(tag_map);
  }

  // Push the mesh onto the stack.
  lua_pushmesh(lua, mesh);
  return 1;
}

int mesh_factory_cvt(lua_State* lua)
{
return luaL_error(lua, "CURRENTLY NOT SUPPORTED.");
  // Check the arguments.
  int num_args = lua_gettop(lua);
  if ((num_args != 2) || !lua_istable(lua, 1) || !lua_istable(lua, 2))
  {
    return luaL_error(lua, "Invalid argument(s). Usage:\n"
                      "mesh = mesh_factory.cvt(surface, options) OR\n"
                      "mesh = mesh_factory.cvt(surfaces, options)");
  }

  // Get the options for the Voronoi tessellation.
  lua_getfield(lua, 2, "merge_distance");
  if (!lua_isnil(lua, -1) && !lua_isnumber(lua, -1))
    return luaL_error(lua, "Option 'merge_distance' must be a number.");

  real_t merge_distance = 1e-12;
  if (!lua_isnil(lua, -1))
    merge_distance = (real_t)lua_tonumber(lua, -1);
  lua_pop(lua, 1);
  if (merge_distance <= 0.0)
    return luaL_error(lua, "Option 'merge_distance' must be a number.");

  lua_getfield(lua, 2, "num_interior_points");
  if (!lua_isnumber(lua, -1))
    return luaL_error(lua, "Option 'num_interior_points' must be a number.");

  int num_interior_points = (int)lua_tonumber(lua, -1);
  lua_pop(lua, 1);

  lua_getfield(lua, 2, "iteration_method");
  if (!lua_isnil(lua, -1) && !lua_isstring(lua, -1))
    return luaL_error(lua, "Option 'iteration_method' must be a string.");

  char iteration_method[1024];
  strcpy(iteration_method, "lloyd");
  if (!lua_isnil(lua, -1))
    strcpy(iteration_method, lua_tostring(lua, -1));
  lua_pop(lua, 1);
  if (strcasecmp(iteration_method, "lloyd"))
    return luaL_error(lua, "Invalid iteration method: '%s'\nSupported methods are 'lloyd'.", iteration_method);

  lua_getfield(lua, 2, "num_iterations");
  if (!lua_isnumber(lua, -1))
    return luaL_error(lua, "Option 'num_iterations' must be a number.");

  int num_iterations = (int)lua_tonumber(lua, -1);
  if (num_iterations < 1)
    return luaL_error(lua, "Option 'num_iterations' must be positive.");

  lua_pop(lua, 1);

  // Determine whether the first argument is a surface or a table containing
  // several surfaces. We do this by traversing the table and looking for 
  // certain fields. In the case of a list of surfaces, we make a list of 
  // surface names.
  bool is_surface = false, found_points = false, found_normals = false;
  bool is_surface_list = false, is_invalid = false;

  string_array_t* surface_names = string_array_new();
  int_array_t* num_surface_points = int_array_new();
  ptr_array_t* surface_points = ptr_array_new();
  ptr_array_t* surface_normals = ptr_array_new();
  lua_pushnil(lua);
  while (lua_next(lua, 1))
  {
    static const int key_index = -2;
    static const int val_index = -1;
    const char* key = lua_tostring(lua, key_index);
    if (!is_surface_list && !strcmp(key, "points"))
    {
      if (!lua_ispointlist(lua, val_index))
      {
        is_invalid = true;
        break;
      }
      int num_points;
      ptr_array_append_with_dtor(surface_points, lua_topointlist(lua, val_index, &num_points), DTOR(free));
      if (num_surface_points->size == surface_points->size)
      {
        if (num_surface_points->data[num_surface_points->size-1] != num_points)
        {
          is_invalid = true;
          break;
        }
      }
      else
        int_array_append(num_surface_points, num_points);
      found_points = true;
    }
    else if (!is_surface_list && !strcmp(key, "normals"))
    {
      if (!lua_isvectorlist(lua, val_index))
      {
        is_invalid = true;
        break;
      }
      int num_vectors;
      ptr_array_append_with_dtor(surface_normals, lua_tovectorlist(lua, val_index, &num_vectors), DTOR(free));
      if (num_surface_points->size == surface_normals->size)
      {
        if (num_surface_points->data[num_surface_points->size-1] != num_vectors)
        {
          is_invalid = true;
          break;
        }
      }
      else
        int_array_append(num_surface_points, num_vectors);
      found_normals = true;
    }
    else
    {
      // Something else is in there. It must be a table.
      if (!lua_istable(lua, val_index))
      {
        is_invalid = true;
        break;
      }
      bool found_points = false, found_normals = false;
      lua_pushnil(lua);
      int table_index = val_index - 1;
      while (lua_next(lua, table_index))
      {
        int key_index = -2;
        int val_index = -1;
        const char* key = lua_tostring(lua, key_index);
        if (!is_surface_list && !strcmp(key, "points"))
        {
          if (!lua_ispointlist(lua, val_index))
          {
            is_invalid = true;
            break;
          }
          int num_points;
          ptr_array_append_with_dtor(surface_points, lua_topointlist(lua, val_index, &num_points), DTOR(free));
          if (num_surface_points->size == surface_points->size)
          {
            if (num_surface_points->data[num_surface_points->size-1] != num_points)
            {
              is_invalid = true;
              break;
            }
          }
          else
            int_array_append(num_surface_points, num_points);
          found_points = true;
        }
        else if (!is_surface_list && !strcmp(key, "normals"))
        {
          if (!lua_isvectorlist(lua, val_index))
          {
            is_invalid = true;
            break;
          }
          int num_points;
          ptr_array_append_with_dtor(surface_normals, lua_tovectorlist(lua, val_index, &num_points), DTOR(free));
          if (num_surface_points->size == surface_normals->size)
          {
            if (num_surface_points->data[num_surface_points->size-1] != num_points)
            {
              is_invalid = true;
              break;
            }
          }
          else
            int_array_append(num_surface_points, num_points);
          found_normals = true;
        }
        lua_pop(lua, 1);
      }
      if (is_invalid) break;
      if (found_points && found_normals)
      {
        string_array_append_with_dtor(surface_names, string_dup(key), free_string);
        is_surface_list = true;
      }
    }
    if (!is_surface_list && (found_points && found_normals)) 
    {
      string_array_append_with_dtor(surface_names, string_dup("all"), free_string);
      is_surface = true;
    }
    lua_pop(lua, 1);
  }

  if (is_invalid || (!is_surface && !is_surface_list)) 
  {
    string_array_free(surface_names);
    int_array_free(num_surface_points);
    ptr_array_free(surface_points);
    ptr_array_free(surface_normals);
    return luaL_error(lua, "Argument 1 must be a surface (a table with points and normals fields) or a list of surfaces.");
  }

  // Now we organize the surface or surfaces into a list of points with 
  // normals and tags. We begin by computing the total number of points 
  // on the surfaces, and allocating some containers.
  int num_surfaces = surface_names->size, num_surf_points = 0;
  for (int i = 0; i < num_surfaces; ++i)
    num_surf_points += num_surface_points->data[i];
  point_t* all_surf_points = malloc(sizeof(point_t) * num_surf_points);
  vector_t* all_normals = malloc(sizeof(vector_t) * num_surf_points);
  int* all_surf_indices = malloc(sizeof(int) * num_surf_points);

  // Now toss all the points into a kd-tree.
  int offset = 0;
  for (int i = 0; i < num_surfaces; ++i)
  {
    point_t* spoints = surface_points->data[i];
    memcpy(&all_surf_points[offset], spoints, sizeof(point_t) * num_surface_points->data[i]);
    vector_t* snormals = surface_normals->data[i];
    memcpy(&all_normals[offset], snormals, sizeof(vector_t) * num_surface_points->data[i]);
    for (int j = 0; j < num_surface_points->data[i]; ++j)
      all_surf_indices[offset+j] = i;
    offset += num_surface_points->data[i];
  }
  kd_tree_t* tree = kd_tree_new(all_surf_points, num_surf_points);

  // Now merge the points and collect their normals and tags. Also compute 
  // a bounding box that contains all of the surface points.
  ptr_array_t* merged_points = ptr_array_new();
  ptr_array_t* merged_normals = ptr_array_new();
  ptr_array_t* merged_tags = ptr_array_new();
  bbox_t bbox = {.x1 = FLT_MAX, .x2 = -FLT_MAX, .y1 = FLT_MAX, .y2 = -FLT_MAX, .z1 = FLT_MAX, .z2 = -FLT_MAX};
  for (int i = 0; i < num_surf_points; ++i)
  {
    point_t* point = &all_surf_points[i];
    vector_t* normal = &all_normals[i];
    int surf_index = all_surf_indices[i];
    char* tag = surface_names->data[surf_index];

    // Alter the bounding box if needed.
    bbox_grow(&bbox, point);

    // Search for the points within the merging distance.
    int_slist_t* coincident_points = kd_tree_within_radius(tree, point, merge_distance);
    ASSERT(coincident_points != NULL);

    // If i is the minimum index of all these neighbors, append this point
    // to the list of merged points.
    int min_index = INT_MAX;
    int_slist_node_t* s_iter = coincident_points->front;
    while (s_iter != NULL)
    {
      min_index = MIN(min_index, s_iter->value);
      s_iter = s_iter->next;
    }
    if (i == min_index)
    {
      ptr_array_append(merged_points, point);

      ptr_slist_t* normals = ptr_slist_new();
      ptr_slist_append(normals, normal);
      ptr_array_append_with_dtor(merged_normals, normals, DTOR(ptr_slist_free));

      string_slist_t* tags = string_slist_new();
      string_slist_append(tags, surface_names->data[all_surf_indices[i]]);
      ptr_array_append_with_dtor(merged_tags, tags, DTOR(string_slist_free));
    }
    else
    {
      // Associate the normal vector and the tag with the merged point.
      ptr_slist_t* normals = merged_normals->data[min_index];
      ptr_slist_append(normals, normal);

      string_slist_t* tags = merged_tags->data[min_index];
      string_slist_append(tags, tag);
    }
  }

  // Clean up a bit.
  kd_tree_free(tree);
  ptr_array_free(surface_normals);
  ptr_array_free(surface_points);
  int_array_free(num_surface_points);

  // Initialize a set of stationary points that represent the above surfaces.
  int num_boundary_points, num_tags;
  char** tag_names;
  int_array_t** tags;
  point_t* boundary_points;
  create_boundary_generators(merged_points, merged_normals, merged_tags,
                             &boundary_points, &num_boundary_points,
                             &tag_names, &tags, &num_tags);

  // Initialize a set of interior points that is inside the bounding box 
  // that we computed.
  point_t* interior_points = malloc(sizeof(point_t) * num_interior_points);
  for (int i = 0; i < num_interior_points; ++i)
    point_randomize(&interior_points[i], rand, &bbox);

  // Construct a centroidal voronoi tessellation using Lloyd iteration.
//  mesh_t* mesh = create_cvt_with_lloyd_iteration(boundary_points, num_boundary_points,
//                                                 interior_points, num_interior_points,
//                                                 tag_names, tags, num_tags, num_iterations);

  // Clean up the rest.
  free(interior_points);
  free(boundary_points);
  free(tags);
  free(tag_names);
  string_array_free(surface_names);
  ptr_array_free(merged_tags);
  ptr_array_free(merged_normals);
  ptr_array_free(merged_points);
  free(all_surf_indices);
  free(all_normals);
  free(all_surf_points);

  // Push the mesh onto the stack.
//  lua_pushmesh(lua, mesh);
  return 1;
}
#endif

int mesh_factory_pebi(lua_State* lua)
{
  // Check the arguments.
  int num_args = lua_gettop(lua);
  if ((num_args != 3) && (num_args != 4))
  {
    return luaL_error(lua, "Invalid arguments. Usage:\n"
                      "mesh = mesh_factory.pebi(cell_centers, cell_volumes, faces) OR\n"
                      "mesh = mesh_factory.pebi(cell_centers, cell_volumes, faces, face_centers)");
  }
  if (!lua_ispointlist(lua, 1))
    return luaL_error(lua, "cell_centers must be a list of points.");
  if (!lua_issequence(lua, 2))
    return luaL_error(lua, "cell_volumes must be a sequence of cell volumes.");
  if (!lua_istable(lua, 3))
    return luaL_error(lua, "faces must be a table of 3-tuples containing (cell1, cell2, area).");
  if ((num_args == 4) && !lua_ispointlist(lua, 4))
    return luaL_error(lua, "face_centers must be a list of points.");

  // Get the arguments.

  // Cell center point list.
  int num_cells;
  point_t* cell_centers = lua_topointlist(lua, 1, &num_cells);

  // Cell volume list.
  int num_cell_volumes;
  real_t* cell_volumes = lua_tosequence(lua, 1, &num_cell_volumes);
  if (num_cell_volumes != num_cells)
    return luaL_error(lua, "Number of cell volumes (%d) does not match number of cells (%d).", num_cell_volumes, num_cells);

  // Mine the faces table for all its 3-tuples.
  int num_faces = luaL_len(lua, 3);
  real_t** faces_table_entries = malloc(sizeof(real_t*)*num_faces);
  lua_pushnil(lua);
  int face = 0;
  while (lua_next(lua, 2))
  {
    // Key is at index -2, value is at -1.
    static const int key_index = -2;
    static const int val_index = -1;
    bool key_is_number = lua_isnumber(lua, key_index);
    bool val_is_sequence = lua_issequence(lua, val_index);
    if (!key_is_number || !val_is_sequence)
    {
      lua_pop(lua, 1);
      return luaL_error(lua, "Found non-numeric entries in faces table.");
    }
    int tuple_len;
    faces_table_entries[face] = lua_tosequence(lua, val_index, &tuple_len);
    if (tuple_len != 3)
    {
      lua_pop(lua, 1);
      return luaL_error(lua, "Tuple at index %d of faces table has %d values (should be 3).", key_index, tuple_len);
    }
    ++face;
    lua_pop(lua, 1);
  }
  ASSERT(face == num_faces);

  // Check the faces data.
  for (int f = 0; f < num_faces; ++f)
  {
    real_t* face_tuple = faces_table_entries[f];
    if (face_tuple[0] < 0.0)
      return luaL_error(lua, "Invalid first cell for face %d: %d (must be non-negative).", f, (int)face_tuple[0]);
    if ((face_tuple[1] < 0.0) && (face_tuple[1] != -1.0))
      return luaL_error(lua, "Invalid second cell for face %d: %d (must be non-negative or -1).", f, (int)face_tuple[1]);
    if (face_tuple[2] < 0.0)
      return luaL_error(lua, "Invalid area for face %d: %g.", f, face_tuple[2]);
  }

  // Face centers?
  int num_face_centers = 0;
  point_t* face_centers = NULL;
  if (num_args == 4)
  {
    face_centers = lua_topointlist(lua, 1, &num_face_centers);
    if (num_face_centers != num_faces)
      return luaL_error(lua, "Number of face centers (%d) does not match number of faces (%d).", num_face_centers, num_faces);
  }

  // Shuffle the faces data into canonical form.
  int* faces = malloc(2 * sizeof(int) * num_faces);
  real_t* face_areas = malloc(sizeof(real_t) * num_faces);
  for (int f = 0; f < num_faces; ++f)
  {
    faces[2*f] = faces_table_entries[f][0];
    faces[2*f+1] = faces_table_entries[f][1];
    face_areas[f] = faces_table_entries[f][2];
    free(faces_table_entries[f]);
  }
  free(faces_table_entries);

  // Create the mesh.
  mesh_t* mesh = create_pebi_mesh(MPI_COMM_WORLD, cell_centers, cell_volumes, num_cells, 
                                  faces, face_areas, face_centers, num_faces);
  free(face_areas);
  free(faces);
  free(cell_centers);

  // Pop all the previous arguments off the stack.
  lua_pop(lua, lua_gettop(lua));

  // Push the mesh onto the stack.
  lua_pushmesh(lua, mesh);
  return 1;
}

