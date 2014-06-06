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
#include "geometry/create_tetgen_mesh.h"
#include "geometry/create_dual_mesh.h"
#include "geometry/create_boundary_generators.h"
#include "geometry/rect_prism.h"

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

int mesh_factory_tetgen(lua_State* lua)
{
  // Check the arguments.
  int num_args = lua_gettop(lua);
  if (((num_args != 1) && !lua_isstring(lua, 1)) || 
      ((num_args != 2) && !lua_isstring(lua, 1) && !lua_isnumber(lua, 2)))
  {
    return luaL_error(lua, "Invalid argument(s). Usage:\n"
                      "mesh = mesh_factory.tetgen(mesh_prefix) OR\n"
                      "mesh = mesh_factory.tetgen(mesh_prefix, num_reader_groups).");
  }

  // Make sure the number of reader groups makes sense.
  int num_reader_groups = (num_args == 1) ? -1 : (int)lua_tonumber(lua, 2);
  int nproc;
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  if (num_reader_groups > nproc)
    return luaL_error(lua, "Invalid number of reader groups given: %d", num_reader_groups);

  // Use the mesh prefix to generate filenames.
  const char* mesh_prefix = lua_tostring(lua, 1);
  char node_file[512], ele_file[512], face_file[512], neigh_file[512];
  snprintf(node_file, 512, "%s.node", mesh_prefix);
  snprintf(ele_file, 512, "%s.ele", mesh_prefix);
  snprintf(face_file, 512, "%s.face", mesh_prefix);
  snprintf(neigh_file, 512, "%s.neigh", mesh_prefix);
  mesh_t* mesh = create_tetgen_mesh(MPI_COMM_WORLD, node_file, ele_file, 
                                    face_file, neigh_file, num_reader_groups);

  // Push the mesh onto the stack.
  lua_pushmesh(lua, mesh);
  return 1;
}

int mesh_factory_dual(lua_State* lua)
{
  // Check the arguments.
  int num_args = lua_gettop(lua);
  if ((num_args != 4 ) || (num_args != 5) || 
      !lua_ismesh(lua, 1) || !lua_isstringlist(lua, 2) || 
      !lua_isstringlist(lua, 3) || !lua_isstringlist(lua, 4) || 
      ((num_args == 5) && !lua_isstringlist(lua, 5)))
  {
    return luaL_error(lua, "Invalid argument(s). Usage:\n"
                      "mesh = mesh_factory.dual(original_mesh, external_model_face_tags, model_edge_tags, model_vertex_tags) OR"
                      "mesh = mesh_factory.dual(original_mesh, external_model_face_tags, internal_model_face_tags, model_edge_tags, model_vertex_tags) OR");
  }

  mesh_t* orig_mesh = lua_tomesh(lua, 1);
  int num_external_model_face_tags, num_internal_model_face_tags = 0,
      num_model_edge_tags, num_model_vertex_tags;
  char** external_model_face_tags = lua_tostringlist(lua, 2, &num_external_model_face_tags);
  char** internal_model_face_tags = NULL;
  char** model_edge_tags; 
  char** model_vertex_tags;
  if (num_args == 4)
  {
    model_edge_tags = lua_tostringlist(lua, 3, &num_model_edge_tags);
    model_vertex_tags = lua_tostringlist(lua, 4, &num_model_vertex_tags);
  }
  else
  {
    internal_model_face_tags = lua_tostringlist(lua, 3, &num_external_model_face_tags);
    model_edge_tags = lua_tostringlist(lua, 4, &num_model_edge_tags);
    model_vertex_tags = lua_tostringlist(lua, 5, &num_model_vertex_tags);
  }

  // Make sure the mesh contains the given tags.
  for (int i = 0; i < num_external_model_face_tags; ++i)
  {
    if (!mesh_has_tag(orig_mesh->face_tags, external_model_face_tags[i]))
      return luaL_error(lua, "mesh_factory.dual: Original mesh does not contain face tag '%s'.", external_model_face_tags[i]);
  }
  for (int i = 0; i < num_internal_model_face_tags; ++i)
  {
    if (!mesh_has_tag(orig_mesh->face_tags, internal_model_face_tags[i]))
      return luaL_error(lua, "mesh_factory.dual: Original mesh does not contain face tag '%s'.", internal_model_face_tags[i]);
  }
  for (int i = 0; i < num_model_edge_tags; ++i)
  {
    if (!mesh_has_tag(orig_mesh->edge_tags, model_edge_tags[i]))
      return luaL_error(lua, "mesh_factory.dual: Original mesh does not contain edge tag '%s'.", model_edge_tags[i]);
  }
  for (int i = 0; i < num_model_vertex_tags; ++i)
  {
    if (!mesh_has_tag(orig_mesh->node_tags, model_vertex_tags[i]))
      return luaL_error(lua, "mesh_factory.dual: Original mesh does not contain node tag '%s'.", model_vertex_tags[i]);
  }

  // For now, we only support duals of tet meshes.
  if (!mesh_has_feature(orig_mesh, TETRAHEDRAL))
    return luaL_error(lua, "mesh_factory.dual: A dual mesh can only be created from a tetrahedral mesh.");

  mesh_t* mesh = create_dual_mesh(MPI_COMM_WORLD, orig_mesh, 
                                  external_model_face_tags, num_external_model_face_tags,
                                  internal_model_face_tags, num_internal_model_face_tags,
                                  model_edge_tags, num_model_edge_tags,
                                  model_vertex_tags, num_model_vertex_tags);

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
  real_t* cell_volumes = lua_tosequence(lua, 2, &num_cell_volumes);
  if (num_cell_volumes != num_cells)
    return luaL_error(lua, "Number of cell volumes (%d) does not match number of cells (%d).", num_cell_volumes, num_cells);

  // Mine the faces table for all its 3-tuples.
  int num_faces = luaL_len(lua, 3);
  ASSERT(num_faces >= 4);
  real_t** faces_table_entries = polymec_malloc(sizeof(real_t*)*num_faces);
  lua_pushnil(lua);
  int face = 0;
  while (lua_next(lua, 3))
  {
    // Key is at index -2, value is at -1.
    static const int key_index = -2;
    static const int val_index = -1;
    bool key_is_number = lua_isnumber(lua, key_index);
    bool val_is_sequence = lua_issequence(lua, val_index);
    if (!key_is_number || !val_is_sequence)
    {
      lua_pop(lua, 1);
      polymec_free(faces_table_entries);
      return luaL_error(lua, "Found non-numeric entries in faces table.");
    }
    int tuple_len;
    faces_table_entries[face] = lua_tosequence(lua, val_index, &tuple_len);
    if (tuple_len != 3)
    {
      lua_pop(lua, 1);
      polymec_free(faces_table_entries);
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
  int* faces = polymec_malloc(2 * sizeof(int) * num_faces);
  real_t* face_areas = polymec_malloc(sizeof(real_t) * num_faces);
  for (int f = 0; f < num_faces; ++f)
  {
    faces[2*f] = faces_table_entries[f][0];
    faces[2*f+1] = faces_table_entries[f][1];
    face_areas[f] = faces_table_entries[f][2];
    polymec_free(faces_table_entries[f]);
  }
  polymec_free(faces_table_entries);

  // Create the mesh.
  mesh_t* mesh = create_pebi_mesh(MPI_COMM_WORLD, cell_centers, cell_volumes, num_cells, 
                                  faces, face_areas, face_centers, num_faces);
  polymec_free(face_areas);
  polymec_free(faces);
  polymec_free(cell_centers);

  // Pop all the previous arguments off the stack.
  lua_pop(lua, lua_gettop(lua));

  // Push the mesh onto the stack.
  lua_pushmesh(lua, mesh);
  return 1;
}

