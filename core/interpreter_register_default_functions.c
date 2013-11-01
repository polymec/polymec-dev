// Copyright 2012-2013 Jeffrey Johnson.
// 
// This file is part of Polymec, and is licensed under the Apache License, 
// Version 2.0 (the "License"); you may not use this file except in 
// compliance with the License. You may may find the text of the license in 
// the LICENSE file at the top-level source directory, or obtain a copy of 
// it at
// 
// http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "core/interpreter.h"
#include "core/constant_st_func.h"
#include "core/periodic_bc.h"
#include "core/write_silo.h"

// Lua stuff.
#include "lua.h"
#include "lualib.h"
#include "lauxlib.h"

// This file contains definitions for functions that are bundled with any 
// interpreter instance.

// Creates a point from a set of coordinates.
static int point(lua_State* lua)
{
  // Check the arguments.
  int num_args = lua_gettop(lua);
  if ((num_args != 3) || 
      !lua_isnumber(lua, 1) || !lua_isnumber(lua, 2) || !lua_isnumber(lua, 3))
  {
    return luaL_error(lua, "Arguments must be x, y, z coordinates.");
  }

  double x = lua_tonumber(lua, 1);
  double y = lua_tonumber(lua, 2);
  double z = lua_tonumber(lua, 3);
  point_t* point = point_new(x, y, z);
  lua_pushpoint(lua, point);
  return 1;
}

// Creates a vector from a set of components.
static int vector(lua_State* lua)
{
  // Check the arguments.
  int num_args = lua_gettop(lua);
  if ((num_args != 3) || 
      !lua_isnumber(lua, 1) || !lua_isnumber(lua, 2) || !lua_isnumber(lua, 3))
  {
    return luaL_error(lua, "Arguments must be x, y, z components.");
  }

  double x = lua_tonumber(lua, 1);
  double y = lua_tonumber(lua, 2);
  double z = lua_tonumber(lua, 3);
  vector_t* vec = vector_new(x, y, z);
  lua_pushvector(lua, vec);
  return 1;
}

// Creates a bounding box from a table.
static int bounding_box(lua_State* lua)
{
  // Check the arguments.
  int num_args = lua_gettop(lua);
  if (num_args != 1)
  {
    if (!lua_istable(lua, 1))
      return luaL_error(lua, "Argument must be a table containing x1, x2, y1, y2, z1, z2 values.");
  }

  // Look for x1, x2, y1, y2, z1, z2 in the table.
  bbox_t* bbox = bbox_new(0.0, 1.0, 0.0, 1.0, 0.0, 1.0);
  const char* entries[] = {"x1", "x2", "y1", "y2", "z1", "z2"};
  for (int i = 0; i < 6; ++i)
  {
    lua_pushstring(lua, entries[i]);
    lua_gettable(lua, 1); // Reads name from top, replaces with bounds[name].
    if (!lua_isnumber(lua, -1))
    {
      return luaL_error(lua, "Invalid entry for '%s'.\n"
                        "x1, x2, y1, y2, z1, z2, must all be numbers.", entries[i]);
    }
    switch(i)
    {
      case 0: bbox->x1 = lua_tonumber(lua, -1);
              break;
      case 1: bbox->x2 = lua_tonumber(lua, -1);
              break;
      case 2: bbox->y1 = lua_tonumber(lua, -1);
              break;
      case 3: bbox->y2 = lua_tonumber(lua, -1);
              break;
      case 4: bbox->z1 = lua_tonumber(lua, -1);
              break;
      case 5: bbox->z2 = lua_tonumber(lua, -1);
              break;
      default: break;
    }
    lua_pop(lua, 1); 
  }

  // Push the bounding box onto the stack.
  lua_pushboundingbox(lua, bbox);
  return 1;
}

// Creates a constant function from a number or a 3-tuple.
static int constant_function(lua_State* lua)
{
  // Check the argument.
  int num_args = lua_gettop(lua);
  if (num_args == 1) // Scalar-valued constant.
  {
    if (!lua_isnumber(lua, 1))
      return luaL_error(lua, "Argument must be a number.");
  }
  else if (num_args == 3) // Vector-valued constant.
  {
    for (int i = 0; i < 3; ++i)
    {
      if (!lua_isnumber(lua, i+1))
        return luaL_error(lua, "Argument %d must be a number.", i);
    }
  }
  else
    return luaL_error(lua, "Argument must be a 1 or 3 numbers.");

  // Get the arguments.
  double args[3];
  for (int i = 0; i < num_args; ++i)
    args[i] = lua_tonumber(lua, i+1);

  // Push a constant function onto the stack.
  st_func_t* func = constant_st_func_new(num_args, args);
  if (num_args == 1)
    lua_pushscalarfunction(lua, func);
  else
    lua_pushvectorfunction(lua, func);
  return 1;
}

// Creates a vector-valued function from 3 scalar functions.
static int vector_function(lua_State* lua)
{
  // Check the argument.
  int num_args = lua_gettop(lua);
  if (num_args != 3)
    return luaL_error(lua, "Arguments must be 3 scalar functions.");

  for (int i = 1; i <= 3; ++i)
  {
    if (!lua_isscalarfunction(lua, i))
      return luaL_error(lua, "Argument %d must be a scalar function.", i);
  }

  st_func_t* functions[3];
  for (int i = 1; i <= 3; ++i)
    functions[i-1] = lua_toscalarfunction(lua, i);
  lua_pushvectorfunction(lua, multicomp_st_func_from_funcs("vector function", functions, 3));
  return 1;
}

// Creates a periodic boundary condition from a pair of tags.
static int periodic_bc(lua_State* lua)
{
  // Check the argument.
  int num_args = lua_gettop(lua);
  if (num_args != 2)
    return luaL_error(lua, "Arguments must be 2 boundary mesh (face) tags.");

  for (int i = 1; i <= 3; ++i)
  {
    if (!lua_isstring(lua, i))
      return luaL_error(lua, "Argument %d must be a face tag.", i);
  }

  const char* tag1 = lua_tostring(lua, 1);
  const char* tag2 = lua_tostring(lua, 2);
  periodic_bc_t* bc = periodic_bc_new(tag1, tag2);
  lua_pushuserdefined(lua, bc, NULL);
  return 1;
}

// Evaluates the gradient of a scalar function.
static int grad(lua_State* lua)
{
  // Check the argument.
  int num_args = lua_gettop(lua);
  if ((num_args != 3) && (num_args != 2))
    return luaL_error(lua, "grad: Invalid number of arguments. Must be 2 or 3.");
  
  if ((num_args == 2) && (!lua_isscalarfunction(lua, 1) || 
      (!lua_ispoint(lua, 2) && !lua_ispointlist(lua, 2))))
  {
    return luaL_error(lua, "grad: Arguments must be (F, x).");
  }
  else if ((num_args == 3) && (!lua_isscalarfunction(lua, 1) || 
      (!lua_ispoint(lua, 2) && !lua_ispointlist(lua, 2)) ||
      !lua_isnumber(lua, 3)))
  {
    return luaL_error(lua, "grad: Arguments must be (F, x, t).");
  }

  st_func_t* f = lua_toscalarfunction(lua, 1);
  if (!st_func_has_deriv(f, 1))
    return luaL_error(lua, "grad: Scalar function '%s' has no derivative.", st_func_name(f));

  double t = 0.0;
  if (num_args == 3)
    t = lua_tonumber(lua, 3);

  if (lua_ispoint(lua, 2))
  {
    point_t* x = lua_topoint(lua, 2);
    double val[3];
    st_func_eval_deriv(f, 1, x, t, val);
    vector_t* vec = vector_new(val[0], val[1], val[2]);
    lua_pushvector(lua, vec);
  }
  else
  {
    int num_points;
    point_t* x = lua_topointlist(lua, 2, &num_points);
    vector_t* vecs = malloc(sizeof(vector_t) * num_points);
    for (int i = 0; i < num_points; ++i)
    {
      double val[3];
      st_func_eval_deriv(f, 1, &x[i], t, val);
      vecs[i].x = val[0], vecs[i].y = val[1]; vecs[i].z = val[2];
    }
    lua_pushvectorlist(lua, vecs, num_points);
  }

  return 1;
}

// write_silo_mesh(args) -- This function writes a given mesh to a file 
// on disk. Arguments (passed in a table according to Chapter 5.3 of the 
// Lua reference manual) are:
//
// mesh -> mesh object 
// filename -> name of the file to write (1 file only)
static int lua_write_silo_mesh(lua_State* lua)
{
  // Check the arguments.
  int num_args = lua_gettop(lua);
  if (((num_args == 2) && (!lua_ismesh(lua, 1) || !lua_isstring(lua, 2))) || 
      ((num_args == 3) && (!lua_ismesh(lua, 1) || !lua_istable(lua, 2) || !lua_isstring(lua, 3))) || 
      ((num_args != 2) && (num_args != 3)))
  {
    return luaL_error(lua, "write_silo_mesh: invalid arguments. Usage:\n"
                      "write_silo_mesh(mesh, filename) OR\n"
                      "write_silo_mesh(mesh, fields, filename)");
  }

  // Get the argument(s).
  mesh_t* mesh = lua_tomesh(lua, 1);
  ASSERT(mesh != NULL);
  int N = mesh->num_cells;
  bool has_fields = (num_args == 3) ? lua_istable(lua, 2) : false;
  char* filename = (num_args == 3) ? string_dup(lua_tostring(lua, 3)) : string_dup(lua_tostring(lua, 2));

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
        polymec_error("write_silo_mesh: argument 2 must be a table mapping field names to values.");
      }
      if (lua_issequence(lua, val_index))
      {
        int num_vals;
        double* vals = lua_tosequence(lua, val_index, &num_vals);
        if (num_vals != N)
        {
          lua_pop(lua, 2);
          polymec_error("write_silo_mesh: a scalar field has %d values (should have %d).", num_vals, N);
        }
        vals = NULL;
      }
      else
      {
        int num_vals;
        vector_t* vals = lua_tovectorlist(lua, val_index, &num_vals);
        if (num_vals != N)
        {
          lua_pop(lua, 2);
          polymec_error("write_silo_mesh: a vector field has %d values (should have %d).", num_vals, N);
        }
        vals = NULL;
      }
      lua_pop(lua, 1);
    }
  }

  // Construct a set of fields.
  string_ptr_unordered_map_t* fields = string_ptr_unordered_map_new();
  double* volume = malloc(sizeof(double) * N);
  for (int c = 0; c < N; ++c)
    volume[c] = mesh->cell_volumes[c];
  string_ptr_unordered_map_insert_with_v_dtor(fields, "volume", volume, DTOR(free));

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
        string_ptr_unordered_map_insert(fields, (char*)field_name, field_data);
      }
      else
      {
        ASSERT(lua_isvectorlist(lua, val_index));
        int num_vals;
        vector_t* vector_data = lua_tovectorlist(lua, val_index, &num_vals);
        double* Fx = malloc(sizeof(double) * num_vals);
        double* Fy = malloc(sizeof(double) * num_vals);
        double* Fz = malloc(sizeof(double) * num_vals);
        for (int i = 0; i < num_vals; ++i)
        {
          Fx[i] = vector_data[i].x;
          Fy[i] = vector_data[i].y;
          Fz[i] = vector_data[i].z;
        }
        char Fx_name[128], Fy_name[128], Fz_name[128];
        snprintf(Fx_name, 128, "%s_x", field_name);
        snprintf(Fy_name, 128, "%s_y", field_name);
        snprintf(Fz_name, 128, "%s_z", field_name);
        string_ptr_unordered_map_insert(fields, Fx_name, Fx);
        string_ptr_unordered_map_insert(fields, Fy_name, Fy);
        string_ptr_unordered_map_insert(fields, Fz_name, Fz);
        free(Fx);
        free(Fy);
        free(Fz);
      }
      lua_pop(lua, 1);
    }
  }

  // Write the thing to a mesh file.
  log_info("Writing SILO mesh file with prefix '%s'...", filename);
  write_silo_mesh(mesh, NULL, NULL, NULL, fields, filename, ".", 0, 0.0, 
                  MPI_COMM_SELF, 1, 0);

  // Clean up.
  string_ptr_unordered_map_free(fields);

  return 1;
}

static int lua_write_silo_points(lua_State* lua)
{
  // Check the arguments.
  int num_args = lua_gettop(lua);
  if (((num_args == 2) && (!lua_ispointlist(lua, 1) || !lua_isstring(lua, 2))) || 
      ((num_args == 3) && (!lua_ispointlist(lua, 1) || !lua_istable(lua, 2) || !lua_isstring(lua, 3))) || 
      ((num_args != 2) && (num_args != 3)))
  {
    return luaL_error(lua, "write_silo_points: invalid arguments. Usage:\n"
                      "write_silo_points(points, filename) OR\n"
                      "write_silo_points(points, fields, filename)");
  }

  // Get the argument(s).
  int N;
  point_t* points = lua_topointlist(lua, 1, &N);
  bool has_fields = (num_args == 3) ? lua_istable(lua, 2) : false;
  char* filename = (num_args == 3) ? string_dup(lua_tostring(lua, 3)) : string_dup(lua_tostring(lua, 2));

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
        polymec_error("write_silo_points: argument 2 must be a table mapping field names to values.");
      }
      if (lua_issequence(lua, val_index))
      {
        int num_vals;
        double* vals = lua_tosequence(lua, val_index, &num_vals);
        if (num_vals != N)
        {
          lua_pop(lua, 2);
          polymec_error("write_silo_points: a scalar field has %d values (should have %d).", num_vals, N);
        }
        vals = NULL;
      }
      else
      {
        int num_vals;
        vector_t* vals = lua_tovectorlist(lua, val_index, &num_vals);
        if (num_vals != N)
        {
          lua_pop(lua, 2);
          polymec_error("write_silo_points: a vector field has %d values (should have %d).", num_vals, N);
        }
        vals = NULL;
      }
      lua_pop(lua, 1);
    }
  }

  // Construct a set of fields.
  string_ptr_unordered_map_t* fields = string_ptr_unordered_map_new();
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
        string_ptr_unordered_map_insert(fields, (char*)field_name, field_data);
      }
      else
      {
        ASSERT(lua_isvectorlist(lua, val_index));
        int num_vals;
        vector_t* vector_data = lua_tovectorlist(lua, val_index, &num_vals);
        double* Fx = malloc(sizeof(double) * num_vals);
        double* Fy = malloc(sizeof(double) * num_vals);
        double* Fz = malloc(sizeof(double) * num_vals);
        for (int i = 0; i < num_vals; ++i)
        {
          Fx[i] = vector_data[i].x;
          Fy[i] = vector_data[i].y;
          Fz[i] = vector_data[i].z;
        }
        char Fx_name[128], Fy_name[128], Fz_name[128];
        snprintf(Fx_name, 128, "%s_x", field_name);
        snprintf(Fy_name, 128, "%s_y", field_name);
        snprintf(Fz_name, 128, "%s_z", field_name);
        string_ptr_unordered_map_insert(fields, Fx_name, Fx);
        string_ptr_unordered_map_insert(fields, Fy_name, Fy);
        string_ptr_unordered_map_insert(fields, Fz_name, Fz);
        free(Fx);
        free(Fy);
        free(Fz);
      }
      lua_pop(lua, 1);
    }
  }

  // Write the thing to a file.
  log_info("Writing SILO points file with prefix '%s'...", filename);
  write_silo_points(points, N, fields, filename, ".", 0, 0.0, 
                    MPI_COMM_SELF, 1, 0);

  // Clean up.
  string_ptr_unordered_map_free(fields);

  return 1;
}

static int cell_centers(lua_State* lua)
{
  // Check the arguments.
  int num_args = lua_gettop(lua);
  if ((num_args != 1) || !lua_ismesh(lua, 1))
  {
    return luaL_error(lua, "Invalid argument(s). Usage:\n"
                      "xc = cell_centers(mesh) ->\n"
                      "Returns a list of cell centers for the given mesh.");
  }
  mesh_t* mesh = lua_tomesh(lua, 1);

  // Copy the cell centers to a pointlist and push this list into the interpreter.
  point_t* cc = malloc(sizeof(point_t) * mesh->num_cells);
  memcpy(cc, mesh->cell_centers, sizeof(point_t) * mesh->num_cells);
  lua_pushpointlist(lua, cc, mesh->num_cells);
  return 1;
}

static int cell_tags(lua_State* lua)
{
  // Check the arguments.
  int num_args = lua_gettop(lua);
  if ((num_args != 1) || !lua_ismesh(lua, 1))
  {
    return luaL_error(lua, "Invalid argument(s). Usage:\n"
                      "tags = cell_tags(mesh) ->\n"
                      "Returns a list of names of cell tags for the given mesh.");
  }
  mesh_t* mesh = lua_tomesh(lua, 1);

  lua_newtable(lua);
  int pos = 0, index = 1, *tag_indices, tag_size;
  char* tag_name;
  while (mesh_next_tag(mesh->cell_tags, &pos, &tag_name, &tag_indices, &tag_size))
  {
    lua_pushinteger(lua, index++);
    lua_pushstring(lua, (const char*)tag_name);
    lua_settable(lua, -3);
  }

  return 1;
}

static int cell_tag(lua_State* lua)
{
  // Check the arguments.
  int num_args = lua_gettop(lua);
  if ((num_args != 2) || !lua_ismesh(lua, 1) || !lua_isstring(lua, 2))
  {
    return luaL_error(lua, "Invalid argument(s). Usage:\n"
                      "tag = cell_tag(mesh, tag_name) ->\n"
                      "Returns a list of cell indices associated with the given tag in the mesh.");
  }
  mesh_t* mesh = lua_tomesh(lua, 1);
  const char* tag_name = lua_tostring(lua, 2);

  if (!mesh_has_tag(mesh->cell_tags, tag_name))
    return luaL_error(lua, "The given mesh has no cell tag named '%s'.", tag_name);

  int size;
  int* tag = mesh_tag(mesh->cell_tags, tag_name, &size);
  double* seq = malloc(sizeof(double) * size);
  for (int i = 0; i < size; ++i)
    seq[i] = (double)tag[i];

  lua_pushsequence(lua, seq, size);
  return 1;
}

static int tag_cells(lua_State* lua)
{
  // Check the arguments.
  int num_args = lua_gettop(lua);
  if ((num_args != 3) || !lua_ismesh(lua, 1) || !lua_isstring(lua, 2) || !lua_issequence(lua, 3))
  {
    return luaL_error(lua, "Invalid argument(s). Usage:\n"
                      "tag_cells(mesh, cell_indices, tag_name) ->\n"
                      "Tags a given list of cell indices with the given name within the mesh.");
  }
  mesh_t* mesh = lua_tomesh(lua, 1);
  const char* tag_name = lua_tostring(lua, 2);
  int num_indices = 0;
  double* indices = lua_tosequence(lua, 3, &num_indices);

  // Check the validity of the indices.
  for (int i = 0; i < num_indices; ++i)
  {
    if ((indices[i] < 0) || (indices[i] >= mesh->num_cells))
      return luaL_error(lua, "tag_cells: invalid cell index at %d: %d", i, indices[i]);
  }

  // Overwrite any existing tag.
  if (mesh_has_tag(mesh->cell_tags, tag_name))
    mesh_delete_tag(mesh->cell_tags, tag_name);

  int* tag = mesh_create_tag(mesh->cell_tags, tag_name, num_indices);
  for (int i = 0; i < num_indices; ++i)
    tag[i] = indices[i];

  return 0;
}

static int face_tags(lua_State* lua)
{
  // Check the arguments.
  int num_args = lua_gettop(lua);
  if ((num_args != 1) || !lua_ismesh(lua, 1))
  {
    return luaL_error(lua, "Invalid argument(s). Usage:\n"
                      "tags = face_tags(mesh) ->\n"
                      "Returns a list of names of face tags for the given mesh.");
  }
  mesh_t* mesh = lua_tomesh(lua, 1);

  lua_newtable(lua);
  int pos = 0, index = 1, *tag_indices, tag_size;
  char* tag_name;
  while (mesh_next_tag(mesh->face_tags, &pos, &tag_name, &tag_indices, &tag_size))
  {
    lua_pushinteger(lua, index++);
    lua_pushstring(lua, (const char*)tag_name);
    lua_settable(lua, -3);
  }

  return 1;
}

static int face_tag(lua_State* lua)
{
  // Check the arguments.
  int num_args = lua_gettop(lua);
  if ((num_args != 2) || !lua_ismesh(lua, 1) || !lua_isstring(lua, 2))
  {
    return luaL_error(lua, "Invalid argument(s). Usage:\n"
                      "tag = face_tag(mesh, tag_name) ->\n"
                      "Returns a list of face indices associated with the given tag in the mesh.");
  }
  mesh_t* mesh = lua_tomesh(lua, 1);
  const char* tag_name = lua_tostring(lua, 2);

  if (!mesh_has_tag(mesh->face_tags, tag_name))
    return luaL_error(lua, "The given mesh has no face tag named '%s'.", tag_name);

  int size;
  int* tag = mesh_tag(mesh->face_tags, tag_name, &size);
  double* seq = malloc(sizeof(double) * size);
  for (int i = 0; i < size; ++i)
    seq[i] = (double)tag[i];

  lua_pushsequence(lua, seq, size);
  return 1;
}

static int tag_faces(lua_State* lua)
{
  // Check the arguments.
  int num_args = lua_gettop(lua);
  if ((num_args != 3) || !lua_ismesh(lua, 1) || !lua_isstring(lua, 2) || !lua_issequence(lua, 3))
  {
    return luaL_error(lua, "Invalid argument(s). Usage:\n"
                      "tag_faces(mesh, face_indices, tag_name) ->\n"
                      "Tags a given list of face indices with the given name within the mesh.");
  }
  mesh_t* mesh = lua_tomesh(lua, 1);
  const char* tag_name = lua_tostring(lua, 2);
  int num_indices = 0;
  double* indices = lua_tosequence(lua, 3, &num_indices);

  // Check the validity of the indices.
  for (int i = 0; i < num_indices; ++i)
  {
    if ((indices[i] < 0) || (indices[i] >= mesh->num_faces))
      return luaL_error(lua, "tag_faces: invalid face index at %d: %d", i, indices[i]);
  }

  // Overwrite any existing tag.
  if (mesh_has_tag(mesh->face_tags, tag_name))
    mesh_delete_tag(mesh->face_tags, tag_name);

  int* tag = mesh_create_tag(mesh->face_tags, tag_name, num_indices);
  for (int i = 0; i < num_indices; ++i)
    tag[i] = indices[i];

  return 0;
}

static int edge_tags(lua_State* lua)
{
  // Check the arguments.
  int num_args = lua_gettop(lua);
  if ((num_args != 1) || !lua_ismesh(lua, 1))
  {
    return luaL_error(lua, "Invalid argument(s). Usage:\n"
                      "tags = edge_tags(mesh) ->\n"
                      "Returns a list of names of edge tags for the given mesh.");
  }
  mesh_t* mesh = lua_tomesh(lua, 1);

  lua_newtable(lua);
  int pos = 0, index = 1, *tag_indices, tag_size;
  char* tag_name;
  while (mesh_next_tag(mesh->edge_tags, &pos, &tag_name, &tag_indices, &tag_size))
  {
    lua_pushinteger(lua, index++);
    lua_pushstring(lua, (const char*)tag_name);
    lua_settable(lua, -3);
  }

  return 1;
}

static int edge_tag(lua_State* lua)
{
  // Check the arguments.
  int num_args = lua_gettop(lua);
  if ((num_args != 2) || !lua_ismesh(lua, 1) || !lua_isstring(lua, 2))
  {
    return luaL_error(lua, "Invalid argument(s). Usage:\n"
                      "tag = edge_tag(mesh, tag_name) ->\n"
                      "Returns a list of edge indices associated with the given tag in the mesh.");
  }
  mesh_t* mesh = lua_tomesh(lua, 1);
  const char* tag_name = lua_tostring(lua, 2);

  if (!mesh_has_tag(mesh->edge_tags, tag_name))
    return luaL_error(lua, "The given mesh has no edge tag named '%s'.", tag_name);

  int size;
  int* tag = mesh_tag(mesh->edge_tags, tag_name, &size);
  double* seq = malloc(sizeof(double) * size);
  for (int i = 0; i < size; ++i)
    seq[i] = (double)tag[i];

  lua_pushsequence(lua, seq, size);
  return 1;
}

static int tag_edges(lua_State* lua)
{
  // Check the arguments.
  int num_args = lua_gettop(lua);
  if ((num_args != 3) || !lua_ismesh(lua, 1) || !lua_isstring(lua, 2) || !lua_issequence(lua, 3))
  {
    return luaL_error(lua, "Invalid argument(s). Usage:\n"
                      "tag_edges(mesh, edge_indices, tag_name) ->\n"
                      "Tags a given list of edge indices with the given name within the mesh.");
  }
  mesh_t* mesh = lua_tomesh(lua, 1);
  const char* tag_name = lua_tostring(lua, 2);
  int num_indices = 0;
  double* indices = lua_tosequence(lua, 3, &num_indices);

  // Check the validity of the indices.
  for (int i = 0; i < num_indices; ++i)
  {
    if ((indices[i] < 0) || (indices[i] >= mesh->num_edges))
      return luaL_error(lua, "tag_edges: invalid edge index at %d: %d", i, indices[i]);
  }

  // Overwrite any existing tag.
  if (mesh_has_tag(mesh->edge_tags, tag_name))
    mesh_delete_tag(mesh->edge_tags, tag_name);

  int* tag = mesh_create_tag(mesh->edge_tags, tag_name, num_indices);
  for (int i = 0; i < num_indices; ++i)
    tag[i] = indices[i];

  return 0;
}

static int node_positions(lua_State* lua)
{
  // Check the arguments.
  int num_args = lua_gettop(lua);
  if ((num_args != 1) || !lua_ismesh(lua, 1))
  {
    return luaL_error(lua, "Invalid argument(s). Usage:\n"
                      "xn = node_positions(mesh) ->\n"
                      "Returns a list of node positions for the given mesh.");
  }
  mesh_t* mesh = lua_tomesh(lua, 1);

  // Copy the node positions to a pointlist and push this list into the interpreter.
  point_t* np = malloc(sizeof(point_t) * mesh->num_nodes);
  memcpy(np, mesh->nodes, sizeof(point_t) * mesh->num_nodes);
  lua_pushpointlist(lua, np, mesh->num_nodes);
  return 1;
}

static int node_tags(lua_State* lua)
{
  // Check the arguments.
  int num_args = lua_gettop(lua);
  if ((num_args != 1) || !lua_ismesh(lua, 1))
  {
    return luaL_error(lua, "Invalid argument(s). Usage:\n"
                      "tags = node_tags(mesh) ->\n"
                      "Returns a list of names of node tags for the given mesh.");
  }
  mesh_t* mesh = lua_tomesh(lua, 1);

  lua_newtable(lua);
  int pos = 0, index = 1, *tag_indices, tag_size;
  char* tag_name;
  while (mesh_next_tag(mesh->node_tags, &pos, &tag_name, &tag_indices, &tag_size))
  {
    lua_pushinteger(lua, index++);
    lua_pushstring(lua, (const char*)tag_name);
    lua_settable(lua, -3);
  }

  return 1;
}

static int node_tag(lua_State* lua)
{
  // Check the arguments.
  int num_args = lua_gettop(lua);
  if ((num_args != 2) || !lua_ismesh(lua, 1) || !lua_isstring(lua, 2))
  {
    return luaL_error(lua, "Invalid argument(s). Usage:\n"
                      "tag = node_tag(mesh, tag_name) ->\n"
                      "Returns a list of node indices associated with the given tag in the mesh.");
  }
  mesh_t* mesh = lua_tomesh(lua, 1);
  const char* tag_name = lua_tostring(lua, 2);

  if (!mesh_has_tag(mesh->node_tags, tag_name))
    return luaL_error(lua, "The given mesh has no node tag named '%s'.", tag_name);

  int size;
  int* tag = mesh_tag(mesh->node_tags, tag_name, &size);
  double* seq = malloc(sizeof(double) * size);
  for (int i = 0; i < size; ++i)
    seq[i] = (double)tag[i];

  lua_pushsequence(lua, seq, size);
  return 1;
}

static int tag_nodes(lua_State* lua)
{
  // Check the arguments.
  int num_args = lua_gettop(lua);
  if ((num_args != 3) || !lua_ismesh(lua, 1) || !lua_isstring(lua, 2) || !lua_issequence(lua, 3))
  {
    return luaL_error(lua, "Invalid argument(s). Usage:\n"
                      "tag_nodes(mesh, node_indices, tag_name) ->\n"
                      "Tags a given list of node indices with the given name within the mesh.");
  }
  mesh_t* mesh = lua_tomesh(lua, 1);
  const char* tag_name = lua_tostring(lua, 2);
  int num_indices = 0;
  double* indices = lua_tosequence(lua, 3, &num_indices);

  // Check the validity of the indices.
  for (int i = 0; i < num_indices; ++i)
  {
    if ((indices[i] < 0) || (indices[i] >= mesh->num_nodes))
      return luaL_error(lua, "tag_nodes: invalid node index at %d: %d", i, indices[i]);
  }

  // Overwrite any existing tag.
  if (mesh_has_tag(mesh->node_tags, tag_name))
    mesh_delete_tag(mesh->node_tags, tag_name);

  int* tag = mesh_create_tag(mesh->node_tags, tag_name, num_indices);
  for (int i = 0; i < num_indices; ++i)
    tag[i] = indices[i];

  return 0;
}

void interpreter_register_default_functions(interpreter_t* interp)
{
  interpreter_register_function(interp, "point", point);
  interpreter_register_function(interp, "vector", vector);
  interpreter_register_function(interp, "bounding_box", bounding_box);
  interpreter_register_function(interp, "constant_function", constant_function);
  interpreter_register_function(interp, "vector_function", vector_function);
  interpreter_register_function(interp, "periodic_bc", periodic_bc);
  interpreter_register_function(interp, "grad", grad);
  interpreter_register_function(interp, "write_silo_mesh", lua_write_silo_mesh);
  interpreter_register_function(interp, "write_silo_points", lua_write_silo_points);

  // Mesh query functions.
  interpreter_register_function(interp, "cell_centers", cell_centers);
  interpreter_register_function(interp, "cell_tags", cell_tags);
  interpreter_register_function(interp, "cell_tag", cell_tag);
  interpreter_register_function(interp, "tag_cells", tag_cells);
  interpreter_register_function(interp, "face_tags", face_tags);
  interpreter_register_function(interp, "face_tag", face_tag);
  interpreter_register_function(interp, "tag_faces", tag_faces);
  interpreter_register_function(interp, "edge_tags", edge_tags);
  interpreter_register_function(interp, "edge_tag", edge_tag);
  interpreter_register_function(interp, "tag_edges", tag_edges);
  interpreter_register_function(interp, "node_positions", node_positions);
  interpreter_register_function(interp, "node_tags", node_tags);
  interpreter_register_function(interp, "node_tag", node_tag);
  interpreter_register_function(interp, "tag_nodes", tag_nodes);
}

