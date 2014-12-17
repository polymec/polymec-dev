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

// point_factory.c - Implementations of interpreter functions for generating
// sets of points.

#include <strings.h>
#include "core/polymec.h"
#include "core/point.h"
#include "core/interpreter.h"
#include "core/slist.h"
#include "core/array.h"
#include "core/kd_tree.h"
#include "geometry/generate_random_points.h"
#include "geometry/create_point_lattice.h"

// Lua stuff.
#include "lua.h"
#include "lualib.h"
#include "lauxlib.h"

int point_factory_cubic_lattice(lua_State* lua)
{
  // The argument should be a single table of named values.
  int num_args = lua_gettop(lua);
  if ((num_args != 1) || (!lua_istable(lua, 1)))
  {
    return luaL_error(lua, "Invalid arguments. Usage:\n"
                      "points = point_factory.cubic_lattice{nx, ny, nz, bounds}");
  }

  // Extract arguments.
  const char* entries[] = {"nx", "ny", "nz", "bounds"};
  int nx = 0, ny = 0, nz = 0;
  bbox_t* bbox = NULL;
  for (int i = 0; i < 4; ++i)
  {
    lua_pushstring(lua, entries[i]);
    lua_gettable(lua, 1);
    if (i < 3)
    {
      if (!lua_isnumber(lua, -1))
        return luaL_error(lua, "Missing integer argument: %s", entries[i]);

      switch(i) 
      {
        case 1: nx = (int)lua_tonumber(lua, -1); break;
        case 2: ny = (int)lua_tonumber(lua, -1); break;
        case 3: nz = (int)lua_tonumber(lua, -1); break;
        default: break;
      }
    }
    else 
    {
      if (!lua_isboundingbox(lua, -1))
        return luaL_error(lua, "bounds must be a bounding box.");

      bbox = lua_toboundingbox(lua, -1);
    }
  }

  // Validate arguments.
  if ((nx <= 0) || (ny <= 0) || (nz <= 0))
    return luaL_error(lua, "nx, ny, and nz must all be positive.");

  if (bbox->x1 >= bbox->x2)
    return luaL_error(lua, "In bounding_box: x1 must be less than x2.");

  if (bbox->y1 >= bbox->y2)
    return luaL_error(lua, "In bounding_box: y1 must be less than y2.");

  if (bbox->z1 >= bbox->z2)
    return luaL_error(lua, "In bounding_box: z1 must be less than z2.");

  // Pop all the previous arguments off the stack.
  lua_pop(lua, lua_gettop(lua));

  // Create the lattice of points.
  int num_points = nx * ny * nz, offset = 0;
  point_t* points = polymec_malloc(sizeof(point_t) * num_points);
  real_t dx = (bbox->x2 - bbox->x1) / nx;
  real_t dy = (bbox->y2 - bbox->y1) / ny;
  real_t dz = (bbox->z2 - bbox->z1) / nz;
  for (int i = 0; i < nx; ++i)
  {
    real_t xi = bbox->x1 + (i+0.5) * dx;
    for (int j = 0; j < ny; ++j)
    {
      real_t yj = bbox->y1 + (j+0.5) * dy;
      for (int k = 0; k < nz; ++k, ++offset)
      {
        real_t zk = bbox->z1 + (k+0.5) * dz;
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

int point_factory_cylinder(lua_State* lua)
{
  // The argument should be a single table of named values.
  int num_args = lua_gettop(lua);
  if ((num_args != 1) || (!lua_istable(lua, 1)))
  {
    return luaL_error(lua, "Invalid arguments. Usage:\n"
                      "points = point_factory.cylinder{radius, length, center = {0, 0, 0}, nr, nz, axis = {0, 0, 1}, radial_spacing = 'linear', log_spacing_factor = 1.1, num_ghost = 1}");
  }

  // Extract arguments.
  const char* entries[] = {"radius", "length", "center", "nr", "nz", "axis", "radial_spacing", "log_spacing_factor", "num_ghost"};
  real_t radius = 0.0, length = 0.0;
  real_t log_spacing_factor = 1.1;
  point_t* x0 = NULL;
  vector_t* axis = NULL;
  int nr = 0, nz = 0, ng = 0;
  const char* radial_spacing = NULL;
  for (int i = 0; i < 8; ++i)
  {
    lua_pushstring(lua, entries[i]);
    lua_gettable(lua, 1);
    if ((i > 2) && (i != 5) && (i != 6))
    {
      if ((i == 8) && lua_isnil(lua, -1)) continue;

      if ((i < 5) && !lua_isnumber(lua, -1))
        return luaL_error(lua, "Missing integer argument: %s", entries[i]);

      switch(i) 
      {
        case 3: nr = (int)lua_tonumber(lua, -1); break;
        case 4: nz = (int)lua_tonumber(lua, -1); break;
        case 8: ng = (int)lua_tonumber(lua, -1); break;
        default: break;
      }
    }
    else if (i == 0)
    {
      if (!lua_isnumber(lua, -1))
        return luaL_error(lua, "radius must be a positive number.");

      radius = lua_tonumber(lua, -1);
    }
    else if (i == 1)
    {
      if (!lua_isnumber(lua, -1))
        return luaL_error(lua, "length must be a positive number.");

      length = lua_tonumber(lua, -1);
    }
    else if (i == 2)
    {
      if (!lua_isnil(lua, -1))
      {
        if (!lua_ispoint(lua, -1))
          return luaL_error(lua, "center must be a point.");
        x0 = lua_topoint(lua, -1);
      }
    }
    else if (i == 5)
    {
      if (!lua_isnil(lua, -1))
      {
        if (!lua_isvector(lua, -1))
          return luaL_error(lua, "axis must be a vector.");
        axis = lua_tovector(lua, -1);
        if (vector_mag(axis) == 0.0)
          return luaL_error(lua, "axis must not be the zero vector.");
        vector_normalize(axis);
      }
    }
    else if (i == 6)
    {
      if (!lua_isnil(lua, -1))
      {
        if (!lua_isstring(lua, -1))
          return luaL_error(lua, "radial_spacing must be 'linear' or 'log'.");
        radial_spacing = lua_tostring(lua, -1);
      }
    }
    else if (i == 7)
    {
      if (!lua_isnil(lua, -1))
      {
        if (!lua_isnumber(lua, -1))
          return luaL_error(lua, "log_spacing_factor must be positive.");
        log_spacing_factor = lua_tonumber(lua, -1);
      }
    }
  }
  if (x0 == NULL)
    x0 = point_new(0.0, 0.0, 0.0);
  if (axis == NULL)
    axis = vector_new(0.0, 0.0, 1.0);

  // Validate inputs.
  if (radius <= 0)
    return luaL_error(lua, "radius must be positive.");

  if (length <= 0)
    return luaL_error(lua, "length must be positive.");

  if ((nr <= 0) || (nz <= 0))
    return luaL_error(lua, "nr and nz must all be positive.");

  if ((radial_spacing != NULL) && 
      !strcasecmp(radial_spacing, "linear") &&
      !strcasecmp(radial_spacing, "log"))
  {
    return luaL_error(lua, "radial_spacing must be 'linear' or 'log'.");
  }

  if (log_spacing_factor <= 0.0)
    return luaL_error(lua, "log_spacing_factor must be positive.");

  if (ng < 0)
    return luaL_error(lua, "num_ghost must be non-negative.");

  // Pop all the previous arguments off the stack.
  lua_pop(lua, lua_gettop(lua));

  // Determine the radial spacing.
  real_t linear_dr = radius / (1.0*nr - 0.5);
  real_t r[nr-1+ng], dr[nr-1+ng];
  if ((radial_spacing == NULL) || (!strcasecmp(radial_spacing, "linear")))
  {
    for (int j = 0; j < nr-1+ng; ++j)
    {
      dr[j] = linear_dr;
      r[j] = (j+1) * linear_dr;
    }
  }
  else
  {
    // Figure out logarithmic spacing.
    real_t sum = 0.0;
    for (int j = 0; j < nr; ++j)
      sum += pow(log_spacing_factor, 1.0*j);
    real_t dr0 = radius / sum;
    for (int j = 0; j < nr-1+ng; ++j)
    {
      if (j < nr-1)
        dr[j] = pow(log_spacing_factor, 1.0*j) * dr0;
      else
        dr[j] = dr[j-1];
      r[j] = (j == 0) ? dr[j] : r[j-1] + dr[j];
    }
  }

  // Create the points. They should be a set of concentric points with equal
  // angular spacing, with a single point at the center.
//  int num_points_in_disk = 1 + (nr - 1) * ntheta;
  int num_points_in_disk = 1;
  for (int i = 0; i < nr-1+ng; ++i)
  {
    real_t dtheta = dr[i]/ r[i];
    int ntheta = (int)(2.0 * M_PI / dtheta);
    num_points_in_disk += ntheta;
  }
  int num_disks = nz;
  int num_points = num_points_in_disk * (nz + 2*ng);
  point_t* points = polymec_malloc(sizeof(point_t) * num_points);
  real_t dz = length / nz;

  // Set up an orthonormal basis for the given axis.
  vector_t e1, e2;
  compute_orthonormal_basis(axis, &e1, &e2);

  // Find the point at the "bottom" of the cylinder.
  point_t x_bottom;
  x_bottom.x = x0->x - 0.5*length*axis->x;
  x_bottom.y = x0->y - 0.5*length*axis->y;
  x_bottom.z = x0->z - 0.5*length*axis->z;
  int offset = 0;
  for (int i = -ng; i < num_disks+ng; ++i)
  {
    // Find the location of the center of the disk.
    point_t x_center = {.x = x_bottom.x + (i+0.5) * dz * axis->x,
                        .y = x_bottom.y + (i+0.5) * dz * axis->y,
                        .z = x_bottom.z + (i+0.5) * dz * axis->z};

    // Plant a point in the center of the disk.
    points[offset++] = x_center;

    // Construct the other points in the disk.
    for (int j = 0; j < nr-1+ng; ++j)
    {
      real_t rj = r[j];
      real_t dtheta = dr[j]/ rj;
      int ntheta = (int)(2.0 * M_PI / dtheta);
      dtheta = 2.0 * M_PI / ntheta; // Re-adjust.
      for (int k = 0; k < ntheta; ++k, ++offset)
      {
        real_t thetak = k * dtheta;
        real_t cos_thetak = cos(thetak);
        real_t sin_thetak = sin(thetak);
        points[offset].x = x_center.x + rj * (cos_thetak*e1.x + sin_thetak*e2.x);
        points[offset].y = x_center.y + rj * (cos_thetak*e1.y + sin_thetak*e2.y);
        points[offset].z = x_center.z + rj * (cos_thetak*e1.z + sin_thetak*e2.z);
      }
    }
  }
  ASSERT(offset == num_points);

  // Push the points onto the stack.
  lua_pushpointlist(lua, points, num_points);
  return 1;
}

static int read_ascii_stl_file(FILE* stl_file, 
                               ptr_array_t* all_vertices, 
                               ptr_array_t* all_normals,
                               char* error_message)
{
  int status = 0;
  vector_t* n;
  point_t *v1, *v2, *v3;
  do
  {
    n = polymec_malloc(sizeof(vector_t));
    v1 = polymec_malloc(sizeof(point_t));
    v2 = polymec_malloc(sizeof(point_t));
    v3 = polymec_malloc(sizeof(point_t));
    status = fscanf(stl_file, "facet normal %le %le %le\n", &n->x, &n->y, &n->z);
    if (status != 3)
    {
      snprintf(error_message, 1024, "Problem reading facet %zd.", all_normals->size/3);
      goto exit_on_error;
    }
    status = fscanf(stl_file, "outer loop");
    if (status != 0)
    {
      snprintf(error_message, 1024, "Problem reading outer loop header for facet %zd.", all_normals->size/3);
      goto exit_on_error;
    }
    fscanf(stl_file, "\n");
    status = fscanf(stl_file, "vertex %le %le %le\n", &v1->x, &v1->y, &v1->z);
    if (status != 3)
    {
      snprintf(error_message, 1024, "Problem reading vertex 1 for facet %zd.", all_normals->size/3);
      goto exit_on_error;
    }
    status = fscanf(stl_file, "vertex %le %le %le\n", &v2->x, &v2->y, &v2->z);
    if (status != 3)
    {
      snprintf(error_message, 1024, "Problem reading vertex 2 for facet %zd.", all_normals->size/3);
      goto exit_on_error;
    }
    status = fscanf(stl_file, "vertex %le %le %le\n", &v3->x, &v3->y, &v3->z);
    if (status != 3)
    {
      snprintf(error_message, 1024, "Problem reading vertex 3 for facet %zd.", all_normals->size/3);
      goto exit_on_error;
    }
    status = fscanf(stl_file, "endloop");
    if (status != 0)
    {
      snprintf(error_message, 1024, "Problem reading outer loop footer for facet %zd.", all_normals->size/3);
      goto exit_on_error;
    }
    fscanf(stl_file, "\n");
    status = fscanf(stl_file, "endfacet\n");
    if (status != 0)
    {
      snprintf(error_message, 1024, "Problem reading footer for facet %zd.", all_normals->size/3);
      goto exit_on_error;
    }

    // Add the entries.
    ptr_array_append(all_normals, n);
    ptr_array_append(all_normals, n);
    // FIXME: The following line results in a branch from an uninitialized value.
    ptr_array_append_with_dtor(all_normals, n, DTOR(free)); 
    ptr_array_append_with_dtor(all_vertices, v1, DTOR(free));
    ptr_array_append_with_dtor(all_vertices, v2, DTOR(free));
    ptr_array_append_with_dtor(all_vertices, v3, DTOR(free));

    // Have we reached the end?
    char solid_name[1024];
    status = fscanf(stl_file, "endsolid%s", solid_name);
    if (status == 1) break;
  }
  while (status != EOF);
  ASSERT(all_normals->size == all_vertices->size);
  return 0;

exit_on_error:
  polymec_free(n);
  polymec_free(v1);
  polymec_free(v2);
  polymec_free(v3);
  return -1;
}

static int read_binary_stl_file(FILE* stl_file, 
                                ptr_array_t* all_vertices, 
                                ptr_array_t* all_normals,
                                char* error_message)
{
  // Skip the 80-byte header.
  fseek(stl_file, 80*sizeof(char), SEEK_SET);
  
  // Read the number of triangles.
  int num_triangles;
  fread(&num_triangles, sizeof(int), 1, stl_file);
  if (num_triangles <= 0)
  {
    snprintf(error_message, 1024, "Non-positive number of triangles read from binary STL file.");
    return -1;
  }

  // Read all of the triangle information.
  for (int i = 0; i < num_triangles; ++i)
  {
    // Each triangle consists of 12 floats: 3 for a normal vector, and 3x3 
    // for each of the coordinates of the vertices.
    float floats[12];
    int num_floats = fread(floats, sizeof(float), 12, stl_file);
    if (num_floats != 12)
    {
      snprintf(error_message, 1024, "Error reading data for facet %d (%d/12 entries read).", i, num_floats);
      return -1;
    }

    // There should be 2 bytes of "attribute byte count" data that is usually 
    // not used.
    short int attr_byte_count;
    int num_shorts = fread(&attr_byte_count, sizeof(short), 1, stl_file);
    if (num_shorts != 1)
    {
      snprintf(error_message, 1024, "Error reading attribute byte count for facet %d.", i);
      return -1;
    }

    vector_t* n = polymec_malloc(sizeof(vector_t));
    n->x = floats[0], n->y = floats[1], n->z = floats[2];
    vector_t* v1 = polymec_malloc(sizeof(vector_t));
    v1->x = floats[3], v1->y = floats[4], v1->z = floats[5];
    vector_t* v2 = polymec_malloc(sizeof(vector_t));
    v2->x = floats[6], v2->y = floats[7], v2->z = floats[8];
    vector_t* v3 = polymec_malloc(sizeof(vector_t));
    v3->x = floats[9], v3->y = floats[10], v3->z = floats[11];

    // Add the entries.
    ptr_array_append(all_normals, n);
    ptr_array_append(all_normals, n);
    ptr_array_append_with_dtor(all_normals, n, DTOR(free));
    ptr_array_append_with_dtor(all_vertices, v1, DTOR(free));
    ptr_array_append_with_dtor(all_vertices, v2, DTOR(free));
    ptr_array_append_with_dtor(all_vertices, v3, DTOR(free));
  }
  ASSERT(all_normals->size == all_vertices->size);

  return 0;
}

// Import a set of points on a surface from an STL file.
static void import_points_from_stl(const char* stl_file_name, int* num_points, point_t** points, vector_t** normals, char* error_message)
{
  // Notice that realloc and garbage collection don't interact nicely 
  // with each other, so we don't use garbage-collected points and vectors
  // in this method.
  ptr_array_t* all_vertices = ptr_array_new();
  ptr_array_t* all_normals = ptr_array_new();

  // Read the header to determine whether it's an ASCII or binary STL file, 
  // or whether it's not an STL file at all.
  FILE* stl_file = fopen(stl_file_name, "r");
  int status = 0;
  bool is_ascii = false;
  {
    char solid_name[1024];
    status = fscanf(stl_file, "solid%s\n", solid_name);
    if (status == 1)
      is_ascii = true;
    else
    {
      // Definitely not an ASCII STL file, so try binary.
      fclose(stl_file);
      stl_file = fopen(stl_file_name, "rb");
    }
  }

  // Read all of the triangular facets.
  if (is_ascii)
    status = read_ascii_stl_file(stl_file, all_vertices, all_normals, error_message);
  else
    status = read_binary_stl_file(stl_file, all_vertices, all_normals, error_message);
  if (status != 0)
    goto exit_on_error;

  // Dump the points into a kd-tree so that we can see which ones coincide.
  kd_tree_t* point_tree;
  {
    point_t* points_with_duplicates = polymec_malloc(sizeof(point_t) * all_vertices->size);
    for (int i = 0; i < all_vertices->size; ++i)
    {
      point_t* v = all_vertices->data[i];
      points_with_duplicates[i] = *v;
    }
    point_tree = kd_tree_new(points_with_duplicates, all_vertices->size);
    polymec_free(points_with_duplicates);
  }

  // Now make a list of unique point indices (the indices of those points 
  // with the lowest index that are coincident with other points in the tree).
  int_array_t* unique_point_indices = int_array_new();
  ptr_array_t* averaged_normals = ptr_array_new();
  for (int i = 0; i < all_vertices->size; ++i)
  {
    point_t* point = all_vertices->data[i];
    int_array_t* coincident_points = kd_tree_within_radius(point_tree, point, 1e-12);

    // NOTE: To find the normal vector at this point, we average 
    // NOTE: the normals of the triangles that include this point.
    // WARNING: This assumes that the surface is sufficiently smooth!

    if (coincident_points == NULL)
    {
      // I think this can only happen in open surfaces.
      int_array_append(unique_point_indices, i);
      vector_t* n_avg = polymec_malloc(sizeof(vector_t));
      ptr_array_append_with_dtor(averaged_normals, n_avg, DTOR(free));
    }
    else
    {
      ASSERT(coincident_points->size > 0);
      int min_index = INT_MAX;
      vector_t* n_avg = polymec_malloc(sizeof(vector_t)); 
      for (int j = 0; j < coincident_points->size; ++j)
      {
        int value = coincident_points->data[j];
        min_index = MIN(value, min_index);
        vector_t* n = all_normals->data[value];
        n_avg->x = n->x;
        n_avg->y = n->y;
        n_avg->z = n->z;
      }
      if (i == min_index)
      {
        int_array_append(unique_point_indices, i);
        n_avg->x /= coincident_points->size;
        n_avg->y /= coincident_points->size;
        n_avg->z /= coincident_points->size;
        ptr_array_append_with_dtor(averaged_normals, n_avg, DTOR(free));
      }
      else
        polymec_free(n_avg);
      int_array_free(coincident_points);
    }
  }
  ASSERT(unique_point_indices->size > 0);
  ASSERT(unique_point_indices->size == averaged_normals->size);

  // Now build our list of points and normals.
  *num_points = unique_point_indices->size;
  *points = polymec_malloc(sizeof(point_t) * unique_point_indices->size);
  *normals = polymec_malloc(sizeof(vector_t) * unique_point_indices->size);
  for (int i = 0; i < unique_point_indices->size; ++i)
  {
    point_t* p = all_vertices->data[unique_point_indices->data[i]];
    (*points)[i] = *p;
    vector_t* n = averaged_normals->data[i];
    (*normals)[i] = *n; 
  }

  // Clean up.
  ptr_array_free(averaged_normals);
  int_array_free(unique_point_indices);
  kd_tree_free(point_tree);
  ptr_array_free(all_vertices);
  ptr_array_free(all_normals);
  return;

exit_on_error:
  ptr_array_free(all_vertices);
  ptr_array_free(all_normals);
  *points = NULL;
  *normals = NULL;
  *num_points = 0;
  return;
}

int point_factory_import_from_cad(lua_State* lua)
{
  int num_args = lua_gettop(lua);
  if ((num_args != 1) || (!lua_isstring(lua, 1)))
  {
    return luaL_error(lua, "Invalid arguments. Usage:\n"
                      "points, normals = point_factory.import_from_cad(cad_file_name)");
  }
  const char* cad_file_name = lua_tostring(lua, 1);

  // Determine the type of file from its suffix.
  char* suffix = strstr(cad_file_name, ".");
  char* next;
  while ((next = strstr(&suffix[1], ".")) != NULL)
    suffix = next;
  if (suffix == NULL)
    return luaL_error(lua, "Argument must be a filename ending in a suffix that indicates its CAD format.");

  // Make sure the suffix indicates a supported format.
  static const char* supported_suffixes[] = {".stl", NULL};
  int i = 0;
  while ((supported_suffixes[i] != NULL) && strcasecmp(suffix, supported_suffixes[i]))
    ++i;
  if (supported_suffixes[i] == NULL)
    return luaL_error(lua, "Unsupported file format: %s", suffix);
  
  // Check for the existence of the CAD file on disk.
  FILE* cad_file = fopen(cad_file_name, "r");
  if (cad_file == NULL)
    return luaL_error(lua, "Could not open file '%s'", cad_file_name);
  fclose(cad_file);

  // Read the data and generate the list of points.
  point_t* points = NULL;
  vector_t* normals = NULL;
  char error_message[1024];
  int num_points;
  if (!strcasecmp(suffix, ".stl"))
    import_points_from_stl(cad_file_name, &num_points, &points, &normals, error_message);

  if (points == NULL)
    return luaL_error(lua, "Error reading %s: %s", cad_file_name, error_message);

  // Close up shop and return the points and normals as fields in a table.
  lua_createtable(lua, 0, 2);
  lua_pushpointlist(lua, points, num_points);
  lua_setfield(lua, -2, "points");
  lua_pushvectorlist(lua, normals, num_points);
  lua_setfield(lua, -2, "normals");
  return 1;
}

int point_factory_random_points(lua_State* lua)
{
  // Check the arguments.
  int num_args = lua_gettop(lua);
  if ((num_args != 2) && (num_args != 3))
    return luaL_error(lua, "Invalid arguments. Usage:\npoints = random_points(N, bounding_box) OR\npoints = random_points(N, density, bounding_box)");

  // Get the arguments.
  int N = (int)lua_tonumber(lua, 1);
  if (N <= 0)
    return luaL_error(lua, "Invalid (nonpositive) number of points.");

  sp_func_t* density = NULL;
  bbox_t* bbox = NULL;
  if (num_args == 2)
  {
    if (!lua_isboundingbox(lua, 2))
      return luaL_error(lua, "Second argument must be a bounding box.");

    bbox = lua_toboundingbox(lua, 2);
    ASSERT(bbox != NULL);
    real_t one = 1.0;
    density = constant_sp_func_new(1, &one);
  }
  else
  {
    if (!lua_isscalarfunction(lua, 2))
      return luaL_error(lua, "Second argument must be a scalar function.");

    st_func_t* density_t = lua_toscalarfunction(lua, 2);
    density = st_func_freeze(density_t, 0.0);
    if (!lua_isboundingbox(lua, 3))
      return luaL_error(lua, "Third argument must be a bounding box.");

    bbox = lua_toboundingbox(lua, 3);
  }

  point_t* points = polymec_malloc(sizeof(point_t) * N);
  rng_t* rng = host_rng_new();
  generate_random_points(rng, density, bbox, N, points);

  // Return the point list.
  lua_pushpointlist(lua, points, N);
  return 1;
}

int point_factory_ccp_points(lua_State* lua)
{
  // Check the arguments.
  int num_args = lua_gettop(lua);
  if (num_args != 4)
  {
    return luaL_error(lua, "Invalid arguments. Usage:\n"
                      "points = ccp_points(Nx, Ny, Nz, bounding_box)");
  }

  // Get the arguments.
  int Nx = (int)lua_tonumber(lua, 1);
  int Ny = (int)lua_tonumber(lua, 2);
  int Nz = (int)lua_tonumber(lua, 3);
  if ((Nx <= 0) || (Ny <= 0) || (Nz <= 0))
    return luaL_error(lua, "Nx, Ny, and Nz must all be positive.");

  if (!lua_isboundingbox(lua, 4))
    return luaL_error(lua, "Fourth argument must be a bounding box.");

  bbox_t* bbox = lua_toboundingbox(lua, 4);

  // Create the point list.
  ptr_slist_t* point_list = ptr_slist_new();
  real_t dx = (bbox->x2 - bbox->x1) / Nx,
         dy = (bbox->y2 - bbox->y1) / Ny,
         dz = (bbox->z2 - bbox->z1) / Nz;
  for (int i = 0; i < Nx; ++i)
  {
    bool x1face = (i > 0);
    bool x2face = (i < Nx-1);
    real_t x1 = bbox->x1 + i * dx;
    real_t x2 = x1 + dx;
    for (int j = 0; j < Ny; ++j)
    {
      bool y1face = (j > 0);
      bool y2face = (j < Ny-1);
      real_t y1 = bbox->x1 + j * dy;
      real_t y2 = y1 + dy;
      for (int k = 0; k < Nz; ++k)
      {
        bool z1face = (k > 0);
        bool z2face = (k < Ny-1);
        real_t z1 = bbox->x1 + k * dz;
        real_t z2 = z1 + dz;

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
  ASSERT(num_points > 0);
  point_t* points = polymec_malloc(sizeof(point_t) * num_points);
  for (int i = 0; i < num_points; ++i)
    point_copy(&points[i], ptr_slist_pop(point_list, NULL));
  ptr_slist_free(point_list);
  lua_pushpointlist(lua, points, num_points);
  return 1;
}

