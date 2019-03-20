// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// This file implements unimesh-related methods for the silo_file class.

#include "core/polymec.h"

#include "silo.h"
#include "core/arch.h"
#include "core/logging.h"
#include "core/array.h"
#include "core/array_utils.h"
#include "core/timer.h"
#include "io/silo_file.h"

#if POLYMEC_HAVE_DOUBLE_PRECISION
#define SILO_FLOAT_TYPE DB_DOUBLE
#else
#define SILO_FLOAT_TYPE DB_FLOAT
#endif

// These functions are implemented in io/silo_file.c and used
// here, even though they are not part of polymec's API.
extern MPI_Comm silo_file_comm(silo_file_t* file);
extern DBfile* silo_file_dbfile(silo_file_t* file);
extern DBoptlist* optlist_from_metadata(field_metadata_t* metadata, int component);
extern void optlist_free(DBoptlist* optlist);
extern void silo_file_add_subdomain_mesh(silo_file_t* file, const char* mesh_name, int silo_mesh_type, DBoptlist* optlist);
extern void silo_file_add_subdomain_field(silo_file_t* file, const char* mesh_name, const char* field_name, int silo_field_type, DBoptlist* optlist);
extern void silo_file_push_domain_dir(silo_file_t* file);
extern void silo_file_pop_dir(silo_file_t* file);
extern string_ptr_unordered_map_t* silo_file_scratch(silo_file_t* file);
extern void silo_file_write_field_metadata(silo_file_t* file, const char* md_name, field_metadata_t* md);
extern void silo_file_read_field_metadata(silo_file_t* file, const char* md_name, field_metadata_t* md);
extern bool silo_file_contains_field_metadata(silo_file_t* file, const char* md_name);

static void write_unimesh_patch_grid(silo_file_t* file,
                                     const char* patch_grid_name,
                                     int nx, int ny, int nz,
                                     bbox_t* bbox,
                                     coord_mapping_t* mapping)
{
  ASSERT(nx > 0);
  ASSERT(ny > 0);
  ASSERT(nz > 0);

  silo_file_push_domain_dir(file);

  DBfile* dbfile = silo_file_dbfile(file);

  // Name the coordinate axes.
  const char* const coord_names[3] = {"x", "y", "z"};

  // Provide coordinates.
  int N = (nx+1)*(ny+1)*(nz+1);
  real_t* x_node = polymec_malloc(sizeof(real_t) * N);
  real_t* y_node = polymec_malloc(sizeof(real_t) * N);
  real_t* z_node = polymec_malloc(sizeof(real_t) * N);
  int dimensions[3] = {nx+1, ny+1, nz+1};
  real_t dx = (bbox->x2 - bbox->x1) / nx,
         dy = (bbox->y2 - bbox->y1) / ny,
         dz = (bbox->z2 - bbox->z1) / nz;
  int coord_type;
  if (mapping != NULL)
  {
    coord_type = DB_NONCOLLINEAR;
    point_t x;
    int l = 0;
    for (int i = 0; i <= nx; ++i)
    {
      x.x = bbox->x1 + i*dx;
      for (int j = 0; j <= ny; ++j)
      {
        x.y = bbox->y1 + j*dy;
        for (int k = 0; k <= nz; ++k, ++l)
        {
          x.z = bbox->z1 + k*dz;
          point_t y;
          coord_mapping_map_point(mapping, &x, &y);
          x_node[l] = y.x;
          y_node[l] = y.y;
          z_node[l] = y.z;
        }
      }
    }
  }
  else
  {
    coord_type = DB_COLLINEAR;
    for (int i = 0; i <= nx; ++i)
      x_node[i] = bbox->x1 + i*dx;
    for (int j = 0; j <= ny; ++j)
      y_node[j] = bbox->y1 + j*dy;
    for (int k = 0; k <= nz; ++k)
      z_node[k] = bbox->z1 + k*dz;
  }

  real_t* coords[3];
  coords[0] = x_node;
  coords[1] = y_node;
  coords[2] = z_node;

  // Write the patch grid.
  DBoptlist* optlist = DBMakeOptlist(1);
  int one = 1;
  DBAddOption(optlist, DBOPT_HIDE_FROM_GUI, &one);
  DBPutQuadmesh(dbfile, patch_grid_name, coord_names, coords, dimensions, 3,
                SILO_FLOAT_TYPE, coord_type, optlist);
  DBFreeOptlist(optlist);
  polymec_free(x_node);
  polymec_free(y_node);
  polymec_free(z_node);

  silo_file_pop_dir(file);
}

void silo_file_write_unimesh(silo_file_t* file,
                             const char* mesh_name,
                             unimesh_t* mesh,
                             coord_mapping_t* mapping)
{
  START_FUNCTION_TIMER();

  silo_file_push_domain_dir(file);

  DBfile* dbfile = silo_file_dbfile(file);

  // This mesh is really just a grouping of patches, as far as SILO is
  // concerned.
  DBSetDir(dbfile, "/");
  int num_local_patches = unimesh_num_patches(mesh);
  int npx, npy, npz, nx, ny, nz;
  unimesh_get_extents(mesh, &npx, &npy, &npz);
  unimesh_get_patch_size(mesh, &nx, &ny, &nz);

  // Write the bounding box.
  {
    char bbox_name[FILENAME_MAX+1];
    snprintf(bbox_name, FILENAME_MAX, "%s_bbox", mesh_name);
    bbox_t* bbox = unimesh_bbox(mesh);
    silo_file_write_real_array(file, bbox_name, (real_t*)bbox, 6);
  }

  // Write extent, patch size, and periodicity information for this mesh.
  {
    char array_name[FILENAME_MAX+1];
    snprintf(array_name, FILENAME_MAX, "%s_extents", mesh_name);
    int extents[3] = {npx, npy, npz};
    silo_file_write_int_array(file, array_name, extents, 3);
    snprintf(array_name, FILENAME_MAX, "%s_patch_sizes", mesh_name);
    int patch_sizes[3] = {nx, ny, nz};
    silo_file_write_int_array(file, array_name, patch_sizes, 3);
    snprintf(array_name, FILENAME_MAX, "%s_periodicity", mesh_name);
    bool x_periodic, y_periodic, z_periodic;
    unimesh_get_periodicity(mesh, &x_periodic, &y_periodic, &z_periodic);
    int periodicity[3] = {(int)x_periodic, (int)y_periodic, (int)z_periodic};
    silo_file_write_int_array(file, array_name, periodicity, 3);
  }

  char* patch_grid_names[num_local_patches];
  int patch_grid_types[num_local_patches];
  int pos = 0, i, j, k, l = 0;
  int_array_t* patch_indices = int_array_new();
  bbox_t bbox;
  while (unimesh_next_patch(mesh, &pos, &i, &j, &k, &bbox))
  {
    // Write out the grid for the patch itself.
    char patch_grid_name[FILENAME_MAX+1];
    snprintf(patch_grid_name, FILENAME_MAX, "%s_%d_%d_%d", mesh_name, i, j, k);
    patch_grid_names[l] = string_dup(patch_grid_name);
    patch_grid_types[l] = DB_QUAD_RECT;
    write_unimesh_patch_grid(file, patch_grid_names[l], nx, ny, nz,
                             &bbox, mapping);

    // Jot down this (i, j, k) triple.
    int_array_append(patch_indices, i);
    int_array_append(patch_indices, j);
    int_array_append(patch_indices, k);

    ++l;
  }
  ASSERT(l == num_local_patches);

  // Record the indices of the patches in the mesh.
  {
    char array_name[FILENAME_MAX+1];
    snprintf(array_name, FILENAME_MAX, "%s_patch_indices", mesh_name);
    silo_file_write_int_array(file, array_name, patch_indices->data, patch_indices->size);
    int_array_free(patch_indices);
  }

  // Group all the patches together.
  DBPutMultimesh(dbfile, mesh_name, num_local_patches, (const char* const*)patch_grid_names, patch_grid_types, NULL);

  // Clean up.
  for (int p = 0; p < num_local_patches; ++p)
    string_free(patch_grid_names[p]);

  // Add subdomain information for this mesh.
  silo_file_add_subdomain_mesh(file, mesh_name, DB_QUAD_RECT, NULL);

  silo_file_pop_dir(file);

  STOP_FUNCTION_TIMER();
}

static void read_unimesh_extent_and_patch_size(silo_file_t* file,
                                               const char* mesh_name,
                                               int* npx, int* npy, int* npz,
                                               int* nx, int* ny, int* nz)
{
  char array_name[FILENAME_MAX+1];
  snprintf(array_name, FILENAME_MAX, "%s_extents", mesh_name);
  size_t size;
  int* extents = silo_file_read_int_array(file, array_name, &size);
  if (size != 3)
    polymec_error("silo_file_read_unimesh: Invalid extent data.");
  *npx = extents[0];
  *npy = extents[1];
  *npz = extents[2];
  polymec_free(extents);

  snprintf(array_name, FILENAME_MAX, "%s_patch_sizes", mesh_name);
  int* patch_sizes = silo_file_read_int_array(file, array_name, &size);
  if (size != 3)
    polymec_error("silo_file_read_unimesh: Invalid patch size data.");
  *nx = patch_sizes[0];
  *ny = patch_sizes[1];
  *nz = patch_sizes[2];
  polymec_free(patch_sizes);
}

static void read_unimesh_periodicity(silo_file_t* file,
                                     const char* mesh_name,
                                     bool* x_periodic, bool* y_periodic, bool* z_periodic)
{
  char array_name[FILENAME_MAX+1];
  snprintf(array_name, FILENAME_MAX, "%s_periodicity", mesh_name);
  size_t size;
  int* periodicity = silo_file_read_int_array(file, array_name, &size);
  if (size != 3)
    polymec_error("silo_file_read_unimesh: Invalid periodicity data.");
  *x_periodic = (periodicity[0] != 0);
  *y_periodic = (periodicity[1] != 0);
  *z_periodic = (periodicity[2] != 0);
  polymec_free(periodicity);
}

static void read_unimesh_patch_indices(silo_file_t* file,
                                       const char* mesh_name,
                                       int_array_t* i_array,
                                       int_array_t* j_array,
                                       int_array_t* k_array)
{
  // Clear the arrays.
  int_array_clear(i_array);
  int_array_clear(j_array);
  int_array_clear(k_array);

  char array_name[FILENAME_MAX+1];
  snprintf(array_name, FILENAME_MAX, "%s_patch_indices", mesh_name);
  size_t size;
  int* indices = silo_file_read_int_array(file, array_name, &size);
  if ((size % 3) != 0)
    polymec_error("silo_file_read_unimesh: Invalid patch index data.");
  int num_patches = (int)(size/3);
  for (int p = 0; p < num_patches; ++p)
  {
    int_array_append(i_array, indices[3*p]);
    int_array_append(j_array, indices[3*p+1]);
    int_array_append(k_array, indices[3*p+2]);
  }
  polymec_free(indices);
}

unimesh_t* silo_file_read_unimesh(silo_file_t* file,
                                  const char* mesh_name)
{
  START_FUNCTION_TIMER();

  silo_file_push_domain_dir(file);

  // Read the bounding box.
  char bbox_name[FILENAME_MAX+1];
  snprintf(bbox_name, FILENAME_MAX, "%s_bbox", mesh_name);
  bbox_t bbox;
  size_t six;
  real_t* bounds = silo_file_read_real_array(file, bbox_name, &six);
  memcpy(&bbox, bounds, sizeof(real_t) * six);
  polymec_free(bounds);

  // Read the extent/patch size information for the mesh.
  int npx, npy, npz, nx, ny, nz;
  read_unimesh_extent_and_patch_size(file, mesh_name, &npx, &npy, &npz, &nx, &ny, &nz);

  // Read the periodicity information for the mesh.
  bool x_periodic, y_periodic, z_periodic;
  read_unimesh_periodicity(file, mesh_name, &x_periodic, &y_periodic, &z_periodic);

  // Create the mesh.
#if POLYMEC_HAVE_MPI
  MPI_Comm comm = silo_file_comm(file);
#else
  MPI_Comm comm = MPI_COMM_WORLD;
#endif
  unimesh_t* mesh = create_empty_unimesh(comm, &bbox,
                                         npx, npy, npz,
                                         nx, ny, nz,
                                         x_periodic, y_periodic, z_periodic);

  // Fill it with patches whose indices we read from the file.
  int_array_t* i_array = int_array_new();
  int_array_t* j_array = int_array_new();
  int_array_t* k_array = int_array_new();
  read_unimesh_patch_indices(file, mesh_name, i_array, j_array, k_array);
  for (size_t p = 0; p < i_array->size; ++p)
    unimesh_insert_patch(mesh, i_array->data[p], j_array->data[p], k_array->data[p]);
  int_array_free(i_array);
  int_array_free(j_array);
  int_array_free(k_array);

  unimesh_finalize(mesh);

  silo_file_pop_dir(file);

  STOP_FUNCTION_TIMER();
  return mesh;
}

bool silo_file_contains_unimesh(silo_file_t* file,
                                const char* mesh_name)
{
  silo_file_push_domain_dir(file);
  DBfile* dbfile = silo_file_dbfile(file);
  bool exists = (DBInqVarExists(dbfile, mesh_name) &&
                 (DBInqVarType(dbfile, mesh_name) == DB_MULTIMESH));
  if (exists)
  {
    static const char* data_fields[4] =
      {"bbox_real_array",
       "extents_int_array",
       "patch_sizes_int_array",
       "patch_indices_int_array"};
    for (int i = 0; i < 4; ++i)
    {
      char data_name[FILENAME_MAX+1];
      snprintf(data_name, FILENAME_MAX, "%s_%s", mesh_name, data_fields[i]);
      exists = DBInqVarExists(dbfile, data_name);
      if (!exists) break;
    }
  }
  silo_file_pop_dir(file);
  return exists;
}

static void query_unimesh_vector_comps(unimesh_patch_t* patch,
                                       field_metadata_t* md,
                                       coord_mapping_t* mapping,
                                       bool* is_vector_comp,
                                       int* first_vector_comp)
{
  if ((mapping != NULL) && (field_metadata_has_vectors(md)))
  {
    int pos = 0;
    while (field_metadata_next_vector(md, &pos, first_vector_comp))
    {
      is_vector_comp[*first_vector_comp] = true;
      is_vector_comp[*first_vector_comp+1] = true;
      is_vector_comp[*first_vector_comp+2] = true;
    }
  }
  else
  {
    memset(is_vector_comp, 0, sizeof(bool) * patch->nc);
    *first_vector_comp = -1;
  }
}

static void copy_out_unimesh_node_component(unimesh_patch_t* patch,
                                            field_metadata_t* md,
                                            int c,
                                            bbox_t* bbox,
                                            coord_mapping_t* mapping,
                                            real_t* data)
{
  // If we're given a mapping and there are vector fields, we need to treat
  // them specially.
  bool is_vector_comp[patch->nc];
  int first_vector_comp;
  query_unimesh_vector_comps(patch, md, mapping, is_vector_comp, &first_vector_comp);

  // Now copy the data.
  DECLARE_UNIMESH_NODE_ARRAY(a, patch);
  if ((mapping != NULL) && is_vector_comp[c])
  {
    // We need to map this vector field before we write it out.
    int c1 = first_vector_comp,
        c2 = first_vector_comp+1,
        c3 = first_vector_comp+2;
    int which_component = c - c1;
    int l = 0;
    point_t x;
    real_t dx = (bbox->x2 - bbox->x1)/patch->nx,
           dy = (bbox->y2 - bbox->y1)/patch->ny,
           dz = (bbox->z2 - bbox->z1)/patch->nz;
    for (int i = 0; i <= patch->nx; ++i)
    {
      x.x = bbox->x1 + i * dx;
      for (int j = 0; j <= patch->ny; ++j)
      {
        x.y = bbox->y1 + j * dy;
        for (int k = 0; k <= patch->nz; ++k, ++l)
        {
          x.z = bbox->z1 + k * dz;
          vector_t v = {.x = a[i][j][k][c1], a[i][j][k][c2], a[i][j][k][c3]};
          vector_t v1;
          coord_mapping_map_vector(mapping, &x, &v, &v1);
          switch (which_component)
          {
            case 0: data[l] = v1.x; break;
            case 1: data[l] = v1.y; break;
            default: data[l] = v1.z;
          }
        }
      }
    }
  }
  else
  {
    // Copy the field data verbatim.
    int l = 0;
    for (int i = 0; i <= patch->nx; ++i)
      for (int j = 0; j <= patch->ny; ++j)
        for (int k = 0; k <= patch->nz; ++k, ++l)
          data[l] = a[i][j][k][c];
  }
}

static void copy_out_unimesh_xedge_component(unimesh_patch_t* patch,
                                             field_metadata_t* md,
                                             int c,
                                             bbox_t* bbox,
                                             coord_mapping_t* mapping,
                                             real_t* data)
{
  // If we're given a mapping and there are vector fields, we need to treat
  // them specially.
  bool is_vector_comp[patch->nc];
  int first_vector_comp;
  query_unimesh_vector_comps(patch, md, mapping, is_vector_comp, &first_vector_comp);

  // Now copy the data.
  DECLARE_UNIMESH_XEDGE_ARRAY(a, patch);
  int l = 0;
  if ((mapping != NULL) && is_vector_comp[c])
  {
    // We need to map this vector field before we write it out.
    int c1 = first_vector_comp,
        c2 = first_vector_comp+1,
        c3 = first_vector_comp+2;
    int which_component = c - c1;
    point_t x;
    real_t dx = (bbox->x2 - bbox->x1)/patch->nx,
           dy = (bbox->y2 - bbox->y1)/patch->ny,
           dz = (bbox->z2 - bbox->z1)/patch->nz;
    for (int i = 0; i < patch->nx; ++i)
    {
      x.x = bbox->x1 + (i+0.5) * dx;
      for (int j = 0; j <= patch->ny; ++j)
      {
        x.y = bbox->y1 + j * dy;
        for (int k = 0; k <= patch->nz; ++k, ++l)
        {
          x.z = bbox->z1 + k * dz;
          vector_t v = {.x = a[i][j][k][c1], a[i][j][k][c2], a[i][j][k][c3]};
          vector_t v1;
          coord_mapping_map_vector(mapping, &x, &v, &v1);
          switch (which_component)
          {
            case 0: data[l] = v1.x; break;
            case 1: data[l] = v1.y; break;
            default: data[l] = v1.z;
          }
        }
      }
    }
  }
  else
  {
    // Copy the field data verbatim.
    for (int i = 0; i < patch->nx; ++i)
      for (int j = 0; j <= patch->ny; ++j)
        for (int k = 0; k <= patch->nz; ++k, ++l)
          data[l] = a[i][j][k][c];
  }
}

static void copy_out_unimesh_yedge_component(unimesh_patch_t* patch,
                                             field_metadata_t* md,
                                             int c,
                                             bbox_t* bbox,
                                             coord_mapping_t* mapping,
                                             real_t* data)
{
  // If we're given a mapping and there are vector fields, we need to treat
  // them specially.
  bool is_vector_comp[patch->nc];
  int first_vector_comp;
  query_unimesh_vector_comps(patch, md, mapping, is_vector_comp, &first_vector_comp);

  // Now copy the data.
  int l = (patch->nx+1)*(patch->ny+1)*(patch->nz+1);
  DECLARE_UNIMESH_YEDGE_ARRAY(a, patch);
  if ((mapping != NULL) && is_vector_comp[c])
  {
    // We need to map this vector field before we write it out.
    int c1 = first_vector_comp,
        c2 = first_vector_comp+1,
        c3 = first_vector_comp+2;
    int which_component = c - c1;
    point_t x;
    real_t dx = (bbox->x2 - bbox->x1)/patch->nx,
           dy = (bbox->y2 - bbox->y1)/patch->ny,
           dz = (bbox->z2 - bbox->z1)/patch->nz;
    for (int i = 0; i <= patch->nx; ++i)
    {
      x.x = bbox->x1 + i * dx;
      for (int j = 0; j < patch->ny; ++j)
      {
        x.y = bbox->y1 + (j+0.5) * dy;
        for (int k = 0; k <= patch->nz; ++k, ++l)
        {
          x.z = bbox->z1 + k * dz;
          vector_t v = {.x = a[i][j][k][c1], a[i][j][k][c2], a[i][j][k][c3]};
          vector_t v1;
          coord_mapping_map_vector(mapping, &x, &v, &v1);
          switch (which_component)
          {
            case 0: data[l] = v1.x; break;
            case 1: data[l] = v1.y; break;
            default: data[l] = v1.z;
          }
        }
      }
    }
  }
  else
  {
    // Copy the field data verbatim.
    for (int i = 0; i <= patch->nx; ++i)
      for (int j = 0; j < patch->ny; ++j)
        for (int k = 0; k <= patch->nz; ++k, ++l)
          data[l] = a[i][j][k][c];
  }
}

static void copy_out_unimesh_zedge_component(unimesh_patch_t* patch,
                                             field_metadata_t* md,
                                             int c,
                                             bbox_t* bbox,
                                             coord_mapping_t* mapping,
                                             real_t* data)
{
  // If we're given a mapping and there are vector fields, we need to treat
  // them specially.
  bool is_vector_comp[patch->nc];
  int first_vector_comp;
  query_unimesh_vector_comps(patch, md, mapping, is_vector_comp, &first_vector_comp);

  // Now copy the data.
  int l = 2*(patch->nx+1)*(patch->ny+1)*(patch->nz+1);
  DECLARE_UNIMESH_ZEDGE_ARRAY(a, patch);
  if ((mapping != NULL) && is_vector_comp[c])
  {
    // We need to map this vector field before we write it out.
    int c1 = first_vector_comp,
        c2 = first_vector_comp+1,
        c3 = first_vector_comp+2;
    int which_component = c - c1;
    point_t x;
    real_t dx = (bbox->x2 - bbox->x1)/patch->nx,
           dy = (bbox->y2 - bbox->y1)/patch->ny,
           dz = (bbox->z2 - bbox->z1)/patch->nz;
    for (int i = 0; i <= patch->nx; ++i)
    {
      x.x = bbox->x1 + i * dx;
      for (int j = 0; j <= patch->ny; ++j)
      {
        x.y = bbox->y1 + j * dy;
        for (int k = 0; k < patch->nz; ++k, ++l)
        {
          x.z = bbox->z1 + (k+0.5) * dz;
          vector_t v = {.x = a[i][j][k][c1], a[i][j][k][c2], a[i][j][k][c3]};
          vector_t v1;
          coord_mapping_map_vector(mapping, &x, &v, &v1);
          switch (which_component)
          {
            case 0: data[l] = v1.x; break;
            case 1: data[l] = v1.y; break;
            default: data[l] = v1.z;
          }
        }
      }
    }
  }
  else
  {
    // Copy the field data verbatim.
    for (int i = 0; i <= patch->nx; ++i)
      for (int j = 0; j <= patch->ny; ++j)
        for (int k = 0; k < patch->nz; ++k, ++l)
          data[l] = a[i][j][k][c];
  }
}

static void copy_out_unimesh_xface_component(unimesh_patch_t* patch,
                                             field_metadata_t* md,
                                             int c,
                                             bbox_t* bbox,
                                             coord_mapping_t* mapping,
                                             real_t* data)
{
  // If we're given a mapping and there are vector fields, we need to treat
  // them specially.
  bool is_vector_comp[patch->nc];
  int first_vector_comp;
  query_unimesh_vector_comps(patch, md, mapping, is_vector_comp, &first_vector_comp);

  // Now copy the data.
  int l = 0;
  DECLARE_UNIMESH_XFACE_ARRAY(a, patch);
  if ((mapping != NULL) && is_vector_comp[c])
  {
    // We need to map this vector field before we write it out.
    int c1 = first_vector_comp,
        c2 = first_vector_comp+1,
        c3 = first_vector_comp+2;
    int which_component = c - c1;
    point_t x;
    real_t dx = (bbox->x2 - bbox->x1)/patch->nx,
           dy = (bbox->y2 - bbox->y1)/patch->ny,
           dz = (bbox->z2 - bbox->z1)/patch->nz;
    for (int i = 0; i <= patch->nx; ++i)
    {
      x.x = bbox->x1 + i * dx;
      for (int j = 0; j < patch->ny; ++j)
      {
        x.y = bbox->y1 + (j+0.5) * dy;
        for (int k = 0; k < patch->nz; ++k, ++l)
        {
          x.z = bbox->z1 + (k+0.5) * dz;
          vector_t v = {.x = a[i][j][k][c1], a[i][j][k][c2], a[i][j][k][c3]};
          vector_t v1;
          coord_mapping_map_vector(mapping, &x, &v, &v1);
          switch (which_component)
          {
            case 0: data[l] = v1.x; break;
            case 1: data[l] = v1.y; break;
            default: data[l] = v1.z;
          }
        }
      }
    }
  }
  else
  {
    // Copy the field data verbatim.
    for (int i = 0; i <= patch->nx; ++i)
      for (int j = 0; j < patch->ny; ++j)
        for (int k = 0; k < patch->nz; ++k, ++l)
          data[l] = a[i][j][k][c];
  }
}

static void copy_out_unimesh_yface_component(unimesh_patch_t* patch,
                                             field_metadata_t* md,
                                             int c,
                                             bbox_t* bbox,
                                             coord_mapping_t* mapping,
                                             real_t* data)
{
  // If we're given a mapping and there are vector fields, we need to treat
  // them specially.
  bool is_vector_comp[patch->nc];
  int first_vector_comp;
  query_unimesh_vector_comps(patch, md, mapping, is_vector_comp, &first_vector_comp);

  // Now copy the data.
  int l = (patch->nx+1)*(patch->ny+1)*(patch->nz+1);
  DECLARE_UNIMESH_YFACE_ARRAY(a, patch);
  if ((mapping != NULL) && is_vector_comp[c])
  {
    // We need to map this vector field before we write it out.
    int c1 = first_vector_comp,
        c2 = first_vector_comp+1,
        c3 = first_vector_comp+2;
    int which_component = c - c1;
    point_t x;
    real_t dx = (bbox->x2 - bbox->x1)/patch->nx,
           dy = (bbox->y2 - bbox->y1)/patch->ny,
           dz = (bbox->z2 - bbox->z1)/patch->nz;
    for (int i = 0; i < patch->nx; ++i)
    {
      x.x = bbox->x1 + (i+0.5) * dx;
      for (int j = 0; j <= patch->ny; ++j)
      {
        x.y = bbox->y1 + j * dy;
        for (int k = 0; k < patch->nz; ++k, ++l)
        {
          x.z = bbox->z1 + (k+0.5) * dz;
          vector_t v = {.x = a[i][j][k][c1], a[i][j][k][c2], a[i][j][k][c3]};
          vector_t v1;
          coord_mapping_map_vector(mapping, &x, &v, &v1);
          switch (which_component)
          {
            case 0: data[l] = v1.x; break;
            case 1: data[l] = v1.y; break;
            default: data[l] = v1.z;
          }
        }
      }
    }
  }
  else
  {
    // Copy the field data verbatim.
    for (int i = 0; i < patch->nx; ++i)
      for (int j = 0; j <= patch->ny; ++j)
        for (int k = 0; k < patch->nz; ++k, ++l)
          data[l] = a[i][j][k][c];
  }
}

static void copy_out_unimesh_zface_component(unimesh_patch_t* patch,
                                             field_metadata_t* md,
                                             int c,
                                             bbox_t* bbox,
                                             coord_mapping_t* mapping,
                                             real_t* data)
{
  // If we're given a mapping and there are vector fields, we need to treat
  // them specially.
  bool is_vector_comp[patch->nc];
  int first_vector_comp;
  query_unimesh_vector_comps(patch, md, mapping, is_vector_comp, &first_vector_comp);

  // Now copy the data from the patch.
  int l = 2*(patch->nx+1)*(patch->ny+1)*(patch->nz+1);
  DECLARE_UNIMESH_ZFACE_ARRAY(a, patch);
  if ((mapping != NULL) && is_vector_comp[c])
  {
    // We need to map this vector field before we write it out.
    int c1 = first_vector_comp,
        c2 = first_vector_comp+1,
        c3 = first_vector_comp+2;
    int which_component = c - c1;
    point_t x;
    real_t dx = (bbox->x2 - bbox->x1)/patch->nx,
           dy = (bbox->y2 - bbox->y1)/patch->ny,
           dz = (bbox->z2 - bbox->z1)/patch->nz;
    for (int i = 0; i < patch->nx; ++i)
    {
      x.x = bbox->x1 + (i+0.5) * dx;
      for (int j = 0; j < patch->ny; ++j)
      {
        x.y = bbox->y1 + (j+0.5) * dy;
        for (int k = 0; k <= patch->nz; ++k, ++l)
        {
          x.z = bbox->z1 + k * dz;
          vector_t v = {.x = a[i][j][k][c1], a[i][j][k][c2], a[i][j][k][c3]};
          vector_t v1;
          coord_mapping_map_vector(mapping, &x, &v, &v1);
          switch (which_component)
          {
            case 0: data[l] = v1.x; break;
            case 1: data[l] = v1.y; break;
            default: data[l] = v1.z;
          }
        }
      }
    }
  }
  else
  {
    // Copy the field data verbatim.
    for (int i = 0; i < patch->nx; ++i)
      for (int j = 0; j < patch->ny; ++j)
        for (int k = 0; k <= patch->nz; ++k, ++l)
          data[l] = a[i][j][k][c];
  }
}

static void copy_out_unimesh_cell_component(unimesh_patch_t* patch,
                                            field_metadata_t* md,
                                            int c,
                                            bbox_t* bbox,
                                            coord_mapping_t* mapping,
                                            real_t* data)
{
  // If we're given a mapping and there are vector fields, we need to treat
  // them specially.
  bool is_vector_comp[patch->nc];
  int first_vector_comp;
  query_unimesh_vector_comps(patch, md, mapping, is_vector_comp, &first_vector_comp);

  // Now copy the data.
  DECLARE_UNIMESH_CELL_ARRAY(a, patch);
  if ((mapping != NULL) && is_vector_comp[c])
  {
    // We need to map this vector field before we write it out.
    int c1 = first_vector_comp,
        c2 = first_vector_comp+1,
        c3 = first_vector_comp+2;
    int which_component = c - c1;
    int l = 0;
    point_t x;
    real_t dx = (bbox->x2 - bbox->x1)/patch->nx,
           dy = (bbox->y2 - bbox->y1)/patch->ny,
           dz = (bbox->z2 - bbox->z1)/patch->nz;
    for (int i = 1; i <= patch->nx; ++i)
    {
      x.x = bbox->x1 + (i+0.5) * dx;
      for (int j = 1; j <= patch->ny; ++j)
      {
        x.y = bbox->y1 + (j+0.5) * dy;
        for (int k = 1; k <= patch->nz; ++k, ++l)
        {
          x.z = bbox->z1 + (k+0.5) * dz;
          vector_t v = {.x = a[i][j][k][c1], a[i][j][k][c2], a[i][j][k][c3]};
          vector_t v1;
          coord_mapping_map_vector(mapping, &x, &v, &v1);
          switch (which_component)
          {
            case 0: data[l] = v1.x; break;
            case 1: data[l] = v1.y; break;
            default: data[l] = v1.z;
          }
        }
      }
    }
  }
  else
  {
    // Copy the field data verbatim.
    int l = 0;
    for (int i = 1; i <= patch->nx; ++i)
      for (int j = 1; j <= patch->ny; ++j)
        for (int k = 1; k <= patch->nz; ++k, ++l)
          data[l] = a[i][j][k][c];
  }
}

// This guy copies out the other centerings in an edge- or face-centered field,
// based on the centering in the given patch, and sets the ready_to_write
// flag based on whether all centerings (x, y, and z) will be in the data
// array after the patch's data is copied to it. This is very delicate logic
// and so I've stuffed it all into this function for concentrated head scratching.
static void copy_out_other_centerings(silo_file_t* file,
                                      unimesh_patch_t* patch,
                                      const char* field_component_name,
                                      int c,
                                      real_t* data,
                                      bool* ready_to_write)
{
  string_ptr_unordered_map_t* scratch = silo_file_scratch(file);
  char scratch_name[FILENAME_MAX+1];
  *ready_to_write = true;

  // Handle edges.
  if ((patch->centering == UNIMESH_XEDGE) ||
      (patch->centering == UNIMESH_YEDGE) ||
      (patch->centering == UNIMESH_ZEDGE))
  {
    // First we try to copy the centerings other than the one in the
    // patch.
    if (patch->centering != UNIMESH_XEDGE)
    {
      snprintf(scratch_name, FILENAME_MAX, "%s_x", field_component_name);
      real_t** other_p = (real_t**)string_ptr_unordered_map_get(scratch, scratch_name);
      if (other_p != NULL)
        memcpy(data, *other_p, sizeof(real_t) * patch->nx*(patch->ny+1)*(patch->nz+1));
      else
        *ready_to_write = false;
    }
    if (patch->centering != UNIMESH_YEDGE)
    {
      snprintf(scratch_name, FILENAME_MAX, "%s_y", field_component_name);
      real_t** other_p = (real_t**)string_ptr_unordered_map_get(scratch, scratch_name);
      if (other_p != NULL)
      {
        size_t offset = (patch->nx+1)*(patch->ny+1)*(patch->nz+1);
        memcpy(&data[offset], *other_p, sizeof(real_t) * (patch->nx+1)*patch->ny*(patch->nz+1));
      }
      else
        *ready_to_write = false;
    }
    if (patch->centering != UNIMESH_ZEDGE)
    {
      snprintf(scratch_name, FILENAME_MAX, "%s_z", field_component_name);
      real_t** other_p = (real_t**)string_ptr_unordered_map_get(scratch, scratch_name);
      if (other_p != NULL)
      {
        size_t offset = 2*(patch->nx+1)*(patch->ny+1)*(patch->nz+1);
        memcpy(&data[offset], *other_p, sizeof(real_t) * (patch->nx+1)*(patch->ny+1)*patch->nz);
      }
      else
        *ready_to_write = false;
    }

    // Write the patch data to scratch if we're not going to write it
    // immediately to the file.
    if (!(*ready_to_write))
    {
      real_t* this_one;
      if (patch->centering == UNIMESH_XEDGE)
      {
        snprintf(scratch_name, FILENAME_MAX, "%s_x", field_component_name);
        this_one = polymec_malloc(sizeof(real_t) * patch->nx*(patch->ny+1)*(patch->nz+1));
        DECLARE_UNIMESH_XEDGE_ARRAY(a, patch);
        int l = 0;
        for (int i = 0; i < patch->nx; ++i)
          for (int j = 0; j <= patch->ny; ++j)
            for (int k = 0; k <= patch->nz; ++k, ++l)
              this_one[l] = a[i][j][k][c];
      }
      else if (patch->centering == UNIMESH_YEDGE)
      {
        snprintf(scratch_name, FILENAME_MAX, "%s_y", field_component_name);
        this_one = polymec_malloc(sizeof(real_t) * (patch->nx+1)*patch->ny*(patch->nz+1));
        DECLARE_UNIMESH_YEDGE_ARRAY(a, patch);
        int l = 0;
        for (int i = 0; i <= patch->nx; ++i)
          for (int j = 0; j < patch->ny; ++j)
            for (int k = 0; k <= patch->nz; ++k, ++l)
              this_one[l] = a[i][j][k][c];
      }
      else
      {
        snprintf(scratch_name, FILENAME_MAX, "%s_z", field_component_name);
        this_one = polymec_malloc(sizeof(real_t) * (patch->nx+1)*(patch->ny+1)*patch->nz);
        DECLARE_UNIMESH_ZEDGE_ARRAY(a, patch);
        int l = 0;
        for (int i = 0; i <= patch->nx; ++i)
          for (int j = 0; j <= patch->ny; ++j)
            for (int k = 0; k < patch->nz; ++k, ++l)
              this_one[l] = a[i][j][k][c];
      }
      string_ptr_unordered_map_insert_with_kv_dtors(scratch,
                                                    string_dup(scratch_name),
                                                    this_one,
                                                    string_free,
                                                    polymec_free);
    }
  }

  // Handle faces.
  else if ((patch->centering == UNIMESH_XFACE) ||
           (patch->centering == UNIMESH_YFACE) ||
           (patch->centering == UNIMESH_ZFACE))
  {
    // First we try to copy the centerings other than the one in the
    // patch.
    if (patch->centering != UNIMESH_XFACE)
    {
      snprintf(scratch_name, FILENAME_MAX, "%s_x", field_component_name);
      real_t** other_p = (real_t**)string_ptr_unordered_map_get(scratch, scratch_name);
      if (other_p != NULL)
        memcpy(data, *other_p, sizeof(real_t) * (patch->nx+1)*patch->ny*patch->nz);
      else
        *ready_to_write = false;
    }
    if (patch->centering != UNIMESH_YFACE)
    {
      snprintf(scratch_name, FILENAME_MAX, "%s_y", field_component_name);
      real_t** other_p = (real_t**)string_ptr_unordered_map_get(scratch, scratch_name);
      if (other_p != NULL)
      {
        size_t offset = (patch->nx+1)*(patch->ny+1)*(patch->nz+1);
        memcpy(&data[offset], *other_p, sizeof(real_t) * patch->nx*(patch->ny+1)*patch->nz);
      }
      else
        *ready_to_write = false;
    }
    if (patch->centering != UNIMESH_ZFACE)
    {
      snprintf(scratch_name, FILENAME_MAX, "%s_z", field_component_name);
      real_t** other_p = (real_t**)string_ptr_unordered_map_get(scratch, scratch_name);
      if (other_p != NULL)
      {
        size_t offset = 2*(patch->nx+1)*(patch->ny+1)*(patch->nz+1);
        memcpy(&data[offset], *other_p, sizeof(real_t) * patch->nx*patch->ny*(patch->nz+1));
      }
      else
        *ready_to_write = false;
    }

    // Write the patch data to scratch if we're not going to write it
    // immediately to the file.
    if (!(*ready_to_write))
    {
      real_t* this_one;
      if (patch->centering == UNIMESH_XFACE)
      {
        snprintf(scratch_name, FILENAME_MAX, "%s_x", field_component_name);
        this_one = polymec_malloc(sizeof(real_t) * (patch->nx+1)*patch->ny*patch->nz);
        DECLARE_UNIMESH_XFACE_ARRAY(a, patch);
        int l = 0;
        for (int i = 0; i <= patch->nx; ++i)
          for (int j = 0; j < patch->ny; ++j)
            for (int k = 0; k < patch->nz; ++k, ++l)
              this_one[l] = a[i][j][k][c];
      }
      else if (patch->centering == UNIMESH_YFACE)
      {
        snprintf(scratch_name, FILENAME_MAX, "%s_y", field_component_name);
        this_one = polymec_malloc(sizeof(real_t) * patch->nx*(patch->ny+1)*patch->nz);
        DECLARE_UNIMESH_YFACE_ARRAY(a, patch);
        int l = 0;
        for (int i = 0; i < patch->nx; ++i)
          for (int j = 0; j <= patch->ny; ++j)
            for (int k = 0; k < patch->nz; ++k, ++l)
              this_one[l] = a[i][j][k][c];
      }
      else
      {
        snprintf(scratch_name, FILENAME_MAX, "%s_z", field_component_name);
        this_one = polymec_malloc(sizeof(real_t) * patch->nx*patch->ny*(patch->nz+1));
        DECLARE_UNIMESH_ZFACE_ARRAY(a, patch);
        int l = 0;
        for (int i = 0; i < patch->nx; ++i)
          for (int j = 0; j < patch->ny; ++j)
            for (int k = 0; k <= patch->nz; ++k, ++l)
              this_one[l] = a[i][j][k][c];
      }
      string_ptr_unordered_map_insert_with_kv_dtors(scratch,
                                                    string_dup(scratch_name),
                                                    this_one,
                                                    string_free,
                                                    polymec_free);
    }
  }

  // If we're ready to write to the file, we can clean up our scratch data.
  if (*ready_to_write)
  {
    snprintf(scratch_name, FILENAME_MAX, "%s_x", field_component_name);
    string_ptr_unordered_map_delete(scratch, scratch_name);
    snprintf(scratch_name, FILENAME_MAX, "%s_y", field_component_name);
    string_ptr_unordered_map_delete(scratch, scratch_name);
    snprintf(scratch_name, FILENAME_MAX, "%s_z", field_component_name);
    string_ptr_unordered_map_delete(scratch, scratch_name);
  }
}

static void write_unimesh_patch_data(silo_file_t* file,
                                     const char** field_component_names,
                                     const char* patch_grid_name,
                                     unimesh_patch_t* patch,
                                     field_metadata_t* md,
                                     bbox_t* bbox,
                                     coord_mapping_t* mapping)
{
  DBfile* dbfile = silo_file_dbfile(file);

  // Allocate an array to use for writing Silo data.
  int dimensions[3] = {patch->nx+1, patch->ny+1, patch->nz+1};
  if (patch->centering == UNIMESH_CELL)
  {
    dimensions[0] = patch->nx;
    dimensions[1] = patch->ny;
    dimensions[2] = patch->nz;
  }
  size_t data_size;
  if ((patch->centering == UNIMESH_CELL) || (patch->centering == UNIMESH_NODE))
    data_size = dimensions[0] * dimensions[1] * dimensions[2];
  else
    data_size = 3 * dimensions[0] * dimensions[1] * dimensions[2];
  real_t* data = polymec_calloc(data_size, sizeof(real_t));
  int centering;

  // Now write each component.
  for (int c = 0; c < patch->nc; ++c)
  {
    // Copy the data in the component into our array.
    DBoptlist* optlist = optlist_from_metadata(md, c);
    int one = 1;
    DBAddOption(optlist, DBOPT_HIDE_FROM_GUI, &one);
    bool ready_to_write = false;
    if (patch->centering == UNIMESH_NODE)
    {
      centering = DB_NODECENT;
      copy_out_unimesh_node_component(patch, md, c, bbox, mapping, data);
      ready_to_write = true;
    }
    else if ((patch->centering == UNIMESH_XEDGE) ||
             (patch->centering == UNIMESH_YEDGE) ||
             (patch->centering == UNIMESH_ZEDGE))
    {
      centering = DB_EDGECENT;
      copy_out_other_centerings(file, patch, field_component_names[c], c, data, &ready_to_write);
      if (patch->centering == UNIMESH_XEDGE)
        copy_out_unimesh_xedge_component(patch, md, c, bbox, mapping, data);
      else if (patch->centering == UNIMESH_YEDGE)
        copy_out_unimesh_yedge_component(patch, md, c, bbox, mapping, data);
      else
        copy_out_unimesh_zedge_component(patch, md, c, bbox, mapping, data);
    }
    else if ((patch->centering == UNIMESH_XFACE) ||
             (patch->centering == UNIMESH_YFACE) ||
             (patch->centering == UNIMESH_ZFACE))
    {
      centering = DB_FACECENT;
      copy_out_other_centerings(file, patch, field_component_names[c], c, data, &ready_to_write);
      if (patch->centering == UNIMESH_XFACE)
        copy_out_unimesh_xface_component(patch, md, c, bbox, mapping, data);
      else if (patch->centering == UNIMESH_YFACE)
        copy_out_unimesh_yface_component(patch, md, c, bbox, mapping, data);
      else
        copy_out_unimesh_zface_component(patch, md, c, bbox, mapping, data);
    }
    else
    {
      ASSERT(patch->centering == UNIMESH_CELL);
      centering = DB_ZONECENT;
      copy_out_unimesh_cell_component(patch, md, c, bbox, mapping, data);
      ready_to_write = true;
    }

    // Write the component to the file if it's ready.
    if (ready_to_write)
    {
      DBPutQuadvar1(dbfile, field_component_names[c], patch_grid_name, data,
                    dimensions, 3, NULL, 0, SILO_FLOAT_TYPE, centering, optlist);
    }
    optlist_free(optlist);
  }

  // Clean up.
  polymec_free(data);
}

void silo_file_write_unimesh_field(silo_file_t* file,
                                   const char* field_name,
                                   const char* mesh_name,
                                   unimesh_field_t* field,
                                   coord_mapping_t* mapping)
{
  START_FUNCTION_TIMER();

  silo_file_push_domain_dir(file);

  int num_local_patches = unimesh_field_num_patches(field);
  int num_components = unimesh_field_num_components(field);

  unimesh_t* mesh = unimesh_field_mesh(field);
  int npx, npy, npz, nx, ny, nz;
  unimesh_get_extents(mesh, &npx, &npy, &npz);
  unimesh_get_patch_size(mesh, &nx, &ny, &nz);

  char* field_names[num_components];
  char* multi_field_names[num_components][num_local_patches];
  int multi_field_types[num_local_patches];

  field_metadata_t* md = unimesh_field_metadata(field);
  char md_name[FILENAME_MAX+1];
  snprintf(md_name, FILENAME_MAX, "%s_%s_md", field_name, mesh_name);
  silo_file_write_field_metadata(file, md_name, md);

  unimesh_patch_t* patch;
  int pos = 0, i, j, k, l = 0;
  bbox_t bbox;
  while (unimesh_field_next_patch(field, &pos, &i, &j, &k, &patch, &bbox))
  {
    // Write out the patch data itself.
    for (int c = 0; c < num_components; ++c)
    {
      char field_comp_name[FILENAME_MAX+1];
      snprintf(field_comp_name, FILENAME_MAX, "%s_%d_%d_%d_%d", field_name, c, i, j, k);
      field_names[c] = string_dup(field_comp_name);
      multi_field_names[c][l] = string_dup(field_comp_name);
    }
    char patch_grid_name[FILENAME_MAX+1];
    snprintf(patch_grid_name, FILENAME_MAX, "%s_%d_%d_%d", mesh_name, i, j, k);
    write_unimesh_patch_data(file, (const char**)field_names, patch_grid_name,
                             patch, md, &bbox, mapping);
    multi_field_types[l] = DB_QUADVAR;
    ++l;

    for (int c = 0; c < num_components; ++c)
      string_free(field_names[c]);
  }
  ASSERT(l == num_local_patches);

  // Finally, place multi-* entries into the Silo file.
  DBfile* dbfile = silo_file_dbfile(file);
  for (int c = 0; c < num_components; ++c)
  {
    // We need to associate this multi-variable with our multi-mesh.
    DBoptlist* optlist = DBMakeOptlist(4);
    DBAddOption(optlist, DBOPT_MMESH_NAME, (void*)mesh_name);
    char field_comp_name[FILENAME_MAX+1];
    snprintf(field_comp_name, FILENAME_MAX, "%s_%d", field_name, c);
    DBPutMultivar(dbfile, field_comp_name, num_local_patches, (const char* const*)multi_field_names[c], multi_field_types, NULL);
    DBFreeOptlist(optlist);

    // Add subdomain information for this component.
    silo_file_add_subdomain_mesh(file, field_comp_name, DB_QUAD_RECT, NULL);
  }

  // Clean up.
  for (int c = 0; c < num_components; ++c)
    for (int p = 0; p < num_local_patches; ++p)
      string_free(multi_field_names[c][p]);

  silo_file_pop_dir(file);

  STOP_FUNCTION_TIMER();
}

static void copy_in_unimesh_node_component(DBquadvar* var,
                                           int c,
                                           unimesh_patch_t* patch)
{
  ASSERT(var->dims[0] == patch->nx + 1);
  ASSERT(var->dims[1] == patch->ny + 1);
  ASSERT(var->dims[2] == patch->nz + 1);
  int l = 0;
  real_t* data = (real_t*)var->vals[0];
  DECLARE_UNIMESH_NODE_ARRAY(a, patch);
  for (int i = 0; i <= patch->nx; ++i)
    for (int j = 0; j <= patch->ny; ++j)
      for (int k = 0; k <= patch->nz; ++k, ++l)
        a[i][j][k][c] = data[l];
}

static void copy_in_unimesh_xedge_component(DBquadvar* var,
                                            int c,
                                            unimesh_patch_t* patch)
{
  ASSERT(var->dims[0] == patch->nx + 1);
  ASSERT(var->dims[1] == patch->ny + 1);
  ASSERT(var->dims[2] == patch->nz + 1);
  int l = 0;
  real_t* data = (real_t*)var->vals[0];
  DECLARE_UNIMESH_XEDGE_ARRAY(a, patch);
  for (int i = 0; i < patch->nx; ++i)
    for (int j = 0; j <= patch->ny; ++j)
      for (int k = 0; k <= patch->nz; ++k, ++l)
        a[i][j][k][c] = data[l];
}

static void copy_in_unimesh_yedge_component(DBquadvar* var,
                                            int c,
                                            unimesh_patch_t* patch)
{
  ASSERT(var->dims[0] == patch->nx + 1);
  ASSERT(var->dims[1] == patch->ny + 1);
  ASSERT(var->dims[2] == patch->nz + 1);
  int l = var->dims[0]*var->dims[1]*var->dims[2];
  real_t* data = (real_t*)var->vals[0];
  DECLARE_UNIMESH_YEDGE_ARRAY(a, patch);
  for (int i = 0; i <= patch->nx; ++i)
    for (int j = 0; j < patch->ny; ++j)
      for (int k = 0; k <= patch->nz; ++k, ++l)
        a[i][j][k][c] = data[l];
}

static void copy_in_unimesh_zedge_component(DBquadvar* var,
                                            int c,
                                            unimesh_patch_t* patch)
{
  ASSERT(var->dims[0] == patch->nx + 1);
  ASSERT(var->dims[1] == patch->ny + 1);
  ASSERT(var->dims[2] == patch->nz + 1);
  int l = 2*var->dims[0]*var->dims[1]*var->dims[2];
  real_t* data = (real_t*)var->vals[0];
  DECLARE_UNIMESH_ZEDGE_ARRAY(a, patch);
  for (int i = 0; i <= patch->nx; ++i)
    for (int j = 0; j <= patch->ny; ++j)
      for (int k = 0; k < patch->nz; ++k, ++l)
        a[i][j][k][c] = data[l];
}

static void copy_in_unimesh_xface_component(DBquadvar* var,
                                            int c,
                                            unimesh_patch_t* patch)
{
  ASSERT(var->dims[0] == patch->nx + 1);
  ASSERT(var->dims[1] == patch->ny + 1);
  ASSERT(var->dims[2] == patch->nz + 1);
  int l = 0;
  real_t* data = (real_t*)var->vals[0];
  DECLARE_UNIMESH_XFACE_ARRAY(a, patch);
  for (int i = 0; i <= patch->nx; ++i)
    for (int j = 0; j < patch->ny; ++j)
      for (int k = 0; k < patch->nz; ++k, ++l)
        a[i][j][k][c] = data[l];
}

static void copy_in_unimesh_yface_component(DBquadvar* var,
                                            int c,
                                            unimesh_patch_t* patch)
{
  ASSERT(var->dims[0] == patch->nx + 1);
  ASSERT(var->dims[1] == patch->ny + 1);
  ASSERT(var->dims[2] == patch->nz + 1);
  int l = var->dims[0]*var->dims[1]*var->dims[2];
  real_t* data = (real_t*)var->vals[0];
  DECLARE_UNIMESH_YFACE_ARRAY(a, patch);
  for (int i = 0; i < patch->nx; ++i)
    for (int j = 0; j <= patch->ny; ++j)
      for (int k = 0; k < patch->nz; ++k, ++l)
        a[i][j][k][c] = data[l];
}

static void copy_in_unimesh_zface_component(DBquadvar* var,
                                            int c,
                                            unimesh_patch_t* patch)
{
  ASSERT(var->dims[0] == patch->nx + 1);
  ASSERT(var->dims[1] == patch->ny + 1);
  ASSERT(var->dims[2] == patch->nz + 1);
  int l = 2*var->dims[0]*var->dims[1]*var->dims[2];
  real_t* data = (real_t*)var->vals[0];
  DECLARE_UNIMESH_ZFACE_ARRAY(a, patch);
  for (int i = 0; i < patch->nx; ++i)
    for (int j = 0; j < patch->ny; ++j)
      for (int k = 0; k <= patch->nz; ++k, ++l)
        a[i][j][k][c] = data[l];
}

static void copy_in_unimesh_cell_component(DBquadvar* var,
                                           int c,
                                           unimesh_patch_t* patch)
{
  ASSERT(var->dims[0] == patch->nx);
  ASSERT(var->dims[1] == patch->ny);
  ASSERT(var->dims[2] == patch->nz);
  int l = 0;
  real_t* data = (real_t*)var->vals[0];
  DECLARE_UNIMESH_CELL_ARRAY(a, patch);
  for (int i = 1; i <= patch->nx; ++i)
    for (int j = 1; j <= patch->ny; ++j)
      for (int k = 1; k <= patch->nz; ++k, ++l)
        a[i][j][k][c] = data[l];
}

static void read_unimesh_patch_data(silo_file_t* file,
                                    const char** field_component_names,
                                    const char* patch_grid_name,
                                    int num_components,
                                    unimesh_patch_t* patch,
                                    field_metadata_t* md)
{
  DBfile* dbfile = silo_file_dbfile(file);

  // Fetch each component from the file.
  for (int c = 0; c < num_components; ++c)
  {
    DBquadvar* var = DBGetQuadvar(dbfile, field_component_names[c]);
    ASSERT(var != NULL);
    ASSERT(var->datatype == SILO_FLOAT_TYPE);
    ASSERT(var->nvals == 1);
    switch (patch->centering)
    {
      // Copy the data in the component into our array.
      case UNIMESH_NODE:
        copy_in_unimesh_node_component(var, c, patch);
        break;
      case UNIMESH_XEDGE:
        copy_in_unimesh_xedge_component(var, c, patch);
        break;
      case UNIMESH_YEDGE:
        copy_in_unimesh_yedge_component(var, c, patch);
        break;
      case UNIMESH_ZEDGE:
        copy_in_unimesh_zedge_component(var, c, patch);
        break;
      case UNIMESH_XFACE:
        copy_in_unimesh_xface_component(var, c, patch);
        break;
      case UNIMESH_YFACE:
        copy_in_unimesh_yface_component(var, c, patch);
        break;
      case UNIMESH_ZFACE:
        copy_in_unimesh_zface_component(var, c, patch);
        break;
      default:
        copy_in_unimesh_cell_component(var, c, patch);
    }

    // Extract metadata.
    field_metadata_set_name(md, c, var->label);
    field_metadata_set_units(md, c, var->units);
    field_metadata_set_conserved(md, c, (var->conserved != 0));
    field_metadata_set_extensive(md, c, (var->extensive != 0));

    DBFreeQuadvar(var);
  }
}

void silo_file_read_unimesh_field(silo_file_t* file,
                                  const char* field_name,
                                  const char* mesh_name,
                                  unimesh_field_t* field)
{
  START_FUNCTION_TIMER();
  silo_file_push_domain_dir(file);

  field_metadata_t* md = unimesh_field_metadata(field);
  char md_name[FILENAME_MAX+1];
  snprintf(md_name, FILENAME_MAX, "%s_%s_md", field_name, mesh_name);
  silo_file_read_field_metadata(file, md_name, md);

  unimesh_patch_t* patch;
  int pos = 0, i, j, k;
  while (unimesh_field_next_patch(field, &pos, &i, &j, &k, &patch, NULL))
  {
    char* field_names[patch->nc];
    for (int c = 0; c < patch->nc; ++c)
    {
      char field_comp_name[FILENAME_MAX+1];
      snprintf(field_comp_name, FILENAME_MAX, "%s_%d_%d_%d_%d", field_name, c, i, j, k);
      field_names[c] = string_dup(field_comp_name);
    }
    read_unimesh_patch_data(file, (const char**)field_names, mesh_name,
                            patch->nc, patch, md);
    for (int c = 0; c < patch->nc; ++c)
      string_free(field_names[c]);
  }
  silo_file_pop_dir(file);
  STOP_FUNCTION_TIMER();
}

bool silo_file_contains_unimesh_field(silo_file_t* file,
                                      const char* field_name,
                                      const char* mesh_name,
                                      unimesh_centering_t centering)
{
  bool result = false;
  if (silo_file_contains_unimesh(file, mesh_name))  // mesh exists...
  {
    // Look for the field's metadata array.
    char md_name[FILENAME_MAX+1];
    snprintf(md_name, FILENAME_MAX, "%s_%s_md", field_name, mesh_name);
    result = silo_file_contains_field_metadata(file, md_name);
  }
  return result;
}

