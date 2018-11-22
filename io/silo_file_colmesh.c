// Copyright (c) 2012-2018, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// This file implements colmesh-related methods for the silo_file class.

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
extern DBoptlist* optlist_from_metadata(silo_field_metadata_t* metadata);
extern void optlist_free(DBoptlist* optlist);
extern void silo_file_add_subdomain_mesh(silo_file_t* file, const char* mesh_name, int silo_mesh_type, DBoptlist* optlist);
extern void silo_file_add_subdomain_field(silo_file_t* file, const char* mesh_name, const char* field_name, int silo_field_type, DBoptlist* optlist);
extern void silo_file_push_domain_dir(silo_file_t* file);
extern void silo_file_pop_dir(silo_file_t* file);
extern string_ptr_unordered_map_t* silo_file_scratch(silo_file_t* file);

static void write_colmesh_chunk_grid(silo_file_t* file,
                                     const char* chunk_grid_name,
                                     colmesh_chunk_t* chunk,
                                     coord_mapping_t* mapping)
{
  // Construct a planar polymesh representing the top face of the chunk, 
  // and write it out.
  planar_polymesh_t* chunk_top = planar_polymesh_new(chunk->num_columns,
                                                     chunk->num_xy_faces,
                                                     chunk->num_xy_nodes);
  memcpy(chunk_top->cell_edge_offsets, chunk->column_xy_face_offsets, sizeof(int) * (chunk->num_columns+1));
  planar_polymesh_reserve_connectivity_storage(chunk_top);
  memcpy(chunk_top->cell_edges, chunk->column_xy_faces, sizeof(int) * chunk->column_xy_face_offsets[chunk->num_columns+1]);
  memcpy(chunk_top->edge_cells, chunk->xy_face_columns, 2 * sizeof(int) * chunk->num_xy_faces);
  memcpy(chunk_top->edge_nodes, chunk->xy_edge_nodes, 2 * sizeof(int) * chunk->num_xy_edges);
  memcpy(chunk_top->nodes, chunk->xy_nodes, sizeof(point2_t) * chunk->num_xy_nodes);
  char pp_name[FILENAME_MAX+1];
  snprintf(pp_name, FILENAME_MAX, "%s_pp", chunk_grid_name);
  silo_file_write_planar_polymesh(file, pp_name, chunk_top);
  planar_polymesh_free(chunk_top);

  // Write chunk metadata.
  char array_name[FILENAME_MAX+1];
  snprintf(array_name, FILENAME_MAX, "%s_sizes", chunk_grid_name);
  int sizes[5] = {chunk->num_columns, chunk->num_z_cells,
                  chunk->num_xy_faces, chunk->num_xy_edges,
                  chunk->num_xy_nodes};
  silo_file_write_int_array(file, array_name, sizes, 5);
  snprintf(array_name, FILENAME_MAX, "%s_endpts", chunk_grid_name);
  real_t endpts[2] = {chunk->z1, chunk->z2};
  silo_file_write_real_array(file, array_name, endpts, 2);
}

void silo_file_write_colmesh(silo_file_t* file, 
                             const char* mesh_name,
                             colmesh_t* mesh,
                             coord_mapping_t* mapping)
{
  START_FUNCTION_TIMER();
  silo_file_push_domain_dir(file);

  // Write z axis information for this mesh.
  {
    real_t z1, z2;
    bool periodic;
    colmesh_get_z_info(mesh, &z1, &z2, &periodic);
    char array_name[FILENAME_MAX+1];
    snprintf(array_name, FILENAME_MAX, "%s_endpts", mesh_name);
    real_t endpts[2] = {z1, z2};
    silo_file_write_real_array(file, array_name, endpts, 2);
    snprintf(array_name, FILENAME_MAX, "%s_periodic", mesh_name);
    int periodic_int = (int)(periodic);
    silo_file_write_int_array(file, array_name, &periodic_int, 1);
  }

  size_t l = 0;
  int pos = 0, xy, z;
  colmesh_chunk_t* chunk;
  int_array_t* chunk_indices = int_array_new();
  while (colmesh_next_chunk(mesh, &pos, &xy, &z, &chunk)) 
  {
    // Write out the grid for the chunk itself.
    char chunk_grid_name[FILENAME_MAX+1];
    snprintf(chunk_grid_name, FILENAME_MAX, "%s_%d", mesh_name, xy);
    write_colmesh_chunk_grid(file, chunk_grid_name, chunk, mapping);

    // Jot down this (xy, z) tuple.
    int_array_append(chunk_indices, xy);
    int_array_append(chunk_indices, z);

    ++l;
  }
  ASSERT(l == colmesh_num_chunks(mesh));

  // Record the indices of the chunks in the mesh.
  {
    char array_name[FILENAME_MAX+1];
    snprintf(array_name, FILENAME_MAX, "%s_chunk_indices", mesh_name);
    silo_file_write_int_array(file, array_name, chunk_indices->data, chunk_indices->size);
    int_array_free(chunk_indices);
  }

  // Write chunk metadata.
  int num_xy_chunks, num_z_chunks, nz_per_chunk;
  colmesh_get_chunk_info(mesh, &num_xy_chunks, &num_z_chunks, &nz_per_chunk);
  {
    char array_name[FILENAME_MAX+1];
    snprintf(array_name, FILENAME_MAX, "%s_chunk_md", mesh_name);
    int chunk_md[3] = {(int)num_xy_chunks, (int)num_z_chunks, (int)nz_per_chunk};
    silo_file_write_int_array(file, array_name, chunk_md, 3);
  }

  silo_file_pop_dir(file);
  STOP_FUNCTION_TIMER();
}

static bool weld_pts(void* context, point2_t* p1, point2_t* p2)
{
  return point2s_coincide(p1, p2);
}

colmesh_t* silo_file_read_colmesh(silo_file_t* file, 
                                  const char* mesh_name)
{
  START_FUNCTION_TIMER();
  silo_file_push_domain_dir(file);

  // Read chunk metadata.
  int num_xy_chunks, num_z_chunks, nz_per_chunk;
  {
    char array_name[FILENAME_MAX+1];
    snprintf(array_name, FILENAME_MAX, "%s_chunk_md", mesh_name);
    size_t size;
    int* chunk_md = silo_file_read_int_array(file, array_name, &size);
    ASSERT(size == 3);
    num_xy_chunks = chunk_md[0];
    num_z_chunks = chunk_md[1];
    nz_per_chunk = chunk_md[2];
    polymec_free(chunk_md);
  }

  // Read in the indices of the locally stored chunks:
  // (xy0, z0), (xy1, z1), ...
  size_t num_chunk_indices;
  int* chunk_indices;
  {
    char array_name[FILENAME_MAX+1];
    snprintf(array_name, FILENAME_MAX, "%s_chunk_indices", mesh_name);
    chunk_indices = silo_file_read_int_array(file, array_name, 
                                             &num_chunk_indices);
  }

  // Read the mesh's column data.
  int num_fragments = (int)(num_chunk_indices/2);
  planar_polymesh_t* fragments[num_fragments];
  for (int f = 0; f < num_fragments; ++f)
  {
    int xy = chunk_indices[f/2];
    char fragment_name[FILENAME_MAX+1];
    snprintf(fragment_name, FILENAME_MAX, "%s_%d_pp", mesh_name, xy);
    fragments[f] = silo_file_read_planar_polymesh(file, fragment_name);
  }

  // Read z axis information for this mesh.
  real_t z1, z2;
  bool periodic;
  {
    char array_name[FILENAME_MAX+1];
    snprintf(array_name, FILENAME_MAX, "%s_endpts", mesh_name);
    size_t size;
    real_t* endpts = silo_file_read_real_array(file, array_name, &size);
    ASSERT(size == 2);
    z1 = endpts[0];
    z2 = endpts[1];
    polymec_free(endpts);

    snprintf(array_name, FILENAME_MAX, "%s_periodic", mesh_name);
    int* per = silo_file_read_int_array(file, array_name, &size);
    ASSERT(size == 1);
    periodic = (bool)per[0];
    polymec_free(per);
  }

  // We'll weld fragments together if the points exactly coincide.
  point2_inspector_t* inspector = point2_inspector_new(NULL, weld_pts, NULL);

  // Create the mesh.
#if POLYMEC_HAVE_MPI
  MPI_Comm comm = silo_file_comm(file);
#else
  MPI_Comm comm = MPI_COMM_WORLD;
#endif
  colmesh_t* mesh = 
    create_empty_colmesh_from_fragments(comm, 
                                        fragments, (size_t)num_fragments, 
                                        z1, z2, num_xy_chunks, num_z_chunks, 
                                        nz_per_chunk, periodic, inspector);

  // Insert our local chunks.
  for (size_t i = 0; i < num_chunk_indices/2; ++i)
  {
    int xy_index = 2*i;
    int z_index  = 2*i+1;
    colmesh_insert_chunk(mesh, xy_index, z_index);
  }

  // Clean up.
  release_ref(inspector);
  for (int f = 0; f < num_fragments; ++f)
    planar_polymesh_free(fragments[f]);
  polymec_free(chunk_indices);

  // Finish constructing the colmesh.
  colmesh_finalize(mesh);

  silo_file_pop_dir(file);
  STOP_FUNCTION_TIMER();
  return mesh;
}

bool silo_file_contains_colmesh(silo_file_t* file, 
                                const char* mesh_name)
{
  silo_file_push_domain_dir(file);
  DBfile* dbfile = silo_file_dbfile(file);
  bool exists = (DBInqVarExists(dbfile, mesh_name) && 
                 (DBInqVarType(dbfile, mesh_name) == DB_MULTIMESH));
  if (exists)
  {
    // Check for the column data.
    char columns_name[FILENAME_MAX+1];
    snprintf(columns_name, FILENAME_MAX, "%s_columns", mesh_name);
    exists = silo_file_contains_planar_polymesh(file, columns_name);
    if (exists)
    {
      static const char* data_fields[2] = 
        {"chunk_sizes_int_array", 
         "chunk_indices_int_array"};
      for (int i = 0; i < 2; ++i)
      {
        char data_name[FILENAME_MAX+1];
        snprintf(data_name, FILENAME_MAX, "%s_%s", mesh_name, data_fields[i]);
        exists = DBInqVarExists(dbfile, data_name);
        if (!exists) break;
      }
    }
  }
  silo_file_pop_dir(file);
  return exists;
}

static void query_colmesh_vector_comps(colmesh_chunk_data_t* chunk_data,
                                       silo_field_metadata_t** field_metadata,
                                       coord_mapping_t* mapping,
                                       bool* is_vector_comp,
                                       int* first_vector_comp)
{
  int num_vector_comps = 0;
  *first_vector_comp = -1;
  if ((mapping != NULL) && (field_metadata != NULL))
  {
    for (int c = 0; c < chunk_data->num_components; ++c)
    {
      if ((field_metadata[c] != NULL) && 
          (field_metadata[c]->vector_component != -1))
      {
        if ((num_vector_comps % 3) == 0)
        {
          ASSERT(chunk_data->num_components >= c + 3); // If you start a vector, you better finish your vector.
          *first_vector_comp = c;
        }
        is_vector_comp[c] = true;
        ++num_vector_comps;
      }
    }
  }
  else
    memset(is_vector_comp, 0, sizeof(bool) * chunk_data->num_components);
  ASSERT((num_vector_comps % 3) == 0); // We should have complete sets of vector triples.
}

static void copy_out_colmesh_node_component(colmesh_chunk_data_t* chunk_data,
                                            silo_field_metadata_t** field_metadata,
                                            int c,
                                            coord_mapping_t* mapping,
                                            real_t* data)
{
  // If we're given a mapping and there are vector fields, we need to treat 
  // them specially.
  bool is_vector_comp[chunk_data->num_components];
  int first_vector_comp;
  query_colmesh_vector_comps(chunk_data, field_metadata, mapping,
                              is_vector_comp, &first_vector_comp);

  // Now copy the data.
  DECLARE_COLMESH_NODE_ARRAY(a, chunk_data);
  if ((mapping != NULL) && is_vector_comp[c])
  {
    // We need to map this vector field before we write it out.
    colmesh_chunk_t* chunk = chunk_data->chunk;
    int c1 = first_vector_comp, 
        c2 = first_vector_comp+1, 
        c3 = first_vector_comp+2;
    int which_component = c - c1;
    int l = 0;
    real_t dz = (chunk->z2 - chunk->z1) / chunk->num_z_cells;
    for (int xy = 0; xy < chunk_data->xy_size; ++xy)
    {
      point_t x;
      x.x = chunk->xy_nodes[xy].x;
      x.y = chunk->xy_nodes[xy].y;
      for (int z = 0; z < chunk_data->z_size; ++z, ++l)
      {
        x.z = chunk->z1 + z * dz;
        vector_t v = {.x = a[xy][z][c1], a[xy][z][c2], a[xy][z][c3]};
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
  else
  {
    // Copy the field data verbatim.
    int l = 0;
    for (int xy = 0; xy < chunk_data->xy_size; ++xy)
      for (int z = 0; z < chunk_data->z_size; ++z, ++l)
        data[l] = a[xy][z][c];
  }
}

static void copy_out_colmesh_xyedge_component(colmesh_chunk_data_t* chunk_data,
                                              silo_field_metadata_t** field_metadata,
                                              int c,
                                              coord_mapping_t* mapping,
                                              real_t* data)
{
  // If we're given a mapping and there are vector fields, we need to treat 
  // them specially.
  bool is_vector_comp[chunk_data->num_components];
  int first_vector_comp;
  query_colmesh_vector_comps(chunk_data, field_metadata, mapping,
                              is_vector_comp, &first_vector_comp);

  // Now copy the data.
  DECLARE_COLMESH_XYEDGE_ARRAY(a, chunk_data);
  int l = 0;
  if ((mapping != NULL) && is_vector_comp[c])
  {
    // We need to map this vector field before we write it out.
    colmesh_chunk_t* chunk = chunk_data->chunk;
    int c1 = first_vector_comp, 
        c2 = first_vector_comp+1, 
        c3 = first_vector_comp+2;
    int which_component = c - c1;
    real_t dz = (chunk->z2 - chunk->z1) / chunk->num_z_cells;
    for (int xy = 0; xy < chunk_data->xy_size; ++xy)
    {
      point_t x;
      int n1 = chunk->xy_edge_nodes[2*xy];
      int n2 = chunk->xy_edge_nodes[2*xy+1];
      x.x = 0.5 * (chunk->xy_nodes[n1].x + chunk->xy_nodes[n2].x);
      x.y = 0.5 * (chunk->xy_nodes[n1].y + chunk->xy_nodes[n2].y);
      for (int z = 0; z < chunk_data->z_size; ++z, ++l)
      {
        x.z = chunk->z1 + z * dz;
        vector_t v = {.x = a[xy][z][c1], a[xy][z][c2], a[xy][z][c3]};
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
  else
  {
    // Copy the field data verbatim.
    for (int xy = 0; xy < chunk_data->xy_size; ++xy)
      for (int z = 0; z < chunk_data->z_size; ++z, ++l)
        data[l] = a[xy][z][c];
  }
}

static void copy_out_colmesh_zedge_component(colmesh_chunk_data_t* chunk_data,
                                             silo_field_metadata_t** field_metadata,
                                             int c,
                                             coord_mapping_t* mapping,
                                             real_t* data)
{
  // If we're given a mapping and there are vector fields, we need to treat 
  // them specially.
  bool is_vector_comp[chunk_data->num_components];
  int first_vector_comp;
  query_colmesh_vector_comps(chunk_data, field_metadata, mapping,
                              is_vector_comp, &first_vector_comp);

  // Now copy the data.
  colmesh_chunk_t* chunk = chunk_data->chunk;
  size_t l = chunk->num_xy_edges;
  DECLARE_COLMESH_ZEDGE_ARRAY(a, chunk_data);
  if ((mapping != NULL) && is_vector_comp[c])
  {
    // We need to map this vector field before we write it out.
    int c1 = first_vector_comp, 
        c2 = first_vector_comp+1, 
        c3 = first_vector_comp+2;
    int which_component = c - c1;
    real_t dz = (chunk->z2 - chunk->z1) / chunk->num_z_cells;
    for (int xy = 0; xy < chunk_data->xy_size; ++xy)
    {
      point_t x;
      x.x = chunk->xy_nodes[xy].x;
      x.y = chunk->xy_nodes[xy].y;
      for (int z = 0; z < chunk_data->z_size; ++z, ++l)
      {
        x.z = chunk->z1 + (z+0.5) * dz;
        vector_t v = {.x = a[xy][z][c1], a[xy][z][c2], a[xy][z][c3]};
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
  else
  {
    // Copy the field data verbatim.
    for (int xy = 0; xy < chunk_data->xy_size; ++xy)
      for (int z = 0; z < chunk_data->z_size; ++z, ++l)
        data[l] = a[xy][z][c];
  }
}

static void copy_out_colmesh_xyface_component(colmesh_chunk_data_t* chunk_data,
                                              silo_field_metadata_t** field_metadata,
                                              int c,
                                              coord_mapping_t* mapping,
                                              real_t* data)
{
  // If we're given a mapping and there are vector fields, we need to treat 
  // them specially.
  bool is_vector_comp[chunk_data->num_components];
  int first_vector_comp;
  query_colmesh_vector_comps(chunk_data, field_metadata, mapping,
                              is_vector_comp, &first_vector_comp);

  // Now copy the data.
  int l = 0;
  DECLARE_COLMESH_XYFACE_ARRAY(a, chunk_data);
  if ((mapping != NULL) && is_vector_comp[c])
  {
    // We need to map this vector field before we write it out.
    colmesh_chunk_t* chunk = chunk_data->chunk;
    int c1 = first_vector_comp, 
        c2 = first_vector_comp+1, 
        c3 = first_vector_comp+2;
    int which_component = c - c1;
    real_t dz = (chunk->z2 - chunk->z1) / chunk->num_z_cells;
    for (int xy = 0; xy < chunk_data->xy_size; ++xy)
    {
      // The face's xy coordinates are an average of the attached nodes' xy coordinates.
      // Recall that an xy face shares an index with its corresponding xy edge, for 
      // a fixed z position. So we can use the same technique we use for computing 
      // an "edge center" for the x and y coordinates.
      point_t x;
      int n1 = chunk->xy_edge_nodes[2*xy];
      int n2 = chunk->xy_edge_nodes[2*xy+1];
      x.x = 0.5 * (chunk->xy_nodes[n1].x + chunk->xy_nodes[n2].x);
      x.y = 0.5 * (chunk->xy_nodes[n1].y + chunk->xy_nodes[n2].y);
      for (int z = 0; z < chunk_data->z_size; ++z, ++l)
      {
        x.z = chunk->z1 + (z+0.5) * dz;
        vector_t v = {.x = a[xy][z][c1], a[xy][z][c2], a[xy][z][c3]};
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
  else
  {
    // Copy the field data verbatim.
    for (int xy = 0; xy < chunk_data->xy_size; ++xy)
      for (int z = 0; z < chunk_data->z_size; ++z, ++l)
        data[l] = a[xy][z][c];
  }
}

static void copy_out_colmesh_zface_component(colmesh_chunk_data_t* chunk_data,
                                             silo_field_metadata_t** field_metadata,
                                             int c,
                                             coord_mapping_t* mapping,
                                             real_t* data)
{
  // If we're given a mapping and there are vector fields, we need to treat 
  // them specially.
  bool is_vector_comp[chunk_data->num_components];
  int first_vector_comp;
  query_colmesh_vector_comps(chunk_data, field_metadata, mapping,
                              is_vector_comp, &first_vector_comp);

  // Now copy the data from the chunk.
  colmesh_chunk_t* chunk = chunk_data->chunk;
  size_t l = chunk->num_xy_faces;
  DECLARE_COLMESH_ZFACE_ARRAY(a, chunk_data);
  if ((mapping != NULL) && is_vector_comp[c])
  {
    // We need to map this vector field before we write it out.
    int c1 = first_vector_comp, 
        c2 = first_vector_comp+1, 
        c3 = first_vector_comp+2;
    int which_component = c - c1;
    real_t dz = (chunk->z2 - chunk->z1) / chunk->num_z_cells;
    for (int xy = 0; xy < chunk_data->xy_size; ++xy)
    {
      // The xy coordinates of the z face are those for the cell, which are the 
      // geometric average of its node positions.
      point_t x = {.x = 0.0, .y = 0.0, .z = 0.0};
      int num_nodes = colmesh_chunk_z_face_num_nodes(chunk, xy);
      int nodes[num_nodes];
      colmesh_chunk_z_face_get_nodes(chunk, xy, nodes);
      for (int n = 0; n < num_nodes; ++n)
      {
        x.x += chunk->xy_nodes[nodes[n]].x;
        x.y += chunk->xy_nodes[nodes[n]].y;
      }
      x.x /= num_nodes;
      x.y /= num_nodes;

      for (int z = 0; z < chunk_data->z_size; ++z, ++l)
      {
        x.z = chunk->z1 + z * dz;
        vector_t v = {.x = a[xy][z][c1], a[xy][z][c2], a[xy][z][c3]};
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
  else
  {
    // Copy the field data verbatim.
    for (int xy = 0; xy < chunk_data->xy_size; ++xy)
      for (int z = 0; z < chunk_data->z_size; ++z, ++l)
        data[l] = a[xy][z][c];
  }
}

static void copy_out_colmesh_cell_component(colmesh_chunk_data_t* chunk_data,
                                            silo_field_metadata_t** field_metadata,
                                            int c,
                                            coord_mapping_t* mapping,
                                            real_t* data)
{
  // If we're given a mapping and there are vector fields, we need to treat 
  // them specially.
  bool is_vector_comp[chunk_data->num_components];
  int first_vector_comp;
  query_colmesh_vector_comps(chunk_data, field_metadata, mapping,
                             is_vector_comp, &first_vector_comp);

  // Now copy the data.
  colmesh_chunk_t* chunk = chunk_data->chunk;
  DECLARE_COLMESH_CELL_ARRAY(a, chunk_data);
  if ((mapping != NULL) && is_vector_comp[c])
  {
    // We need to map this vector field before we write it out.
    int c1 = first_vector_comp, 
        c2 = first_vector_comp+1, 
        c3 = first_vector_comp+2;
    int which_component = c - c1;
    int l = 0;
    real_t dz = (chunk->z2 - chunk->z1) / chunk->num_z_cells;
    for (int xy = 0; xy < chunk_data->xy_size; ++xy)
    {
      // Compute the centroid of the cell in the xy plane.
      point_t x = {.x = 0.0, .y = 0.0, .z = 0.0};
      int num_nodes = colmesh_chunk_z_face_num_nodes(chunk, xy);
      int nodes[num_nodes];
      colmesh_chunk_z_face_get_nodes(chunk, xy, nodes);
      for (int n = 0; n < num_nodes; ++n)
      {
        x.x += chunk->xy_nodes[nodes[n]].x;
        x.y += chunk->xy_nodes[nodes[n]].y;
      }
      x.x /= num_nodes;
      x.y /= num_nodes;

      for (int z = 1; z <= chunk->num_z_cells; ++z, ++l)
      {
        x.z = chunk->z1 + (z+0.5) * dz;
        vector_t v = {.x = a[xy][z][c1], a[xy][z][c2], a[xy][z][c3]};
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
  else
  {
    // Copy the field data verbatim.
    int l = 0;
    for (int xy = 0; xy < chunk_data->xy_size; ++xy)
      for (int z = 1; z <= chunk->num_z_cells; ++z, ++l)
        data[l] = a[xy][z][c];
  }
}

// This guy copies out the other centerings in an edge- or face-centered field, 
// based on the centering in the given chunk, and sets the ready_to_write 
// flag based on whether all centerings (x, y, and z) will be in the data
// array after the chunk's data is copied to it. This is very delicate logic 
// and so I've stuffed it all into this function for concentrated head scratching.
static void copy_out_other_centerings(silo_file_t* file,
                                      colmesh_chunk_data_t* chunk_data,
                                      const char* field_component_name,
                                      int c,
                                      real_t* data,
                                      bool* ready_to_write)
{
  string_ptr_unordered_map_t* scratch = silo_file_scratch(file);
  colmesh_chunk_t* chunk = chunk_data->chunk;
  char scratch_name[FILENAME_MAX+1];
  *ready_to_write = true;

  // Handle edges.
  if ((chunk_data->centering == COLMESH_XYEDGE) ||
      (chunk_data->centering == COLMESH_ZEDGE))
  {
    // First we try to copy the centerings other than the one in the 
    // chunk.
    if (chunk_data->centering != COLMESH_XYEDGE)
    {
      snprintf(scratch_name, FILENAME_MAX, "%s_xy", field_component_name);
      real_t** other_p = (real_t**)string_ptr_unordered_map_get(scratch, scratch_name);
      if (other_p != NULL)
        memcpy(data, *other_p, sizeof(real_t) * chunk_data->xy_size * chunk_data->z_size);
      else
        *ready_to_write = false;
    }
    if (chunk_data->centering != COLMESH_ZEDGE)
    {
      snprintf(scratch_name, FILENAME_MAX, "%s_z", field_component_name);
      real_t** other_p = (real_t**)string_ptr_unordered_map_get(scratch, scratch_name);
      if (other_p != NULL)
      {
        size_t offset = chunk->num_xy_edges;
        memcpy(&data[offset], *other_p, sizeof(real_t) * chunk_data->xy_size * chunk_data->z_size);
      }
      else
        *ready_to_write = false;
    }

    // Write the chunk data to scratch if we're not going to write it 
    // immediately to the file.
    if (!(*ready_to_write))
    {
      real_t* this_one = polymec_malloc(sizeof(real_t) * chunk_data->xy_size * chunk_data->z_size);
      if (chunk_data->centering == COLMESH_XYEDGE)
      {
        snprintf(scratch_name, FILENAME_MAX, "%s_xy", field_component_name);
        DECLARE_COLMESH_XYEDGE_ARRAY(a, chunk_data);
        int l = 0;
        for (int xy = 0; xy < chunk_data->xy_size; ++xy)
          for (int z = 0; z < chunk_data->z_size; ++z, ++l)
            this_one[l] = a[xy][z][c];
      }
      else
      {
        snprintf(scratch_name, FILENAME_MAX, "%s_z", field_component_name);
        DECLARE_COLMESH_ZEDGE_ARRAY(a, chunk_data);
        int l = 0;
        for (int xy = 0; xy < chunk_data->xy_size; ++xy)
          for (int z = 0; z < chunk_data->z_size; ++z, ++l)
            this_one[l] = a[xy][z][c];
      }
      string_ptr_unordered_map_insert_with_kv_dtors(scratch, 
                                                    string_dup(scratch_name),
                                                    this_one,
                                                    string_free,
                                                    polymec_free);
    }
  }

  // Handle faces.
  else if ((chunk_data->centering == COLMESH_XYFACE) ||
           (chunk_data->centering == COLMESH_ZFACE))
  {
    // First we try to copy the centerings other than the one in the 
    // chunk.
    if (chunk_data->centering != COLMESH_XYFACE)
    {
      snprintf(scratch_name, FILENAME_MAX, "%s_xy", field_component_name);
      real_t** other_p = (real_t**)string_ptr_unordered_map_get(scratch, scratch_name);
      if (other_p != NULL)
        memcpy(data, *other_p, sizeof(real_t) * chunk_data->xy_size * chunk_data->z_size);
      else
        *ready_to_write = false;
    }
    if (chunk_data->centering != COLMESH_ZFACE)
    {
      snprintf(scratch_name, FILENAME_MAX, "%s_z", field_component_name);
      real_t** other_p = (real_t**)string_ptr_unordered_map_get(scratch, scratch_name);
      if (other_p != NULL)
      {
        size_t offset = chunk->num_xy_faces;
        memcpy(&data[offset], *other_p, sizeof(real_t) * chunk_data->xy_size * chunk_data->z_size);
      }
      else
        *ready_to_write = false;
    }

    // Write the chunk data to scratch if we're not going to write it 
    // immediately to the file.
    if (!(*ready_to_write))
    {
      real_t* this_one = polymec_malloc(sizeof(real_t) * chunk_data->xy_size * chunk_data->z_size);
      if (chunk_data->centering == COLMESH_XYFACE)
      {
        snprintf(scratch_name, FILENAME_MAX, "%s_xy", field_component_name);
        DECLARE_COLMESH_XYFACE_ARRAY(a, chunk_data);
        int l = 0;
        for (int xy = 0; xy < chunk_data->xy_size; ++xy)
          for (int z = 0; z < chunk_data->z_size; ++z)
            this_one[l] = a[xy][z][c];
      }
      else
      {
        snprintf(scratch_name, FILENAME_MAX, "%s_z", field_component_name);
        DECLARE_COLMESH_ZFACE_ARRAY(a, chunk_data);
        int l = 0;
        for (int xy = 0; xy < chunk_data->xy_size; ++xy)
          for (int z = 0; z < chunk_data->z_size; ++z)
            this_one[l] = a[xy][z][c];
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
    snprintf(scratch_name, FILENAME_MAX, "%s_xy", field_component_name);
    string_ptr_unordered_map_delete(scratch, scratch_name);
    snprintf(scratch_name, FILENAME_MAX, "%s_z", field_component_name);
    string_ptr_unordered_map_delete(scratch, scratch_name);
  }
}

// Names associated with field centerings.
static const char* centering_name[6] = {"cell", "xyface", "zface", "xyedge", "zedge", "node"};

static void write_colmesh_chunk_data(silo_file_t* file,
                                     const char** field_component_names,
                                     const char* chunk_grid_name,
                                     colmesh_chunk_data_t* chunk_data,
                                     silo_field_metadata_t** field_metadata,
                                     coord_mapping_t* mapping)
{
  // Because we can't really represent a colmesh faithfully in a SILO format,
  // we're really just dumping a bunch of data into an array.
  size_t data_size;
  if ((chunk_data->centering == COLMESH_CELL) || 
      (chunk_data->centering == COLMESH_NODE))
    data_size = chunk_data->xy_size * chunk_data->z_size;
  else 
    data_size = 2 * chunk_data->xy_size * chunk_data->z_size;
  real_t* data = polymec_calloc(sizeof(real_t) * data_size);

  // Now write each component.
  for (int c = 0; c < chunk_data->num_components; ++c)
  {
    // Copy the data in the component into our array.
    bool ready_to_write = false;
    if (chunk_data->centering == COLMESH_NODE) 
    {
      copy_out_colmesh_node_component(chunk_data, field_metadata, c, mapping, data);
      ready_to_write = true;
    }
    else if ((chunk_data->centering == COLMESH_XYEDGE) || 
             (chunk_data->centering == COLMESH_ZEDGE))
    {
      copy_out_other_centerings(file, chunk_data, field_component_names[c], c, data, &ready_to_write);
      if (chunk_data->centering == COLMESH_XYEDGE)
        copy_out_colmesh_xyedge_component(chunk_data, field_metadata, c, mapping, data);
      else // (data->centering == COLMESH_ZEDGE)
        copy_out_colmesh_zedge_component(chunk_data, field_metadata, c, mapping, data);
    }
    else if ((chunk_data->centering == COLMESH_XYFACE) || 
             (chunk_data->centering == COLMESH_ZFACE))
    {
      copy_out_other_centerings(file, chunk_data, field_component_names[c], c, data, &ready_to_write);
      if (chunk_data->centering == COLMESH_XYFACE)
        copy_out_colmesh_xyface_component(chunk_data, field_metadata, c, mapping, data);
      else // (data->centering == COLMESH_ZFACE)
        copy_out_colmesh_zface_component(chunk_data, field_metadata, c, mapping, data);
    }
    else
    {
      ASSERT(chunk_data->centering == COLMESH_CELL);
      copy_out_colmesh_cell_component(chunk_data, field_metadata, c, mapping, data);
      ready_to_write = true;
    }
    
    // Write the component to the file if it's ready.
    if (ready_to_write)
    {
      char data_name[FILENAME_MAX+1];
      snprintf(data_name, FILENAME_MAX, "%s_%s_%s", chunk_grid_name, field_component_names[c], centering_name[(int)chunk_data->centering]);
      silo_file_write_real_array(file, data_name, data, data_size);

      // Pack any metadata into an array.
      if ((field_metadata != NULL) && (field_metadata[c] != NULL))
      {
        char md_name[FILENAME_MAX+1];
        snprintf(md_name, FILENAME_MAX, "%s_md", data_name);
        size_t label_size = strlen(field_metadata[c]->label);
        size_t units_size = strlen(field_metadata[c]->units);
        size_t md_size = 5 * sizeof(int) + label_size + units_size;
        int md[md_size];
        md[0] = (int)label_size; 
        for (size_t i = 0; i < label_size; ++i)
          md[1+i] = (int)(field_metadata[c]->label[i]);
        md[label_size+1] = (int)units_size; 
        for (size_t i = 0; i < units_size; ++i)
          md[1+label_size+1+i] = (int)(field_metadata[c]->units[i]);
        md[1+label_size+1+units_size] = (int)(field_metadata[c]->conserved);
        md[1+label_size+1+units_size+1] = (int)(field_metadata[c]->extensive);
        md[1+label_size+1+units_size+2] = (int)(field_metadata[c]->vector_component);
        silo_file_write_int_array(file, md_name, md, md_size);
      }
    }
  }

  // Clean up.
  polymec_free(data);
}

void silo_file_write_colmesh_field(silo_file_t* file, 
                                   const char** field_component_names,
                                   const char* mesh_name,
                                   colmesh_field_t* field,
                                   silo_field_metadata_t** field_metadata,
                                   coord_mapping_t* mapping)
{
  START_FUNCTION_TIMER();

  silo_file_push_domain_dir(file);

  size_t num_components = colmesh_field_num_components(field);

  char* field_names[num_components];
  colmesh_chunk_data_t* data;
  int pos = 0, xy, z, l = 0;
  while (colmesh_field_next_chunk(field, &pos, &xy, &z, &data))
  {
    // Write out the chunk data itself.
    for (int c = 0; c < num_components; ++c)
    {
      char field_name[FILENAME_MAX];
      snprintf(field_name, FILENAME_MAX-1, "%s_%d_%d", field_component_names[c], xy, z);
      field_names[c] = string_dup(field_name);
    }
    char chunk_grid_name[FILENAME_MAX];
    snprintf(chunk_grid_name, FILENAME_MAX-1, "%s_%d_%d", mesh_name, xy, z);
    write_colmesh_chunk_data(file, (const char**)field_names, chunk_grid_name,  
                             data, field_metadata, mapping);
    ++l;

    for (int c = 0; c < num_components; ++c)
      string_free(field_names[c]);
  }
  ASSERT(l == colmesh_field_num_chunks(field));

  silo_file_pop_dir(file);

  STOP_FUNCTION_TIMER();
}

static void copy_in_colmesh_node_component(real_t* data, 
                                            int c, 
                                            colmesh_chunk_data_t* chunk_data)
{
  int l = 0;
  DECLARE_COLMESH_NODE_ARRAY(a, chunk_data);
  for (int xy = 0; xy < chunk_data->xy_size; ++xy)
    for (int z = 0; z < chunk_data->z_size; ++z, ++l)
      a[xy][z][c] = data[l];
}

static void copy_in_colmesh_xyedge_component(real_t* data, 
                                             int c, 
                                             colmesh_chunk_data_t* chunk_data)
{
  int l = 0;
  DECLARE_COLMESH_XYEDGE_ARRAY(a, chunk_data);
  for (int xy = 0; xy < chunk_data->xy_size; ++xy)
    for (int z = 0; z < chunk_data->z_size; ++z, ++l)
      a[xy][z][c] = data[l];
}

static void copy_in_colmesh_zedge_component(real_t* data, 
                                            int c, 
                                            colmesh_chunk_data_t* chunk_data)
{
  size_t l = chunk_data->chunk->num_xy_edges;
  DECLARE_COLMESH_ZEDGE_ARRAY(a, chunk_data);
  for (int xy = 0; xy < chunk_data->xy_size; ++xy)
    for (int z = 0; z < chunk_data->z_size; ++z, ++l)
      a[xy][z][c] = data[l];
}

static void copy_in_colmesh_xyface_component(real_t* data, 
                                             int c, 
                                             colmesh_chunk_data_t* chunk_data)
{
  int l = 0;
  DECLARE_COLMESH_XYFACE_ARRAY(a, chunk_data);
  for (int xy = 0; xy < chunk_data->xy_size; ++xy)
    for (int z = 0; z < chunk_data->z_size; ++z, ++l)
      a[xy][z][c] = data[l];
}

static void copy_in_colmesh_zface_component(real_t* data, 
                                            int c, 
                                            colmesh_chunk_data_t* chunk_data)
{
  size_t l = chunk_data->chunk->num_xy_faces;
  DECLARE_COLMESH_ZFACE_ARRAY(a, chunk_data);
  for (int xy = 0; xy < chunk_data->xy_size; ++xy)
    for (int z = 0; z < chunk_data->z_size; ++z, ++l)
      a[xy][z][c] = data[l];
}

static void copy_in_colmesh_cell_component(real_t* data, 
                                           int c, 
                                           colmesh_chunk_data_t* chunk_data)
{
  int l = 0;
  DECLARE_COLMESH_CELL_ARRAY(a, chunk_data);
  for (int xy = 0; xy < chunk_data->xy_size; ++xy)
    for (int z = 1; z <= chunk_data->chunk->num_z_cells; ++z, ++l)
      a[xy][z][c] = data[l];
}

static void read_colmesh_chunk_data(silo_file_t* file,
                                    const char** field_component_names,
                                    const char* chunk_grid_name,
                                    size_t num_components,
                                    colmesh_chunk_data_t* chunk_data,
                                    silo_field_metadata_t** field_metadata)
{
  // Fetch each component from the file.
  for (int c = 0; c < (int)num_components; ++c)
  {
    char data_name[FILENAME_MAX+1];
    snprintf(data_name, FILENAME_MAX, "%s_%s_%s", chunk_grid_name, field_component_names[c], centering_name[(int)chunk_data->centering]);
    size_t data_size;
    real_t* data = silo_file_read_real_array(file, data_name, &data_size);
    ASSERT(data != NULL);
    switch (chunk_data->centering) 
    {
      // Copy the data in the component into our array.
      case COLMESH_NODE: 
        copy_in_colmesh_node_component(data, c, chunk_data);
        break;
      case COLMESH_XYEDGE:
        copy_in_colmesh_xyedge_component(data, c, chunk_data);
        break;
      case COLMESH_ZEDGE:
        copy_in_colmesh_zedge_component(data, c, chunk_data);
        break;
      case COLMESH_XYFACE:
        copy_in_colmesh_xyface_component(data, c, chunk_data);
        break;
      case COLMESH_ZFACE:
        copy_in_colmesh_zface_component(data, c, chunk_data);
        break;
      default:
        copy_in_colmesh_cell_component(data, c, chunk_data);
    }

    // Extract metadata.
    if ((field_metadata != NULL) && (field_metadata[c] != NULL))
    {
      // Unpack the metadata from an array.
      char md_name[FILENAME_MAX+1];
      snprintf(md_name, FILENAME_MAX, "%s_md", data_name);
      size_t md_size;
      int* md = silo_file_read_int_array(file, md_name, &md_size);
      if (md != NULL)
      {
        size_t label_size = (int)md[0];
        field_metadata[c]->label = string_ndup((char*)(&md[1]), label_size);
        size_t units_size = (int)md[1+label_size];
        field_metadata[c]->units = string_ndup((char*)(&md[1+label_size+1]), units_size);
        field_metadata[c]->conserved = (int)(md[1+label_size+1+units_size]);
        field_metadata[c]->extensive = (int)(md[1+label_size+1+units_size+1]);
        field_metadata[c]->vector_component = (int)(md[1+label_size+1+units_size+2]);
      }
    }

    polymec_free(data);
  }
}

void silo_file_read_colmesh_field(silo_file_t* file, 
                                  const char** field_component_names,
                                  const char* mesh_name,
                                  colmesh_field_t* field,
                                  silo_field_metadata_t** field_metadata)
{
  START_FUNCTION_TIMER();
  silo_file_push_domain_dir(file);
  colmesh_chunk_data_t* data;
  int pos = 0, xy, z;
  while (colmesh_field_next_chunk(field, &pos, &xy, &z, &data))
  {
    char* field_names[data->num_components];
    for (int c = 0; c < data->num_components; ++c)
    {
      char field_name[FILENAME_MAX+1];
      snprintf(field_name, FILENAME_MAX, "%s_%d_%d", field_component_names[c], xy, z);
      field_names[c] = string_dup(field_name);
    }
    read_colmesh_chunk_data(file, (const char**)field_names, mesh_name, 
                             data->num_components, data, field_metadata);
    for (int c = 0; c < data->num_components; ++c)
      string_free(field_names[c]);
  }
  silo_file_pop_dir(file);
  STOP_FUNCTION_TIMER();
}

bool silo_file_contains_colmesh_field(silo_file_t* file, 
                                      const char* field_name,
                                      const char* mesh_name,
                                      colmesh_centering_t centering)
{
  // FIXME
  silo_file_push_domain_dir(file);
  DBfile* dbfile = silo_file_dbfile(file);
  bool exists = (DBInqVarExists(dbfile, mesh_name) && 
                 (DBInqVarType(dbfile, mesh_name) == DB_MULTIMESH));
  silo_file_pop_dir(file);
  return exists;
}

