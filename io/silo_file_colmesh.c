// Copyright (c) 2012-2019, Jeffrey N. Johnson
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

extern exchanger_proc_map_t* colmesh_xy_data_send_map(colmesh_t* mesh, int xy_index);
extern exchanger_proc_map_t* colmesh_xy_data_receive_map(colmesh_t* mesh, int xy_index);

static int_array_t* clone_int_array(int_array_t* array)
{
  return int_array_clone(array, NULL);
}

static void write_colmesh_chunk_grid(silo_file_t* file,
                                     const char* chunk_grid_name,
                                     colmesh_t* mesh,
                                     int xy_index, int z_index,
                                     coord_mapping_t* mapping)
{
  colmesh_chunk_t* chunk = colmesh_chunk(mesh, xy_index, z_index);
  // Construct a planar polymesh fragment representing the top face of the chunk,
  // and write it out.
  planar_polymesh_t* fragment = planar_polymesh_new(chunk->num_columns,
                                                    chunk->num_xy_faces,
                                                    chunk->num_xy_nodes);
  memcpy(fragment->cell_edge_offsets, chunk->column_xy_face_offsets, sizeof(int) * (chunk->num_columns+1));
  planar_polymesh_reserve_connectivity_storage(fragment);
  memcpy(fragment->cell_edges, chunk->column_xy_faces, sizeof(int) * fragment->cell_edge_offsets[fragment->num_cells]);
  memcpy(fragment->edge_cells, chunk->xy_face_columns, 2 * sizeof(int) * chunk->num_xy_faces);
  memcpy(fragment->edge_nodes, chunk->xy_edge_nodes, 2 * sizeof(int) * chunk->num_xy_edges);
  memcpy(fragment->nodes, chunk->xy_nodes, sizeof(point2_t) * chunk->num_xy_nodes);
  char pp_name[FILENAME_MAX+1];
  snprintf(pp_name, FILENAME_MAX, "%s_pp", chunk_grid_name);
  silo_file_write_planar_polymesh(file, pp_name, fragment);
  planar_polymesh_free(fragment);

  // Write an exchanger that holds send and receive maps for the fragment.
  char ex_name[FILENAME_MAX+1];
  snprintf(ex_name, FILENAME_MAX, "%s_ex", chunk_grid_name);
  exchanger_t* fragment_ex = exchanger_new(silo_file_comm(file));
  exchanger_proc_map_t* send_map = exchanger_proc_map_clone(colmesh_xy_data_send_map(mesh, xy_index), NULL, clone_int_array, NULL, int_array_free);
  exchanger_set_sends(fragment_ex, send_map);
  exchanger_proc_map_t* receive_map = exchanger_proc_map_clone(colmesh_xy_data_receive_map(mesh, xy_index), NULL, clone_int_array, NULL, int_array_free);
  exchanger_set_receives(fragment_ex, receive_map);
  silo_file_write_exchanger(file, ex_name, fragment_ex);
  release_ref(fragment_ex);

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
    write_colmesh_chunk_grid(file, chunk_grid_name, mesh, xy, z, mapping);

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
  colmesh_fragment_map_t* fragments = colmesh_fragment_map_new();
  for (int f = 0; f < num_chunk_indices/2; ++f)
  {
    int xy = chunk_indices[2*f];
    char fragment_name[FILENAME_MAX+1];
    snprintf(fragment_name, FILENAME_MAX, "%s_%d_pp", mesh_name, xy);
    planar_polymesh_t* frag_mesh = silo_file_read_planar_polymesh(file, fragment_name);

    char fragment_ex_name[FILENAME_MAX+1];
    snprintf(fragment_ex_name, FILENAME_MAX, "%s_%d_ex", mesh_name, xy);
    exchanger_t* fragment_ex = silo_file_read_exchanger(file, fragment_ex_name, silo_file_comm(file));

    exchanger_proc_map_t* send_map = exchanger_proc_map_new();
    int pos = 0, proc, *indices, num_indices;
    while (exchanger_next_send(fragment_ex, &pos, &proc, &indices, &num_indices))
      for (int i = 0; i < num_indices; ++i)
        exchanger_proc_map_add_index(send_map, proc, indices[i]);

    exchanger_proc_map_t* receive_map = exchanger_proc_map_new();
    pos = 0;
    while (exchanger_next_receive(fragment_ex, &pos, &proc, &indices, &num_indices))
      for (int i = 0; i < num_indices; ++i)
        exchanger_proc_map_add_index(receive_map, proc, indices[i]);
    release_ref(fragment_ex);

    // Add our fragment to the map.
    colmesh_fragment_map_add(fragments, xy, frag_mesh, send_map, receive_map);
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

  // Create the mesh.
#if POLYMEC_HAVE_MPI
  MPI_Comm comm = silo_file_comm(file);
#else
  MPI_Comm comm = MPI_COMM_WORLD;
#endif
  colmesh_t* mesh = create_empty_colmesh_from_fragments(comm, fragments, z1, z2,
                                                        num_xy_chunks, num_z_chunks,
                                                        nz_per_chunk, periodic);

  // Insert our local chunks.
  for (size_t i = 0; i < num_chunk_indices/2; ++i)
  {
    int xy_index = chunk_indices[2*i];
    int z_index  = chunk_indices[2*i+1];
    colmesh_insert_chunk(mesh, xy_index, z_index);
  }
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
  bool exists = (chunk_indices != NULL);
  if (exists)
  {
    // Check for the column data.
    for (size_t i = 0; i < num_chunk_indices/2; ++i)
    {
      int xy = chunk_indices[i/2];
      char columns_name[FILENAME_MAX+1];
      snprintf(columns_name, FILENAME_MAX, "%s_%d_pp", mesh_name, xy);
      exists = silo_file_contains_planar_polymesh(file, columns_name);
      if (!exists) break;
    }
  }
  polymec_free(chunk_indices);
  silo_file_pop_dir(file);
  return exists;
}

static void query_colmesh_vector_comps(colmesh_chunk_data_t* chunk_data,
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
    memset(is_vector_comp, 0, sizeof(bool) * chunk_data->num_components);
    *first_vector_comp = -1;
  }
}

static void copy_out_colmesh_node_component(colmesh_chunk_data_t* chunk_data,
                                            field_metadata_t* md,
                                            int c,
                                            coord_mapping_t* mapping,
                                            real_t* data)
{
  // If we're given a mapping and there are vector fields, we need to treat
  // them specially.
  bool is_vector_comp[chunk_data->num_components];
  int first_vector_comp;
  query_colmesh_vector_comps(chunk_data, md, mapping, is_vector_comp, &first_vector_comp);

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
                                              field_metadata_t* md,
                                              int c,
                                              coord_mapping_t* mapping,
                                              real_t* data)
{
  // If we're given a mapping and there are vector fields, we need to treat
  // them specially.
  bool is_vector_comp[chunk_data->num_components];
  int first_vector_comp;
  query_colmesh_vector_comps(chunk_data, md, mapping, is_vector_comp, &first_vector_comp);

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
                                             field_metadata_t* md,
                                             int c,
                                             coord_mapping_t* mapping,
                                             real_t* data)
{
  // If we're given a mapping and there are vector fields, we need to treat
  // them specially.
  bool is_vector_comp[chunk_data->num_components];
  int first_vector_comp;
  query_colmesh_vector_comps(chunk_data, md, mapping, is_vector_comp, &first_vector_comp);

  // Now copy the data.
  colmesh_chunk_t* chunk = chunk_data->chunk;
  size_t l = chunk->num_xy_edges * (chunk->num_z_cells + 1);
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
                                              field_metadata_t* md,
                                              int c,
                                              coord_mapping_t* mapping,
                                              real_t* data)
{
  // If we're given a mapping and there are vector fields, we need to treat
  // them specially.
  bool is_vector_comp[chunk_data->num_components];
  int first_vector_comp;
  query_colmesh_vector_comps(chunk_data, md, mapping, is_vector_comp, &first_vector_comp);

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
                                             field_metadata_t* md,
                                             int c,
                                             coord_mapping_t* mapping,
                                             real_t* data)
{
  // If we're given a mapping and there are vector fields, we need to treat
  // them specially.
  bool is_vector_comp[chunk_data->num_components];
  int first_vector_comp;
  query_colmesh_vector_comps(chunk_data, md, mapping, is_vector_comp, &first_vector_comp);

  // Now copy the data from the chunk.
  colmesh_chunk_t* chunk = chunk_data->chunk;
  size_t l = chunk->num_xy_faces * chunk->num_z_cells;
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
                                            field_metadata_t* md,
                                            int c,
                                            coord_mapping_t* mapping,
                                            real_t* data)
{
  // If we're given a mapping and there are vector fields, we need to treat
  // them specially.
  bool is_vector_comp[chunk_data->num_components];
  int first_vector_comp;
  query_colmesh_vector_comps(chunk_data, md, mapping, is_vector_comp, &first_vector_comp);

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
      for (int z = 0; z < chunk->num_z_cells+2; ++z, ++l)
        data[l] = a[xy][z][c];
  }
}

// This guy copies out the other centerings in an edge- or face-centered field,
// based on the centering in the given chunk, and sets the ready_to_write
// flag based on whether both centerings (xy and z) will be in the data
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
        memcpy(data, *other_p, sizeof(real_t) * chunk->num_xy_edges * (chunk->num_z_cells+1));
      else
        *ready_to_write = false;
    }
    if (chunk_data->centering != COLMESH_ZEDGE)
    {
      snprintf(scratch_name, FILENAME_MAX, "%s_z", field_component_name);
      real_t** other_p = (real_t**)string_ptr_unordered_map_get(scratch, scratch_name);
      if (other_p != NULL)
      {
        size_t offset = chunk->num_xy_edges * (chunk->num_z_cells + 1);
        memcpy(&data[offset], *other_p, sizeof(real_t) * chunk->num_xy_nodes * chunk->num_z_cells);
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
        memcpy(data, *other_p, sizeof(real_t) * chunk->num_xy_faces * chunk->num_z_cells);
      else
        *ready_to_write = false;
    }
    if (chunk_data->centering != COLMESH_ZFACE)
    {
      snprintf(scratch_name, FILENAME_MAX, "%s_z", field_component_name);
      real_t** other_p = (real_t**)string_ptr_unordered_map_get(scratch, scratch_name);
      if (other_p != NULL)
      {
        size_t offset = chunk->num_xy_faces * chunk->num_z_cells;
        memcpy(&data[offset], *other_p, sizeof(real_t) * chunk->num_columns * (chunk->num_z_cells+1));
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
          for (int z = 0; z < chunk_data->z_size; ++z, ++l)
            this_one[l] = a[xy][z][c];
      }
      else
      {
        snprintf(scratch_name, FILENAME_MAX, "%s_z", field_component_name);
        DECLARE_COLMESH_ZFACE_ARRAY(a, chunk_data);
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

  // If we're ready to write to the file, we can clean up our scratch data.
  if (*ready_to_write)
  {
    snprintf(scratch_name, FILENAME_MAX, "%s_xy", field_component_name);
    string_ptr_unordered_map_delete(scratch, scratch_name);
    snprintf(scratch_name, FILENAME_MAX, "%s_z", field_component_name);
    string_ptr_unordered_map_delete(scratch, scratch_name);
  }
}

static void write_colmesh_chunk_data(silo_file_t* file,
                                     const char** field_component_names,
                                     const char* chunk_grid_name,
                                     colmesh_chunk_data_t* chunk_data,
                                     field_metadata_t* md,
                                     coord_mapping_t* mapping)
{
  // Because we can't really represent a colmesh faithfully in a SILO format,
  // we're really just dumping a bunch of data into an array.
  size_t data_size;
  colmesh_chunk_t* chunk = chunk_data->chunk;
  if ((chunk_data->centering == COLMESH_CELL) ||
      (chunk_data->centering == COLMESH_NODE))
    data_size = chunk_data->xy_size * chunk_data->z_size;
  else if ((chunk_data->centering == COLMESH_XYFACE) ||
           (chunk_data->centering == COLMESH_ZFACE))
  {
    data_size = chunk->num_xy_faces * chunk->num_z_cells +
                chunk->num_columns * (chunk->num_z_cells+1);
  }
  else // edge centerings
  {
    data_size = chunk->num_xy_edges * (chunk->num_z_cells+1) +
                chunk->num_xy_nodes * chunk->num_z_cells;
  }
  real_t* data = polymec_calloc(data_size, sizeof(real_t));

  // Now write each component.
  for (int c = 0; c < chunk_data->num_components; ++c)
  {
    // Copy the data in the component into our array.
    bool ready_to_write = false;
    if (chunk_data->centering == COLMESH_NODE)
    {
      copy_out_colmesh_node_component(chunk_data, md, c, mapping, data);
      ready_to_write = true;
    }
    else if ((chunk_data->centering == COLMESH_XYEDGE) ||
             (chunk_data->centering == COLMESH_ZEDGE))
    {
      copy_out_other_centerings(file, chunk_data, field_component_names[c], c, data, &ready_to_write);
      if (chunk_data->centering == COLMESH_XYEDGE)
        copy_out_colmesh_xyedge_component(chunk_data, md, c, mapping, data);
      else // (data->centering == COLMESH_ZEDGE)
        copy_out_colmesh_zedge_component(chunk_data, md, c, mapping, data);
    }
    else if ((chunk_data->centering == COLMESH_XYFACE) ||
             (chunk_data->centering == COLMESH_ZFACE))
    {
      copy_out_other_centerings(file, chunk_data, field_component_names[c], c, data, &ready_to_write);
      if (chunk_data->centering == COLMESH_XYFACE)
        copy_out_colmesh_xyface_component(chunk_data, md, c, mapping, data);
      else // (data->centering == COLMESH_ZFACE)
        copy_out_colmesh_zface_component(chunk_data, md, c, mapping, data);
    }
    else
    {
      ASSERT(chunk_data->centering == COLMESH_CELL);
      copy_out_colmesh_cell_component(chunk_data, md, c, mapping, data);
      ready_to_write = true;
    }

    // Write the component to the file if it's ready.
    if (ready_to_write)
    {
      char data_name[FILENAME_MAX+1];
      snprintf(data_name, FILENAME_MAX, "%s_%s", chunk_grid_name, field_component_names[c]);
      silo_file_write_real_array(file, data_name, data, data_size);
    }
  }

  // Clean up.
  polymec_free(data);
}

void silo_file_write_colmesh_field(silo_file_t* file,
                                   const char* field_name,
                                   const char* mesh_name,
                                   colmesh_field_t* field,
                                   coord_mapping_t* mapping)
{
  START_FUNCTION_TIMER();
  silo_file_push_domain_dir(file);

  size_t num_components = colmesh_field_num_components(field);
  char* field_names[num_components];

  field_metadata_t* md = colmesh_field_metadata(field);
  char md_name[FILENAME_MAX+1];
  snprintf(md_name, FILENAME_MAX, "%s_%s_md", field_name, mesh_name);
  silo_file_write_field_metadata(file, md_name, md);

  colmesh_chunk_data_t* data;
  int pos = 0, xy, z, l = 0;
  while (colmesh_field_next_chunk(field, &pos, &xy, &z, &data))
  {
    // Write out the chunk data itself.
    for (int c = 0; c < num_components; ++c)
    {
      char field_comp_name[FILENAME_MAX+1];
      snprintf(field_comp_name, FILENAME_MAX, "%s_%d", field_name, c);
      field_names[c] = string_dup(field_comp_name);
    }

    char chunk_grid_name[FILENAME_MAX];
    snprintf(chunk_grid_name, FILENAME_MAX-1, "%s_%d_%d", mesh_name, xy, z);
    write_colmesh_chunk_data(file, (const char**)field_names, chunk_grid_name,
                             data, md, mapping);
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
  size_t l = chunk_data->chunk->num_xy_edges * (chunk_data->chunk->num_z_cells+1);
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
  size_t l = chunk_data->chunk->num_xy_faces * chunk_data->chunk->num_z_cells;
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
    for (int z = 0; z < chunk_data->chunk->num_z_cells+2; ++z, ++l)
      a[xy][z][c] = data[l];
}

static void read_colmesh_chunk_data(silo_file_t* file,
                                    const char** field_component_names,
                                    const char* chunk_grid_name,
                                    size_t num_components,
                                    colmesh_chunk_data_t* chunk_data)
{
  // Fetch each component from the file.
  for (int c = 0; c < (int)num_components; ++c)
  {
    char data_name[FILENAME_MAX+1];
    snprintf(data_name, FILENAME_MAX, "%s_%s", chunk_grid_name, field_component_names[c]);
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
    polymec_free(data);
  }
}

void silo_file_read_colmesh_field(silo_file_t* file,
                                  const char* field_name,
                                  const char* mesh_name,
                                  colmesh_field_t* field)
{
  START_FUNCTION_TIMER();
  silo_file_push_domain_dir(file);

  size_t num_components = colmesh_field_num_components(field);
  char* field_names[num_components];
  field_metadata_t* md = colmesh_field_metadata(field);
  char md_name[FILENAME_MAX+1];
  snprintf(md_name, FILENAME_MAX, "%s_%s_md", field_name, mesh_name);
  silo_file_read_field_metadata(file, md_name, md);

  colmesh_chunk_data_t* data;
  int pos = 0, xy, z;
  while (colmesh_field_next_chunk(field, &pos, &xy, &z, &data))
  {
    for (int c = 0; c < num_components; ++c)
    {
      char field_comp_name[FILENAME_MAX];
      snprintf(field_comp_name, FILENAME_MAX-1, "%s_%d", field_name, c);
      field_names[c] = string_dup(field_comp_name);
    }
    char chunk_grid_name[FILENAME_MAX+1];
    snprintf(chunk_grid_name, FILENAME_MAX, "%s_%d_%d", mesh_name, xy, z);
    read_colmesh_chunk_data(file, (const char**)field_names, chunk_grid_name,
                            data->num_components, data);

    for (int c = 0; c < num_components; ++c)
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
  bool result = false;
  if (silo_file_contains_colmesh(file, mesh_name))  // mesh exists...
  {
    // Look for the field's metadata array.
    char md_name[FILENAME_MAX+1];
    snprintf(md_name, FILENAME_MAX, "%s_%s_md", field_name, mesh_name);
    result = silo_file_contains_field_metadata(file, md_name);
  }
  return result;
}

