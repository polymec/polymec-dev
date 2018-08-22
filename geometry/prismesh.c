// Copyright (c) 2012-2018, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/timer.h"
#include "core/array.h"
#include "geometry/prismesh.h"
#include "geometry/prismesh_field.h"
#include "geometry/polymesh.h"

#if POLYMEC_HAVE_MPI
#include "core/partitioning.h"
#endif

DEFINE_ARRAY(chunk_array, prismesh_chunk_t*)

static prismesh_chunk_t* chunk_from_planar_polymesh(planar_polymesh_t* mesh,
                                                    real_t z1, real_t z2,
                                                    size_t nz)
{
  ASSERT(z1 < z2);
  prismesh_chunk_t* chunk = polymec_malloc(sizeof(prismesh_chunk_t));
  chunk->num_columns = (size_t)mesh->num_cells;
  chunk->num_z_cells = nz;
  chunk->z1 = z1;
  chunk->z2 = z2;

  // cell -> xy face connectivity.
  chunk->column_xy_face_offsets = polymec_malloc(sizeof(int) * (chunk->num_columns+1));
  memcpy(chunk->column_xy_face_offsets, mesh->cell_edge_offsets, sizeof(int) * (chunk->num_columns+1));

  chunk->column_xy_faces = polymec_malloc(sizeof(int) * chunk->column_xy_face_offsets[chunk->num_columns]);
  memcpy(chunk->column_xy_faces, mesh->cell_edges, sizeof(int) * chunk->column_xy_face_offsets[chunk->num_columns]);

  // xy face -> cell connectivity.
  chunk->num_xy_faces = (size_t)mesh->num_edges;
  chunk->xy_face_columns = polymec_malloc(sizeof(int) * 2 * chunk->num_xy_faces);
  memcpy(chunk->xy_face_columns, mesh->edge_cells, sizeof(int) * 2 * chunk->num_xy_faces);

  // Node xy coordinates.
  chunk->num_xy_edges = (size_t)mesh->num_edges;
  chunk->num_xy_nodes = (size_t)mesh->num_nodes;
  chunk->xy_nodes = polymec_malloc(sizeof(point2_t) * chunk->num_xy_nodes);
  memcpy(chunk->xy_nodes, mesh->nodes, sizeof(point2_t) * chunk->num_xy_nodes);

  return chunk;
}

static void free_chunk(prismesh_chunk_t* chunk)
{
  polymec_free(chunk->xy_nodes);
  polymec_free(chunk->xy_face_columns);
  polymec_free(chunk->column_xy_faces);
  polymec_free(chunk->column_xy_face_offsets);
  polymec_free(chunk);
}

struct prismesh_t 
{
  MPI_Comm comm;
  int nproc, rank;

  planar_polymesh_t* columns;
  chunk_array_t* chunks;
  int_array_t* xy_indices;
  int_array_t* z_indices;
  size_t num_xy_chunks, num_z_chunks, nz_per_chunk;
  real_t z1, z2;

  // This flag is set by prismesh_finalize() after a mesh has been assembled.
  bool finalized;
};

prismesh_t* create_empty_prismesh(MPI_Comm comm, 
                                  planar_polymesh_t* columns,
                                  real_t z1, real_t z2,
                                  size_t num_xy_chunks, size_t num_z_chunks,
                                  size_t nz_per_chunk)
{
  ASSERT(columns != NULL);
  ASSERT(z1 < z2);
  ASSERT(num_xy_chunks > 0);
  ASSERT(num_z_chunks > 0);
  ASSERT(nz_per_chunk > 0);

  prismesh_t* mesh = polymec_malloc(sizeof(prismesh_t));
  mesh->comm = comm;
  mesh->columns = columns; 
  mesh->chunks = chunk_array_new();
  mesh->xy_indices = int_array_new();
  mesh->z_indices = int_array_new();
  mesh->num_xy_chunks = num_xy_chunks;
  mesh->num_z_chunks = num_z_chunks;
  mesh->nz_per_chunk = nz_per_chunk;
  mesh->z1 = z1;
  mesh->z2 = z2;
  MPI_Comm_size(comm, &mesh->nproc);
  MPI_Comm_rank(comm, &mesh->rank);
  mesh->finalized = false;

  return mesh;
}

void prismesh_insert_chunk(prismesh_t* mesh, int xy_index, int z_index)
{
  ASSERT(!mesh->finalized);
  ASSERT(xy_index >= 0);
  ASSERT(xy_index < mesh->num_xy_chunks);
  ASSERT(z_index >= 0);
  ASSERT(z_index < mesh->num_z_chunks);

  size_t nz = mesh->nz_per_chunk * mesh->num_z_chunks;
  real_t dz = (mesh->z2 - mesh->z1) / nz;
  real_t z1 = mesh->z1 + z_index * dz;
  real_t z2 = mesh->z1 + (z_index+1) * dz;
  prismesh_chunk_t* chunk = chunk_from_planar_polymesh(mesh->columns, 
                                                       z1, z2, 
                                                       mesh->nz_per_chunk);
  chunk_array_append_with_dtor(mesh->chunks, chunk, free_chunk);
  int_array_append(mesh->xy_indices, xy_index);
  int_array_append(mesh->z_indices, z_index);
}

void prismesh_finalize(prismesh_t* mesh)
{
  ASSERT(!mesh->finalized);
  planar_polymesh_free(mesh->columns);
  mesh->columns = NULL;
  mesh->finalized = true;
}

prismesh_t* prismesh_new(MPI_Comm comm,
                         planar_polymesh_t* columns,
                         real_t z1, real_t z2,
                         size_t nz)
{
  prismesh_t* mesh = create_empty_prismesh(comm, columns, z1, z2, 10, 5, nz);

#if 0
  if (mesh->nproc > 1)
  {
    // Figure out an optimal partitioning "geometry". How many "axial" 
    // processes do we want, vs how many process in the xy plane?
    int num_xy_procs = (int)(floor(pow(mesh->nproc, 2.0/3.0)));
    int num_z_procs = mesh->nproc / num_xy_procs;

    if (num_xy_procs > 1)
    {
      // Create an "xy communicator" that exists across the xy plane.
      MPI_Comm xy_comm;
      int color = mesh->rank / num_z_procs;
      int key = mesh->rank % num_z_procs;
      MPI_Comm_split(comm, color, key, &xy_comm);

      // Construct a graph for the planar polymesh.
      adj_graph_t* graph = graph_from_planar_polymesh_cells(columns);

      // Cut it up in the xy plane and get the partition vector.
      int64_t* P = partition_graph(graph, xy_comm, NULL, 0.05, true);

      // Allocate chunks to this process.
      // FIXME: We could store chunk xy connectivity more efficiently by 
      // FIXME: aliasing all chunks' data to the xy data in the first 
      // FIXME: chunk created in each column.
      real_t dz = (z2 - z1) / num_z_procs;
      for (int ij = 0; ij < num_xy_procs; ++ij)
      {
        for (int k = 0; k < num_z_procs; ++k)
        {
          real_t z1_k = z1 + k*dz;
          real_t z2_k = z1 + (k+1)*dz;
          size_t nz = num_vertical_cells / num_z_procs;
          if ((k*num_vertical_cells/num_z_procs + nz) > num_vertical_cells)
            nz = num_vertical_cells - k*num_vertical_cells/num_z_procs;
          prismesh_chunk_t* chunk = chunk_from_planar_polymesh(columns, nz, z1_k, z2_k);
          chunk_array_append_with_dtor(mesh->chunks, chunk, free_chunk);
        }
      }

      // Clean up.
      MPI_Comm_free(&xy_comm);
      polymec_free(P);
    }
    else
    {
      // Allocate chunks to this process along the z axis.
      // FIXME: We could store chunk xy connectivity more efficiently by 
      // FIXME: aliasing all chunks' data to the xy data in the first 
      // FIXME: chunk created in each column.
      real_t dz = (z2 - z1) / num_z_procs;
      for (int k = 0; k < num_z_procs; ++k)
      {
        real_t z1_k = z1 + k*dz;
        real_t z2_k = z1 + (k+1)*dz;
        size_t nz = num_vertical_cells / num_z_procs;
        if ((k*num_vertical_cells/num_z_procs + nz) > num_vertical_cells)
          nz = num_vertical_cells - k*num_vertical_cells/num_z_procs;
        prismesh_chunk_t* chunk = chunk_from_planar_polymesh(columns, nz, z1_k, z2_k);
        chunk_array_append_with_dtor(mesh->chunks, chunk, free_chunk);
      }
    }
  }
  else
  {
    // Allocate a single chunk to this process with columns corresponding 
    // to the cells of the planar polymesh. 
    prismesh_chunk_t* chunk = chunk_from_planar_polymesh(columns, 
                                                         num_vertical_cells, 
                                                         z1, z2);
    chunk_array_append_with_dtor(mesh->chunks, chunk, free_chunk);
  }
#endif

  prismesh_finalize(mesh);
  return mesh;
}

void prismesh_free(prismesh_t* mesh)
{
  int_array_free(mesh->z_indices);
  int_array_free(mesh->xy_indices);
  chunk_array_free(mesh->chunks);
  if (mesh->columns != NULL)
    planar_polymesh_free(mesh->columns);
  polymec_free(mesh);
}

MPI_Comm prismesh_comm(prismesh_t* mesh)
{
  return mesh->comm;
}

size_t prismesh_num_chunks(prismesh_t* mesh)
{
  return mesh->chunks->size;
}

size_t prismesh_num_xy_chunks(prismesh_t* mesh)
{
  return mesh->num_xy_chunks;
}

size_t prismesh_num_z_chunks(prismesh_t* mesh)
{
  return mesh->num_z_chunks;
}

real_t prismesh_z1(prismesh_t* mesh)
{
  return mesh->z1;
}

real_t prismesh_z2(prismesh_t* mesh)
{
  return mesh->z2;
}

polygon_t* prismesh_chunk_polygon(prismesh_chunk_t* chunk, int column)
{
  ASSERT(column < (int)chunk->num_columns);
  int z_face = column;
  int num_nodes = prismesh_chunk_z_face_num_nodes(chunk, z_face);
  int nodes[num_nodes];
  prismesh_chunk_z_face_get_nodes(chunk, z_face, nodes);
  point2_t vertices[num_nodes];
  for (int n = 0; n < num_nodes; ++n)
    vertices[n] = chunk->xy_nodes[nodes[n]];
  return polygon_new(vertices, num_nodes);
}

bool prismesh_next_chunk(prismesh_t* mesh, int* pos, prismesh_chunk_t** chunk)
{
  if (*pos >= (int)mesh->chunks->size)
    return false;
  else
  {
    *chunk = mesh->chunks->data[*pos];
    ++(*pos);
    return true;
  }
}

static void redistribute_prismesh(prismesh_t** mesh, 
                                  MPI_Comm super_comm,
                                  int64_t* partition)
{
  START_FUNCTION_TIMER();
  // FIXME
  STOP_FUNCTION_TIMER();
}

static int64_t* source_vector(prismesh_t* mesh)
{
#if 0
  // Catalog all the chunks on this process.
  int_array_t* my_chunks = int_array_new();
  for (int i = 0; i < mesh->npx; ++i)
  {
    for (int j = 0; j < mesh->npy; ++j)
    {
      for (int k = 0; k < mesh->npz; ++k)
      {
        if (unimesh_has_patch(mesh, i, j, k))
          int_array_append(my_patches, patch_index(mesh, i, j, k));
      }
    }
  }

  // Gather the numbers of patches owned by each process.
  int num_my_patches = (int)my_patches->size;
  int num_patches_for_proc[mesh->nproc];
  MPI_Allgather(&num_my_patches, 1, MPI_INT, 
                num_patches_for_proc, 1, MPI_INT, mesh->comm);

  // Arrange for the storage of the patch indices for the patches stored 
  // on each process.
  int proc_offsets[mesh->nproc+1];
  proc_offsets[0] = 0;
  for (int p = 0; p < mesh->nproc; ++p)
    proc_offsets[p+1] = proc_offsets[p] + num_patches_for_proc[p];

  // Gæther the indices of the patches owned by all processes into a huge list.
  int num_all_patches = mesh->npx * mesh->npy * mesh->npz;
  ASSERT(num_all_patches == proc_offsets[mesh->nproc]);
  int* all_patches = polymec_malloc(sizeof(int) * num_all_patches);
  MPI_Allgatherv(my_patches->data, num_my_patches, MPI_INT, 
                 all_patches, num_patches_for_proc, proc_offsets,
                 MPI_INT, mesh->comm);

  // Clean up a bit.
  int_array_free(my_patches);

  // Convert the huge list into a source vector.
  int64_t* sources = polymec_malloc(sizeof(int64_t) * num_all_patches);
  for (int p = 0; p < mesh->nproc; ++p)
  {
    for (int offset = proc_offsets[p]; offset < proc_offsets[p+1]; ++offset)
      sources[all_patches[offset]] = (int64_t)p;
  }

  polymec_free(all_patches);
  return sources;
#endif
  return NULL;
}

static void redistribute_prismesh_field(prismesh_field_t** field, 
                                        int64_t* partition,
                                        int64_t* sources,
                                        prismesh_t* new_mesh)
{
  START_FUNCTION_TIMER();

  // Create a new field from the old one.
  prismesh_field_t* old_field = *field;
  prismesh_field_t* new_field = prismesh_field_new(new_mesh,
                                                   prismesh_field_centering(old_field),
                                                   prismesh_field_num_components(old_field));
#if 0

  // Copy all local chunks from one field to the other.
  unimesh_patch_t* patch;
  int pos = 0, i, j, k;
  while (unimesh_field_next_patch(new_field, &pos, &i, &j, &k, &patch, NULL))
  {
    unimesh_patch_t* old_patch = unimesh_field_patch(old_field, i, j, k);
    if (old_patch != NULL)
      unimesh_patch_copy(old_patch, patch);
  }

  // Post receives for each patch in the new field.
  int num_new_local_patches = unimesh_field_num_patches(new_field);
  MPI_Request recv_requests[num_new_local_patches];
  pos = 0;
  int num_recv_reqs = 0;
  while (unimesh_field_next_patch(new_field, &pos, &i, &j, &k, &patch, NULL))
  {
    int p = patch_index(new_mesh, i, j, k);
    if (partition[p] == new_mesh->rank)
    {
      size_t data_size = unimesh_patch_data_size(patch->centering, 
                                                 patch->nx, patch->ny, patch->nz,
                                                 patch->nc) / sizeof(real_t);
      int err = MPI_Irecv(patch->data, (int)data_size, MPI_REAL_T, (int)sources[p],
                          0, new_mesh->comm, &(recv_requests[num_recv_reqs]));
      if (err != MPI_SUCCESS)
        polymec_error("Error receiving field data from rank %d", (int)sources[p]);
      ++num_recv_reqs;
    }
  }
  ASSERT(num_recv_reqs <= num_new_local_patches);

  // Post sends.
  int num_old_local_patches = unimesh_field_num_patches(old_field);
  MPI_Request send_requests[num_old_local_patches];
  pos = 0;
  int num_send_reqs = 0;
  while (unimesh_field_next_patch(old_field, &pos, &i, &j, &k, &patch, NULL))
  {
    int p = patch_index(new_mesh, i, j, k);
    if (sources[p] == new_mesh->rank)
    {
      size_t data_size = unimesh_patch_data_size(patch->centering, 
                                                 patch->nx, patch->ny, patch->nz,
                                                 patch->nc) / sizeof(real_t);
      int err = MPI_Isend(patch->data, (int)data_size, MPI_REAL_T, (int)partition[p],
                          0, new_mesh->comm, &(send_requests[num_send_reqs]));
      if (err != MPI_SUCCESS)
        polymec_error("Error sending field data to rank %d", (int)partition[p]);
      ++num_send_reqs;
    }
  }
  ASSERT(num_send_reqs <= num_old_local_patches);

  // Wait for everything to finish.
  MPI_Waitall(num_send_reqs, send_requests, MPI_STATUSES_IGNORE);
  MPI_Waitall(num_recv_reqs, recv_requests, MPI_STATUSES_IGNORE);
#endif

  // Replace the old field with the new one.
  *field = new_field;
  STOP_FUNCTION_TIMER();
}

void repartition_prismesh(prismesh_t** mesh, 
                          int* weights,
                          real_t imbalance_tol,
                          prismesh_field_t** fields,
                          size_t num_fields)
{
  ASSERT((weights == NULL) || (imbalance_tol > 0.0));
  ASSERT((weights == NULL) || (imbalance_tol <= 1.0));
  ASSERT((fields != NULL) || (num_fields == 0));
#if POLYMEC_HAVE_MPI
  START_FUNCTION_TIMER();

  // On a single process, repartitioning has no meaning.
  prismesh_t* old_mesh = *mesh;
  if (old_mesh->nproc == 1) 
  {
    STOP_FUNCTION_TIMER();
    return;
  }

  // Generate a distributed adjacency graph for the polygonal columns.
//  planar_polymesh_t* columns = planar_polymesh_from_columns(old_mesh);
//  adj_graph_t* graph = graph_from_planar_polymesh_cells(columns);

  // Map the graph to the different domains, producing a partition vector.
  // We need the partition vector on all processes in the communicator, so we 
  // scatter it from rank 0.
  log_debug("repartition_prismesh: Repartitioning mesh on %d subdomains.", old_mesh->nproc);
  int64_t* poly_partition = NULL; //partition_graph_n_ways(graph, num_poly_pieces, weights, imbalance_tol);

  // Translate our polygonal partition vector into the real one (which includes
  // axial decomposition).
  int64_t* partition = NULL; // FIXME
  polymec_free(poly_partition);

  // Redistribute the mesh. 
  log_debug("repartition_peximesh: Redistributing mesh.");
  redistribute_prismesh(mesh, old_mesh->comm, partition);

  // Build a sources vector whose ith component is the rank that used to own 
  // the ith patch.
  int64_t* sources = source_vector(old_mesh);

  // Redistribute the fields.
  if (num_fields > 0)
    log_debug("repartition_unimesh: Redistributing %d fields.", (int)num_fields);
  for (size_t f = 0; f < num_fields; ++f)
  {
    prismesh_field_t* old_field = fields[f];
    redistribute_prismesh_field(&(fields[f]), partition, sources, *mesh);
    prismesh_field_free(old_field);
  }

  // Clean up.
  prismesh_free(old_mesh);
//  adj_graph_free(graph);
  polymec_free(sources);
  polymec_free(partition);

  STOP_FUNCTION_TIMER();
#endif
}

polymesh_t* prismesh_as_polymesh(prismesh_t* mesh)
{
  return NULL; // FIXME
}

