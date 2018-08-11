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

#if 0
// Destroys a layer managed by a mesh.
static void free_layer(prismesh_layer_t* layer)
{
  prism_column_array_free(layer->columns);
  polymec_free(layer);
}
#endif

DEFINE_ARRAY(layer_array, prismesh_layer_t*)

struct prismesh_t 
{
  MPI_Comm comm;
  int nproc, rank;

  layer_array_t* layers;
  size_t num_columns, num_vertical_cells;
  real_t z1, z2;
};

prismesh_t* prismesh_new(planar_polymesh_t* columns,
                         size_t num_vertical_cells,
                         real_t z1, real_t z2)
{
  ASSERT(columns != NULL);
  ASSERT(num_vertical_cells > 0);
  ASSERT(z1 < z2);

  prismesh_t* mesh = polymec_malloc(sizeof(prismesh_t));
  mesh->comm = columns->comm;
  mesh->layers = layer_array_new();
  mesh->num_columns = (size_t)columns->num_cells;
  mesh->num_vertical_cells = num_vertical_cells;
  mesh->z1 = z1;
  mesh->z2 = z2;
  return mesh;
}

void prismesh_free(prismesh_t* mesh)
{
  layer_array_free(mesh->layers);
  polymec_free(mesh);
}

MPI_Comm prismesh_comm(prismesh_t* mesh)
{
  return mesh->comm;
}

size_t prismesh_num_layers(prismesh_t* mesh)
{
  return mesh->layers->size;
}

size_t prismesh_num_columns(prismesh_t* mesh)
{
  return mesh->num_columns;
}

size_t prismesh_num_vertical_cells(prismesh_t* mesh)
{
  return mesh->num_vertical_cells;
}

size_t prismesh_num_cells(prismesh_t* mesh)
{
  return mesh->num_columns * mesh->num_vertical_cells;
}

real_t primesh_z1(prismesh_t* mesh)
{
  return mesh->z1;
}

real_t primesh_z2(prismesh_t* mesh)
{
  return mesh->z2;
}

polygon_t* prismesh_polygon(prismesh_t* mesh, size_t column)
{
  ASSERT(column < mesh->num_columns);
  return NULL; // FIXME
}

bool prismesh_next_layer(prismesh_t* mesh, int* pos, prismesh_layer_t** layer)
{
  if (*pos >= (int)mesh->layers->size)
    return false;
  else
  {
    *layer = mesh->layers->data[*pos];
    ++(*pos);
    return true;
  }
}

static void redistribute_prismesh(prismesh_t** mesh, 
                                  int64_t* partition)
{
  START_FUNCTION_TIMER();
  // FIXME
  STOP_FUNCTION_TIMER();
}

static int64_t* source_vector(prismesh_t* mesh)
{
#if 0
  // Catalog all the layers on this process.
  int_array_t* my_layers = int_array_new();
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

  // Copy all local layers from one field to the other.
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
  ASSERT(imbalance_tol > 0.0);
  ASSERT(imbalance_tol <= 1.0);
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

  // Generate a global adjacency graph for the mesh.
  adj_graph_t* graph = NULL;//graph_from_polygons(old_mesh);

  // Figure out how many ways we want to partition the polygon graph by 
  // chopping the z axis into segments of appropriate length. The number of 
  // segments in the z direction is the factor by which we reduce the 
  // number of ways we partition the polygon graph.
  int num_poly_pieces = 1; // FIXME

  // Map the graph to the different domains, producing a partition vector.
  // We need the partition vector on all processes in the communicator, so we 
  // scatter it from rank 0.
  log_debug("repartition_prismesh: Repartitioning mesh on %d subdomains.", old_mesh->nproc);
  int64_t* poly_partition = partition_graph_n_ways(graph, num_poly_pieces, weights, imbalance_tol);

  // Translate our polygonal partition vector into the real one (which includes
  // axial decomposition).
  int64_t* partition = NULL; // FIXME
  polymec_free(poly_partition);

  // Redistribute the mesh. 
  log_debug("repartition_peximesh: Redistributing mesh.");
  redistribute_prismesh(mesh, partition);

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
  adj_graph_free(graph);
  polymec_free(sources);
  polymec_free(partition);

  STOP_FUNCTION_TIMER();
#endif
}

polymesh_t* prismesh_as_polymesh(prismesh_t* mesh)
{
  return NULL; // FIXME
}

