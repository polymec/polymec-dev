// Copyright (c) 2012-2018, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/timer.h"
#include "core/array.h"
#include "geometry/blockmesh.h"
#include "geometry/blockmesh_field.h"

#if POLYMEC_HAVE_MPI
#include "core/partitioning.h"
#endif

#if POLYMEC_HAVE_OPENMP
#include <omp.h>
#endif

DEFINE_ARRAY(unimesh_array, unimesh_t*)

struct blockmesh_t
{
  // Parallel stuff.
  MPI_Comm comm;
  int nproc, rank;

  // Blocks.
  unimesh_array_t* blocks;

  // Connections.
  // FIXME

  // This flag is set by blockmesh_finalize() after a mesh has been assembled.
  bool finalized;
};

blockmesh_t* blockmesh_new(MPI_Comm comm)
{
  blockmesh_t* mesh = polymec_malloc(sizeof(blockmesh_t));
  mesh->comm = comm;
  MPI_Comm_size(comm, &mesh->nproc);
  MPI_Comm_rank(comm, &mesh->rank);
  mesh->blocks = unimesh_array_new();
  mesh->finalized = false;

  return mesh;
}

int blockmesh_add_block(blockmesh_t* mesh, unimesh_t* block)
{
  ASSERT(!mesh->finalized);
  ASSERT(block != NULL);
  ASSERT(!unimesh_is_finalized(block));
  int index = (int)mesh->blocks->size;
  unimesh_array_append_with_dtor(mesh->blocks, block, unimesh_free);
  return index;
}

int blockmesh_face_for_nodes(blockmesh_t* mesh, int block_nodes[4])
{
  // Only certain combos of block faces are acceptible, so let's make sure
  // no one's doing anything stupid.
  int face_nodes[6][4] = {{0, 4, 3, 7},  // -x
                          {1, 2, 6, 5},  // +x
                          {0, 1, 5, 4},  // -y
                          {2, 3, 7, 6},  // +y
                          {0, 1, 2, 3},  // -z
                          {4, 5, 6, 7}}; // +z

  int face = -1;
  for (int f = 0; f < 6; ++f)
  {
    bool face_matches = true;
    for (int n = 0; n < 4; ++n)
    {
      bool node_matches = false;
      for (int nn = 0; nn < 4; ++nn)
      {
        if (block_nodes[n] == face_nodes[f][nn])
        {
          node_matches = true;
          break;
        }
      }
      if (!node_matches)
        face_matches = false;
    }
    if (face_matches)
    {
      face = f;
      break;
    }
  }
  return face;
}

void blockmesh_connect_blocks(blockmesh_t* mesh, 
                              int block1_index, 
                              int block1_nodes[4],
                              int block2_index,
                              int block2_nodes[4])
{
}

void blockmesh_finalize(blockmesh_t* mesh)
{
  START_FUNCTION_TIMER();

  ASSERT(!mesh->finalized);
  mesh->finalized = true;

  // FIXME

  STOP_FUNCTION_TIMER();
}

bool blockmesh_is_finalized(blockmesh_t* mesh)
{
  return mesh->finalized;
}

void blockmesh_free(blockmesh_t* mesh)
{
  unimesh_array_free(mesh->blocks);
  polymec_free(mesh);
}

MPI_Comm blockmesh_comm(blockmesh_t* mesh)
{
  return mesh->comm;
}

size_t blockmesh_num_blocks(blockmesh_t* mesh)
{
  return mesh->blocks->size;
}

unimesh_t* blockmesh_block(blockmesh_t* mesh, int index)
{
  ASSERT(index >= 0);
  ASSERT((size_t)index < mesh->blocks->size);
  return mesh->blocks->data[index];
}

#if POLYMEC_HAVE_MPI
static adj_graph_t* graph_from_blocks(blockmesh_t* mesh)
{
  START_FUNCTION_TIMER();
  // Create a graph whose vertices are the mesh's blocks. NOTE
  // that we associate this graph with the MPI_COMM_SELF communicator 
  // because it's a global graph.
  adj_graph_t* g = adj_graph_new(MPI_COMM_SELF, mesh->blocks->size);

#if 0
  // Allocate space in the graph for the edges (patch boundaries).
  for (int i = 0; i < mesh->npx; ++i)
  {
    int num_x_edges = (i == 0) ? mesh->periodic_in_x ? 2 : 1
                               : (i == mesh->npx-1) ? mesh->periodic_in_x ? 2 : 1 
                                                    : 2;
    for (int j = 0; j < mesh->npy; ++j)
    {
      int num_y_edges = (j == 0) ? mesh->periodic_in_y ? 2 : 1
                                 : (j == mesh->npy-1) ? mesh->periodic_in_z ? 2 : 1 
                                                      : 2;
      for (int k = 0; k < mesh->npz; ++k)
      {
        int num_z_edges = (k == 0) ? mesh->periodic_in_z ? 2 : 1
                                   : (k == mesh->npz-1) ? mesh->periodic_in_z ? 2 : 1 
                                                        : 2;
        int num_edges = num_x_edges + num_y_edges + num_z_edges;
        int p_index = patch_index(mesh, i, j, k);
        adj_graph_set_num_edges(g, p_index, num_edges);
      }
    }
  }

  // Now fill in the edges.
  for (int i = 0; i < mesh->npx; ++i)
  {
    for (int j = 0; j < mesh->npy; ++j)
    {
      for (int k = 0; k < mesh->npz; ++k)
      {
        int p_index = patch_index(mesh, i, j, k);
        int* edges = adj_graph_edges(g, p_index);
        int offset = 0;

        if ((i == 0) && mesh->periodic_in_x)
          edges[offset++] = patch_index(mesh, mesh->npx-1, j, k);
        else if (i > 0)
          edges[offset++] = patch_index(mesh, i-1, j, k);
        if ((i == mesh->npx-1) && mesh->periodic_in_x)
          edges[offset++] = patch_index(mesh, 0, j, k);
        else if (i < mesh->npx-1)
          edges[offset++] = patch_index(mesh, i+1, j, k);

        if ((j == 0) && mesh->periodic_in_y)
          edges[offset++] = patch_index(mesh, i, mesh->npy-1, k);
        else if (j > 0)
          edges[offset++] = patch_index(mesh, i, j-1, k);
        if ((j == mesh->npy-1) && mesh->periodic_in_y)
          edges[offset++] = patch_index(mesh, i, 0, k);
        else if (j < mesh->npy-1)
          edges[offset++] = patch_index(mesh, i, j+1, k);

        if ((k == 0) && mesh->periodic_in_z)
          edges[offset++] = patch_index(mesh, i, j, mesh->npz-1);
        else if (k > 0)
          edges[offset++] = patch_index(mesh, i, j, k-1);
        if ((k == mesh->npz-1) && mesh->periodic_in_z)
          edges[offset++] = patch_index(mesh, i, j, 0);
        else if (k < mesh->npz-1)
          edges[offset++] = patch_index(mesh, i, j, k+1);
      }
    }
  }
#endif

  STOP_FUNCTION_TIMER();
  return g;
}

static int64_t* source_vector(blockmesh_t* mesh)
{
  START_FUNCTION_TIMER();
  STOP_FUNCTION_TIMER();
  return NULL;
}

static void redistribute_blockmesh(blockmesh_t** mesh, 
                                   int64_t* partition)
{
  START_FUNCTION_TIMER();
  STOP_FUNCTION_TIMER();
}

// Redistributes the given block mesh using the given partition vector, but 
// does not finalize the mesh.
static void redistribute_blockmesh_field(blockmesh_field_t** field, 
                                         int64_t* partition,
                                         int64_t* sources,
                                         blockmesh_t* new_mesh)
{
  START_FUNCTION_TIMER();
  STOP_FUNCTION_TIMER();
}
#endif

void repartition_blockmesh(blockmesh_t** mesh, 
                           int* weights,
                           real_t imbalance_tol,
                           blockmesh_field_t** fields,
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
  blockmesh_t* old_mesh = *mesh;
  if (old_mesh->nproc == 1) 
  {
    STOP_FUNCTION_TIMER();
    return;
  }

  // Generate a global adjacency graph for the mesh.
  adj_graph_t* graph = graph_from_blocks(old_mesh);

  // Map the graph to the different domains, producing a partition vector.
  // We need the partition vector on all processes, so we scatter it 
  // from rank 0.
  log_debug("repartition_blockmesh: Repartitioning mesh on %d subdomains.", old_mesh->nproc);
  int64_t* partition = partition_graph(graph, old_mesh->comm, weights, imbalance_tol, true);

  // Redistribute the mesh. 
  log_debug("repartition_blockmesh: Redistributing mesh.");
  redistribute_blockmesh(mesh, partition);
  blockmesh_finalize(*mesh);

  // Build a sources vector whose ith component is the rank that used to own 
  // the ith patch.
  int64_t* sources = source_vector(old_mesh);

  // Redistribute the fields.
  if (num_fields > 0)
    log_debug("repartition_blockmesh: Redistributing %d fields.", (int)num_fields);
  for (size_t f = 0; f < num_fields; ++f)
  {
    blockmesh_field_t* old_field = fields[f];
    redistribute_blockmesh_field(&(fields[f]), partition, sources, *mesh);
    blockmesh_field_free(old_field);
  }

  // Clean up.
  blockmesh_free(old_mesh);
  adj_graph_free(graph);
  polymec_free(sources);
  polymec_free(partition);

  STOP_FUNCTION_TIMER();
#endif
}
