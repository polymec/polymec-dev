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

#include "geometry/repartition.h"

#if POLYMEC_HAVE_MPI
#include "ptscotch.h"

// This helper sets up the exchanger ex so that it can migrate data from the 
// current process according to the global partition vector.
static exchanger_t* create_migrator(SCOTCH_Dgraph* dist_graph,
                                    SCOTCH_Num* global_partition)
{
  MPI_Comm comm;
  SCOTCH_dgraphData(dist_graph, NULL, NULL, NULL, NULL, NULL, NULL,
                    NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
                    &comm);
  exchanger_t* migrator = exchanger_new(comm);
  return migrator;
}

// This helper creates a migration exchanger that describes a repartitioning 
// given a local graph.
static exchanger_t* repartition_graph(adj_graph_t* local_graph, 
                                      int num_ghost_vertices,
                                      int* weights,
                                      real_t imbalance_tol,
                                      exchanger_t* local_graph_ex)
{
  int nprocs, rank;
  MPI_Comm comm = adj_graph_comm(local_graph);
  MPI_Comm_size(comm, &nprocs);
  MPI_Comm_rank(comm, &rank);
  int num_vertices = adj_graph_num_vertices(local_graph);

for (int p = 0; p < nprocs; ++p)
{
if (rank == p)
  adj_graph_fprintf(local_graph, stdout);
MPI_Barrier(comm);
}

  // Extract the adjacency information.
  SCOTCH_Num* xadj = malloc(sizeof(SCOTCH_Num) * (num_vertices+1));
  int* edge_offsets = adj_graph_edge_offsets(local_graph);
  for (int i = 0; i <= num_vertices; ++i)
    xadj[i] = (SCOTCH_Num)edge_offsets[i];
  SCOTCH_Num num_arcs = xadj[num_vertices];
  SCOTCH_Num* adj = malloc(sizeof(SCOTCH_Num) * num_arcs);
  int* edges = adj_graph_adjacency(local_graph);
  for (int i = 0; i < xadj[num_vertices]; ++i)
    adj[i] = (SCOTCH_Num)edges[i];

  // Replace the ghost entries in adj with global indices.
  index_t* vtx_dist = adj_graph_vertex_dist(local_graph);
  {
    index_t* global_indices = polymec_malloc(sizeof(SCOTCH_Num) * (num_vertices + num_ghost_vertices));
    for (int i = 0; i < num_vertices; ++i)
      global_indices[i] = (index_t)(vtx_dist[rank] + i);
exchanger_fprintf(local_graph_ex, stdout);
exchanger_enable_deadlock_detection(local_graph_ex, 1.0, 0, stdout);
    exchanger_exchange(local_graph_ex, global_indices, 1, 0, MPI_UINT64_T);
for (int i = 0; i < num_vertices + num_ghost_vertices; ++i)
  printf("%d: GCI[%d] = %d\n", rank, i, global_indices[i]);
    for (int i = 0; i < num_arcs; ++i)
      adj[i] = global_indices[adj[i]];
    polymec_free(global_indices);
  }

  // Build a distributed graph.
  SCOTCH_Num* vtx_weights = NULL;
  if (weights != NULL)
  {
    vtx_weights = polymec_malloc(sizeof(SCOTCH_Num) * num_vertices);
    for (int i = 0; i < num_vertices; ++i)
      vtx_weights[i] = (SCOTCH_Num)weights[i];
  }
  SCOTCH_Dgraph dist_graph;
  SCOTCH_dgraphInit(&dist_graph, comm);
  SCOTCH_dgraphBuild(&dist_graph, 0, (SCOTCH_Num)num_vertices, (SCOTCH_Num)num_vertices,
                     xadj, NULL, vtx_weights, NULL, num_arcs, num_arcs,
                     adj, NULL, NULL);

  // Generate the local partition vector by mapping the distributed graph.
  SCOTCH_Num* local_partition = polymec_malloc(sizeof(SCOTCH_Num) * num_vertices);
  {
    SCOTCH_Arch arch;
    SCOTCH_archInit(&arch);
    if (vtx_weights == NULL)
      SCOTCH_archCmplt(&arch, num_vertices);
    else
      SCOTCH_archCmpltw(&arch, num_vertices, vtx_weights);
    SCOTCH_Strat strategy;
    SCOTCH_stratInit(&strategy);
    SCOTCH_Num strat_flags = SCOTCH_STRATDEFAULT;
    SCOTCH_stratDgraphMapBuild(&strategy, strat_flags, nprocs, nprocs, (double)imbalance_tol);
    int result = SCOTCH_dgraphMap(&dist_graph, &arch, &strategy, local_partition);
    if (result != 0)
      polymec_error("Repartitioning failed.");
    SCOTCH_stratExit(&strategy);
    SCOTCH_archExit(&arch);
  }
  if (vtx_weights != NULL)
    polymec_free(vtx_weights);

//  // Produce a redistributed graph using the partition vector.
//  SCOTCH_Dgraph redist_graph;
//  SCOTCH_dgraphInit(&redist_graph, comm);
//  result = SCOTCH_dgraphRedist(&dist_graph, partition, NULL, 0, 0, &redist_graph); 

printf("p = [");
for (int i = 0; i < num_vertices; ++i)
printf("%d ", local_partition[i]);
printf("]\n");

  // Assemble a global partition vector.
  index_t global_num_points = vtx_dist[nprocs];
  SCOTCH_Num* global_partition = polymec_malloc(sizeof(SCOTCH_Num) * global_num_points);
  int recv_counts[nprocs], offsets[nprocs];
  for (int p = 0; p < nprocs; ++p)
  {
    recv_counts[p] = (int)(vtx_dist[p+1] - vtx_dist[p]);
    offsets[p] = (int)vtx_dist[p];
  }
  MPI_Allgatherv(local_partition, num_vertices, MPI_UINT64_T, global_partition, 
                 recv_counts, offsets, MPI_UINT64_T, comm);
  polymec_free(local_partition);

  // Clean up.
  SCOTCH_dgraphExit(&dist_graph);

printf("P = [");
for (int i = 0; i < global_num_points; ++i)
printf("%d ", global_partition[i]);
printf("]\n");

  // Set up the migrator to migrate the data.
  exchanger_t* migrator = create_migrator(&dist_graph, global_partition);
  polymec_free(global_partition);

  return migrator;
}

#endif

exchanger_t* repartition_point_cloud(point_cloud_t* cloud, int* weights, real_t imbalance_tol)
{
  ASSERT(imbalance_tol > 0.0);
  ASSERT(imbalance_tol <= 1.0);
#if POLYMEC_HAVE_MPI
  _Static_assert(sizeof(SCOTCH_Num) == sizeof(index_t), "SCOTCH_Num must be 64-bit.");

  int nprocs, rank;
  MPI_Comm_size(cloud->comm, &nprocs);
  MPI_Comm_rank(cloud->comm, &rank);

  // On a single process, repartitioning has no meaning.
  if (nprocs == 1)
    return exchanger_new(cloud->comm);

  // Generate a local adjacency graph for the point cloud.
  adj_graph_t* local_graph = graph_from_point_cloud(cloud);

  // Get the exchanger for the point cloud.
  exchanger_t* cloud_ex = point_cloud_exchanger(cloud);

  // Map the graph to the different domains, producing a migrator.
  exchanger_t* migrator = repartition_graph(local_graph, cloud->num_ghost_points, weights, imbalance_tol, cloud_ex);
  adj_graph_free(local_graph);

  // Migrate the data.
  point_cloud_migrate(cloud, migrator);

  // Return the migrator.
  return migrator;
#else
  return exchanger_new(cloud->comm);
#endif
}

exchanger_t* repartition_mesh(mesh_t* mesh, int* weights, real_t imbalance_tol)
{
  ASSERT(imbalance_tol > 0.0);
  ASSERT(imbalance_tol <= 1.0);
#if POLYMEC_HAVE_MPI
  _Static_assert(sizeof(SCOTCH_Num) == sizeof(index_t), "SCOTCH_Num must be 64-bit.");

  int nprocs, rank;
  MPI_Comm_size(mesh->comm, &nprocs);
  MPI_Comm_rank(mesh->comm, &rank);

  // On a single process, repartitioning has no meaning.
  if (nprocs == 1)
    return exchanger_new(mesh->comm);

  // Generate a local adjacency graph for the mesh.
  adj_graph_t* local_graph = graph_from_mesh_cells(mesh);

  // Get the exchanger for the mesh.
  exchanger_t* mesh_ex = mesh_exchanger(mesh);

  // Now map the distributed graph to the different domains, producing a 
  // global partition vector.
  exchanger_t* migrator = repartition_graph(local_graph, mesh->num_ghost_cells, weights, imbalance_tol, mesh_ex);
  adj_graph_free(local_graph);

  // Migrate the data.
  mesh_migrate(mesh, migrator);

  // Return the migrator.
  return migrator;
#else
  return exchanger_new(mesh->comm);
#endif
}

