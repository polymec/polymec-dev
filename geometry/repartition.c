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
#endif

exchanger_t* repartition_point_cloud(point_cloud_t* cloud, int* weights, real_t imbalance_tol)
{
  exchanger_t* ex = exchanger_new(cloud->comm);

#if POLYMEC_HAVE_MPI
  _Static_assert(sizeof(SCOTCH_Num) == sizeof(index_t), "SCOTCH_Num must be 64 bits.");

  int nprocs, rank;
  MPI_Comm_size(cloud->comm, &nprocs);
  MPI_Comm_rank(cloud->comm, &rank);

  // On a single process, repartitioning has no meaning.
  if (nprocs == 1)
    return ex;

  SCOTCH_Num* vert_weights = NULL;
  if (weights != NULL)
  {
    vert_weights = polymec_malloc(sizeof(SCOTCH_Num) * cloud->num_points);
    for (int i = 0; i < cloud->num_points; ++i)
      vert_weights[i] = (SCOTCH_Num)weights[i];
  }

  // Generate a local adjacency graph for the point cloud.
  adj_graph_t* local_graph = graph_from_point_cloud(cloud);

  // Extract the adjacency information.
  SCOTCH_Num* xadj = malloc(sizeof(SCOTCH_Num) * (cloud->num_points+1));
  memcpy(xadj, adj_graph_edge_offsets(local_graph), sizeof(SCOTCH_Num) * cloud->num_points+1);
  SCOTCH_Num num_arcs = xadj[cloud->num_points];
  SCOTCH_Num* adj = malloc(sizeof(SCOTCH_Num) * num_arcs);
  memcpy(adj, adj_graph_adjacency(local_graph), sizeof(SCOTCH_Num) * num_arcs);

  // Populate adj with global cell indices.
  index_t* vtx_dist = adj_graph_vertex_dist(local_graph);
  index_t* global_point_indices = polymec_malloc(sizeof(SCOTCH_Num) * (cloud->num_points + cloud->num_ghost_points));
  for (int i = 0; i < cloud->num_points + cloud->num_ghost_points; ++i)
    global_point_indices[i] = (index_t)(vtx_dist[rank] + i);
  exchanger_t* cloud_ex = point_cloud_exchanger(cloud);
  exchanger_exchange(cloud_ex, global_point_indices, 1, 0, MPI_UINT64_T);
  for (int i = 0; i < num_arcs; ++i)
  {
    if (adj[i] >= cloud->num_points)
      adj[i] = global_point_indices[adj[i]];
  }
  polymec_free(global_point_indices);
printf("adj = [");
for (int i = 0; i < num_arcs; ++i)
printf("%d ", adj[i]);
printf("]\n");

  // Build a distributed graph for the point cloud.
  SCOTCH_Num num_points = cloud->num_points;
  SCOTCH_Dgraph dist_graph;
  SCOTCH_dgraphInit(&dist_graph, cloud->comm);
  SCOTCH_dgraphBuild(&dist_graph, 0, num_points, num_points,
                     xadj, NULL, vert_weights, NULL, num_arcs, num_arcs,
                     adj, NULL, NULL);

  // Free the local graph.
  adj_graph_free(local_graph);

  // Now map the distributed graph to the different domains.
  SCOTCH_Arch arch;
  SCOTCH_archInit(&arch);
  if (weights == NULL)
    SCOTCH_archCmplt(&arch, num_points);
  else
    SCOTCH_archCmpltw(&arch, num_points, vert_weights);
  SCOTCH_Strat strategy;
  SCOTCH_stratInit(&strategy);
  SCOTCH_Num strat_flags = SCOTCH_STRATDEFAULT;
  SCOTCH_stratDgraphMapBuild(&strategy, strat_flags, nprocs, nprocs, (double)imbalance_tol);
  SCOTCH_Num* partition = polymec_malloc(sizeof(SCOTCH_Num) * cloud->num_points);
  int result = SCOTCH_dgraphMap(&dist_graph, &arch, &strategy, partition);
  if (result != 0)
    polymec_error("Repartitioning failed.");

printf("P = [");
for (int i = 0; i < cloud->num_points; ++i)
printf("%d ", partition[i]);
printf("]\n");
  // FIXME: Transfer data!

  // Clean up.
  if (weights != NULL)
    polymec_free(vert_weights);
  polymec_free(partition);
  SCOTCH_stratExit(&strategy);
  SCOTCH_archExit(&arch);
  SCOTCH_dgraphExit(&dist_graph);

#endif

  // Return the exchanger.
  return ex;
}

exchanger_t* repartition_mesh(mesh_t* mesh, int* weights, real_t imbalance_tol)
{
  ASSERT(imbalance_tol > 0.0);
  ASSERT(imbalance_tol <= 1.0);
  exchanger_t* ex = exchanger_new(mesh->comm);

#if POLYMEC_HAVE_MPI
  _Static_assert(sizeof(SCOTCH_Num) == sizeof(index_t), "SCOTCH_Num must be 64 bits.");

  int nprocs, rank;
  MPI_Comm_size(mesh->comm, &nprocs);
  MPI_Comm_rank(mesh->comm, &rank);

  // On a single process, repartitioning has no meaning.
  if (nprocs == 1)
    return ex;

  SCOTCH_Num* vert_weights = NULL;
  if (weights != NULL)
  {
    vert_weights = polymec_malloc(sizeof(SCOTCH_Num) * mesh->num_cells);
    for (int i = 0; i < mesh->num_cells; ++i)
      vert_weights[i] = (SCOTCH_Num)weights[i];
  }

  // Generate a local adjacency graph for the mesh.
  adj_graph_t* local_graph = graph_from_mesh_cells(mesh);
for (int p = 0; p < nprocs; ++p)
{
if (rank == p)
  adj_graph_fprintf(local_graph, stdout);
MPI_Barrier(mesh->comm);
}

  // Extract the adjacency information.
  SCOTCH_Num* xadj = malloc(sizeof(SCOTCH_Num) * (mesh->num_cells+1));
  int* edge_offsets = adj_graph_edge_offsets(local_graph);
  for (int i = 0; i <= mesh->num_cells; ++i)
    xadj[i] = (SCOTCH_Num)edge_offsets[i];
  SCOTCH_Num num_arcs = xadj[mesh->num_cells];
  SCOTCH_Num* adj = malloc(sizeof(SCOTCH_Num) * num_arcs);
  int* edges = adj_graph_adjacency(local_graph);
  for (int i = 0; i < xadj[mesh->num_cells]; ++i)
    adj[i] = (SCOTCH_Num)edges[i];

  // Replace the ghost entries in adj with global cell indices.
  index_t* vtx_dist = adj_graph_vertex_dist(local_graph);
  index_t* global_cell_indices = polymec_malloc(sizeof(SCOTCH_Num) * (mesh->num_cells + mesh->num_ghost_cells));
  for (int i = 0; i < mesh->num_cells; ++i)
    global_cell_indices[i] = (index_t)(vtx_dist[rank] + i);
  exchanger_t* mesh_ex = mesh_exchanger(mesh);
exchanger_fprintf(mesh_ex, stdout);
exchanger_enable_deadlock_detection(mesh_ex, 1.0, 0, stdout);
  exchanger_exchange(mesh_ex, global_cell_indices, 1, 0, MPI_UINT64_T);
for (int i = 0; i < mesh->num_cells + mesh->num_ghost_cells; ++i)
printf("%d: GCI[%d] = %d\n", rank, i, global_cell_indices[i]);
  for (int i = 0; i < num_arcs; ++i)
    adj[i] = global_cell_indices[adj[i]];
  polymec_free(global_cell_indices);
printf("%d: adj = [", rank);
for (int i = 0; i < num_arcs; ++i)
printf("%d ", adj[i]);
printf("]\n");

  // Build a distributed graph for the mesh.
  SCOTCH_Num num_cells = mesh->num_cells;
  SCOTCH_Dgraph dist_graph;
  SCOTCH_dgraphInit(&dist_graph, mesh->comm);
  SCOTCH_dgraphBuild(&dist_graph, 0, num_cells, num_cells,
                     xadj, NULL, vert_weights, NULL, num_arcs, num_arcs,
                     adj, NULL, NULL);

  // Now map the distributed graph to the different domains.
  SCOTCH_Arch arch;
  SCOTCH_archInit(&arch);
  if (weights == NULL)
    SCOTCH_archCmplt(&arch, num_cells);
  else
    SCOTCH_archCmpltw(&arch, num_cells, vert_weights);
  SCOTCH_Strat strategy;
  SCOTCH_stratInit(&strategy);
  SCOTCH_Num strat_flags = SCOTCH_STRATDEFAULT;
  SCOTCH_stratDgraphMapBuild(&strategy, strat_flags, nprocs, nprocs, (double)imbalance_tol);
  SCOTCH_Num* partition = polymec_malloc(sizeof(SCOTCH_Num) * mesh->num_cells);
  int result = SCOTCH_dgraphMap(&dist_graph, &arch, &strategy, partition);
  if (result != 0)
    polymec_error("Repartitioning failed.");

printf("P = [");
for (int i = 0; i < mesh->num_cells; ++i)
printf("%d ", partition[i]);
printf("]\n");
  // FIXME: Transfer data!

  // Clean up.
  if (weights != NULL)
    polymec_free(vert_weights);
  polymec_free(partition);
  SCOTCH_stratExit(&strategy);
  SCOTCH_archExit(&arch);
  SCOTCH_dgraphExit(&dist_graph);
  adj_graph_free(local_graph);

#endif

  return ex;
}

