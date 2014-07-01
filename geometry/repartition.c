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

exchanger_t* repartition_point_cloud(point_cloud_t* cloud, int* weights)
{
  exchanger_t* ex = exchanger_new(cloud->comm);

#if POLYMEC_HAVE_MPI
  ASSERT(sizeof(SCOTCH_Num) == sizeof(index_t));

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
  int* vtx_dist = adj_graph_vertex_dist(local_graph);
  index_t* global_cell_indices = polymec_malloc(sizeof(SCOTCH_Num) * (cloud->num_points + cloud->num_ghost_points));
  for (int i = 0; i < cloud->num_points + cloud->num_ghost_points; ++i)
    global_cell_indices[i] = (index_t)(vtx_dist[rank] + i);
  exchanger_t* cloud_ex = point_cloud_exchanger(cloud);
  exchanger_exchange(cloud_ex, global_cell_indices, 1, 0, MPI_INT);
  for (int i = 0; i < num_arcs; ++i)
  {
    if (adj[i] >= cloud->num_points)
      adj[i] = global_cell_indices[adj[i]];
  }
  polymec_free(global_cell_indices);

  // Build a distributed graph for the point cloud.
  SCOTCH_Dgraph* dist_graph = SCOTCH_dgraphAlloc();
  SCOTCH_dgraphInit(dist_graph, cloud->comm);
  SCOTCH_dgraphBuild(dist_graph, 0, cloud->num_points, cloud->num_points,
                     xadj, NULL, vert_weights, NULL, num_arcs, num_arcs,
                     adj, NULL, NULL);

  // Free the local graph.
  adj_graph_free(local_graph);

#if 0
  // Now map the distributed graph to the different domains.
  SCOTCH_Num* partition = polymec_malloc(sizeof(SCOTCH_Num) * cloud->num_points);
  int result = SCOTCH_dgraphMap(dist_graph, arch, strategy, partition);
#endif

  // Clean up.
  if (weights != NULL)
    polymec_free(vert_weights);

  SCOTCH_dgraphExit(dist_graph);
  free(dist_graph); // FIXME: ???

  polymec_not_implemented("repartition_mesh");
#endif

  // Return the exchanger.
  return ex;
}

exchanger_t* repartition_mesh(mesh_t* mesh, int* weights)
{
  exchanger_t* ex = exchanger_new(mesh->comm);

#if POLYMEC_HAVE_MPI
  int nprocs, rank;
  MPI_Comm_size(mesh->comm, &nprocs);
  MPI_Comm_rank(mesh->comm, &rank);

  // On a single process, repartitioning has no meaning.
  if (nprocs == 1)
    return ex;

  polymec_not_implemented("repartition_mesh");
#endif
  return ex;
}

