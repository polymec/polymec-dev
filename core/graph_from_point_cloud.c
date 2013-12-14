// Copyright (c) 2012-2013, Jeffrey N. Johnson
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

#include "core/graph_from_point_cloud.h"

adj_graph_t* graph_from_point_cloud(point_cloud_t* cloud)
{
  // Create a graph whose vertices are the cloud's points.
  adj_graph_t* g = adj_graph_new(cloud->comm, cloud->num_points);

  // Allocate space in the graph for the edges (neighbor relations).
  for (int i = 0; i < cloud->num_points; ++i)
    adj_graph_set_num_edges(g, i, cloud->neighbor_offsets[i+1] - cloud->neighbor_offsets[i]);

  // Now fill in the edges.
  for (int i = 0; i < cloud->num_points; ++i)
  {
    int* edges = adj_graph_edges(g, i);
    for (int j = cloud->neighbor_offsets[i]; j < cloud->neighbor_offsets[i+1]; ++j)
      edges[j-cloud->neighbor_offsets[i]] = cloud->neighbors[j];
  }

  return g;
}

