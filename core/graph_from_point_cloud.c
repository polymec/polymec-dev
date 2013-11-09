// Copyright 2012-2013 Jeffrey Johnson.
// 
// This file is part of Polymec, and is licensed under the Apache License, 
// Version 2.0 (the "License"); you may not use this file except in 
// compliance with the License. You may may find the text of the license in 
// the LICENSE file at the top-level source directory, or obtain a copy of 
// it at
// 
// http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

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

