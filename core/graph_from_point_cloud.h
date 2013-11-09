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

#ifndef POLYMEC_GRAPH_FROM_POINT_CLOUD_H
#define POLYMEC_GRAPH_FROM_POINT_CLOUD_H

#include "core/polymec.h"
#include "core/point_cloud.h"
#include "core/adj_graph.h"

// Creates an adjacency graph representing the topological relationships 
// between the points in the given point cloud.
adj_graph_t* graph_from_point_cloud(point_cloud_t* cloud);

#endif 

