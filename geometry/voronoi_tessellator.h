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

#ifndef POLYMEC_VORONOI_TESSELLATOR_H
#define POLYMEC_VORONOI_TESSELLATOR_H

#include "core/slist.h"
#include "core/point.h"
#include "core/point2.h"

// This represents a Voronoi cell in a tessellation.
typedef struct
{
  int num_faces;
  int* faces;
} voronoi_cell_t;

// This represents a Voronoi face separating two cells.
typedef struct
{
  int num_edges;
  int* edges;
  int cell1, cell2;
} voronoi_face_t;

// This represents a Voronoi edge connecting two nodes.
typedef struct
{
  int node1, node2;
} voronoi_edge_t;

// This type holds a three-dimensional Voronoi tessellation 
// produced by a tessellator.
typedef struct
{
  int num_cells;
  voronoi_cell_t* cells;

  int num_faces;
  voronoi_face_t* faces;

  int num_edges;
  voronoi_edge_t* edges;

  int num_nodes;
  double* nodes;

} voronoi_tessellation_t;

// The tessellator class creates Voronoi tessellations from sets of 
// generator points. Objects of this type are garbage-collected.
typedef struct voronoi_tessellator_t voronoi_tessellator_t;

// Constructs a new Voronoi tessellator.
voronoi_tessellator_t* voronoi_tessellator_new();

// Creates a three-dimensional tessellation with the given set of generators. 
// The tessellation is bounded--infinite cells are removed.
// The last argument is an optionally provided (otherwise NULL) linked list
// that gathers points that are deleted to construct a completely bounded
// tessellation.
voronoi_tessellation_t* 
voronoi_tessellator_tessellate(voronoi_tessellator_t* tessellator,
                               point_t* points, int num_points,
                               int_slist_t* deleted_points);

// Creates a planar tessellation with the given set of generators. 
// The tessellation is bounded by a bounding polygon consisting of a 
// set of vertices traversed in the order given. Note that in a 2D planar
// Voronoi graph, each face corresponds to a single edge.
voronoi_tessellation_t* 
voronoi_tessellator_tessellate_2d(voronoi_tessellator_t* tessellator,
                                  point2_t* points, int num_points,
                                  point2_t* bounding_polygon, int num_bounding_edges);

// Destroys a tessellation that has been created by a tessellator.
void voronoi_tessellation_free(voronoi_tessellation_t* tessellation);

#endif

