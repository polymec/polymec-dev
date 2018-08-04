// Copyright (c) 2012-2018, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_VORONOI2_H
#define POLYMEC_VORONOI2_H

#include "core/adj_graph.h"
#include "geometry/polygon.h"

/// \addtogroup geometry geometry
///@{

/// \class voronoi2
/// This class contains data defining a 2D Vornoi diagram in the plane.
struct voronoi2_t
{
  /// Points used to generate the diagram.
  point2_t* generators;
  /// The number of points used to generate the diagram.
  size_t num_generators;
  /// The polygonal cells in the diagram.
  polygon_t** cells;
  /// An adjacency graph that defines the neighbor relations for the polygons.
  adj_graph_t* graph;
};
typedef struct voronoi2_t voronoi2_t;

/// Creates a 2D Voronoi diagram with the given generators, cells, and graph.
/// This does not compute the diagram--it only allocates storage for the 
/// container and assumes control of all resources passed as input.
/// \memberof voronoi2
voronoi2_t* voronoi2_new(point2_t* generators,
                         size_t num_generators,
                         polygon_t** cells,
                         adj_graph_t* graph);

/// Frees the 2D Voronoi diagram and all its resources.
/// \memberof voronoi2
void voronoi2_free(voronoi2_t* diagram);

#endif

