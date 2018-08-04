// Copyright (c) 2012-2018, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_VORONOI3_H
#define POLYMEC_VORONOI3_H

#include "core/adj_graph.h"
#include "geometry/polyhedron.h"

/// \addtogroup geometry geometry
///@{

/// \struct voronoi3
/// This struct contains data defining a 3D Vornoi diagram.
struct voronoi3_t
{
  /// Points used to generate the diagram.
  point_t* generators;
  /// The number of points used to generate the diagram.
  size_t num_generators;
  /// The polyhedral cells in the diagram.
  polyhedron_t** cells;
  /// An adjacency graph that defines the neighbor relations for the cells,
  /// and also their parallel distribution.
  adj_graph_t* graph;
};
typedef struct voronoi3_t voronoi3_t;

/// Creates a 3D Voronoi diagram with the given generators, cells, and graph.
/// This does not compute the diagram--it only allocates storage for a container
/// and assumes control of all resources passed as input.
/// \memberof voronoi3
voronoi3_t* voronoi3_new(point_t* generators,
                         size_t num_generators,
                         polyhedron_t** cells,
                         adj_graph_t* graph);

/// Frees the 3D Voronoi diagram and all its resources.
/// \memberof voronoi3
void voronoi3_free(voronoi3_t* diagram);

#endif

