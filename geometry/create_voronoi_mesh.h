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

#ifndef POLYMEC_CREATE_VORONOI_MESH_H
#define POLYMEC_CREATE_VORONOI_MESH_H

#include "core/mesh.h"
#include "core/point.h"
#include "core/slist.h"

// This function creates an unbounded Voronoi tessellation of the given points 
// in three-dimensional space. This includes unbounded infinite cells on the 
// boundary of the problem domain.
mesh_t* create_voronoi_mesh(point_t* generators, int num_generators, 
                            point_t* ghost_generators, int num_ghost_generators);

// This function creates an Voronoi tessellation within the given bounding box.
mesh_t* create_voronoi_mesh_in_box(point_t* generators, int num_generators, 
                                   point_t* ghost_generators, int num_ghost_generators,
                                   bbox_t* bounding_box);

#endif

