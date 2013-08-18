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

// This function creates a Voronoi tessellation of the given points in 
// three-dimensional space. The tessellation contains only bounded cells.
// This means that generators corresponding to unbounded cells do not appear 
// in the tessellation, so it's important to include any generators that 
// produce a non-zero set of bounded Voronoi cells. Optionally, a linked list 
// may be provided that will store the indices of any cells that were deleted to construct a 
// completely bounded mesh. If such information is not desired, the final 
// argument can be set to NULL.
mesh_t* create_voronoi_mesh(point_t* generators, int num_generators, 
                            point_t* ghost_generators, int num_ghost_generators,
                            int_slist_t* deleted_cells);

#endif

