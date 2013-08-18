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

#ifndef POLYMEC_REPARTITION_H
#define POLYMEC_REPARTITION_H

#include "core/point.h"
#include "core/exchanger.h"
#include "core/mesh.h"

// This function repartitions the given set of points with the given 
// weights, alloting them to parallel domains to balance their load.
// If weights is NULL, the points are assigned equal weights.
// It creates and returns an exchanger object that can be used to migrate 
// data from the old partition to the new.
exchanger_t* repartition_points(MPI_Comm comm, point_t* points, double* weights, int num_points);

// This function repartitions the given mesh, alloting the cells to parallel 
// domains to balance the load.
// It creates and returns an exchanger object that can be used to migrate 
// data from the old partition to the new.
exchanger_t* repartition_mesh(MPI_Comm comm, mesh_t* mesh);

#endif

