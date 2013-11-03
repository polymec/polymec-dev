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

#ifndef POLYMEC_CREATE_UNIFORM_MESH_H
#define POLYMEC_CREATE_UNIFORM_MESH_H

#include "core/mesh.h"

// This function creates and returns a uniform mesh of nx x ny x nz cells. 
// The mesh spans the rectangular region of space defined by the bounding box.
mesh_t* create_uniform_mesh(int nx, int ny, int nz, bbox_t* bbox);

#endif

