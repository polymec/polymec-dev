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

#ifndef POLYMEC_READ_SILO_H
#define POLYMEC_READ_SILO_H

#include "core/polymec.h"
#include "core/mesh.h"
#include "core/unordered_map.h"

// Reads a silo mesh file, including the given fields (whose centerings are 
// assumed to be known to the reader).
void read_silo_mesh(mesh_t** mesh,
                    string_ptr_unordered_map_t** fields,
                    const char* file_prefix,
                    const char* directory,
                    int cycle,
                    double* time,
                    MPI_Comm comm,
                    int num_files,
                    int mpi_tag);

#endif
