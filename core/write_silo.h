// Copyright (c) 2012-2014, Jeffrey N. Johnson
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without 
// modification, are permitted provided that the following conditions are met:
// 
// 1. Redistributions of source code must retain the above copyright notice, this 
// list of conditions and the following disclaimer.
// 
// 2. Redistributions in binary form must reproduce the above copyright notice, 
// this list of conditions and the following disclaimer in the documentation 
// and/or other materials provided with the distribution.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE 
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#ifndef POLYMEC_WRITE_SILO_H
#define POLYMEC_WRITE_SILO_H

#include "core/polymec.h"
#include "core/mesh.h"
#include "core/unordered_map.h"

// Writes a silo mesh file, with the given named cell-centered fields.
void write_silo_mesh(mesh_t* mesh,
                     string_ptr_unordered_map_t* fields,
                     const char* file_prefix,
                     const char* directory,
                     int cycle,
                     real_t time,
                     MPI_Comm comm,
                     int num_files,
                     int mpi_tag);

// Writes a silo point file, with the given fields defined on the points.
void write_silo_points(point_t* points,
                       int num_points,
                       string_ptr_unordered_map_t* fields,
                       const char* file_prefix,
                       const char* directory,
                       int cycle,
                       real_t time,
                       MPI_Comm comm,
                       int num_files,
                       int mpi_tag);

#endif
