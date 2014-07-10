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

#ifndef POLYMEC_SILO_FILE_H
#define POLYMEC_SILO_FILE_H

#include "core/polymec.h"
#include "core/mesh.h"
#include "core/unordered_map.h"
#include "core/slist.h"

// A Silo file can store various geometries (meshes) and data, using 
// "Poor Man's Parallel I/O" (PMPIO) to achieve scalable throughput.
typedef struct silo_file_t silo_file_t;

// Queries the given directory for a file or set of files matching the given 
// prefix, storing the number of files and the number of MPI processes used to 
// write them in the variables num_files and num_mpi_processes. If the cycles 
// linked list is non-NULL, it is filled with the cycle numbers available in 
// the files.
void silo_file_query(const char* file_prefix,
                     const char* directory,
                     int* num_files,
                     int* num_mpi_processes,
                     int_slist_t* cycles);

// Creates and opens a new Silo file for writing simulation data, 
// returning the Silo file object. If cycle is non-negative, the file associates 
// itself with the given simulation cycle number, which is incorporated into 
// its filename. If directory is the blank string (""), a directory named 
// <prefix>_<nprocs>procs is generated and used.
silo_file_t* silo_file_new(MPI_Comm comm,
                           const char* file_prefix,
                           const char* directory,
                           int num_files,
                           int mpi_tag,
                           int cycle,
                           real_t time);

// Opens an existing Silo file for reading simulation data, returning the 
// Silo file object. If the directory is the blank string(""), the directory 
// is assumed to be the current working directory. If the cycle is -1, the 
// most recent cycle will be loaded, unless no files with cycle information 
// can be found, in which case the single set of files containing no cycle 
// information will be loaded. If time is not NULL, it will store the time 
// found in the file (or 0.0 if it does not exist in the file).
silo_file_t* silo_file_open(MPI_Comm comm,
                            const char* file_prefix,
                            const char* directory,
                            int mpi_tag,
                            int cycle, 
                            real_t* time);

// Closes and destroys the given Silo file, writing all its data to disk.
void silo_file_close(silo_file_t* file);

// Writes a named arbitrary polyhedral mesh to the given Silo file.
void silo_file_write_mesh(silo_file_t* file,
                          const char* mesh_name,
                          mesh_t* mesh);

// Reads a named arbitrary polyhedral mesh from the given Silo file, returning 
// a newly-allocated mesh object.
mesh_t* silo_file_read_mesh(silo_file_t* file,
                            const char* mesh_name);

// Writes a named scalar cell-centered field, understood to exist on 
// the mesh with the given name, to the given Silo file.
void silo_file_write_scalar_cell_field(silo_file_t* file,
                                       const char* field_name,
                                       const char* mesh_name,
                                       real_t* field_data);

// Reads a named scalar cell-centered field from the Silo file.
real_t* silo_file_read_scalar_cell_field(silo_file_t* file,
                                         const char* field_name,
                                         const char* mesh_name);

// Writes a named multicomponent cell-centered field, understood to exist on 
// the mesh with the given name, to the given Silo file. The field data is 
// interpreted to be in component-minor order.
void silo_file_write_cell_field(silo_file_t* file,
                                const char** field_component_names,
                                const char* mesh_name,
                                real_t* field_data,
                                int num_components);

// Reads a named multicomponent cell-centered field from the Silo file.
real_t* silo_file_read_cell_field(silo_file_t* file,
                                  const char** field_component_names,
                                  const char* mesh_name,
                                  int num_components);

// Adds a point mesh to the given Silo file.
void silo_file_write_point_mesh(silo_file_t* file,
                                const char* point_mesh_name,
                                point_t* points,
                                int num_points);

// Reads a point mesh from the Silo file.
point_t* silo_file_read_point_mesh(silo_file_t* file,
                                   const char* point_mesh_name,
                                   int* num_points);

// Writes a named scalar point field to the given Silo file, associated with 
// the given point mesh.
void silo_file_write_scalar_point_field(silo_file_t* file,
                                        const char* field_name,
                                        const char* point_mesh_name,
                                        real_t* field_data);

// Reads a named scalar point field from the Silo file.
real_t* silo_file_read_scalar_point_field(silo_file_t* file,
                                          const char* field_name,
                                          const char* point_mesh_name);

// Writes a named multicomponent point field to the given Silo file, 
// associated with the given point mesh. The field data is interpreted to 
// be in component-minor order.
void silo_file_write_point_field(silo_file_t* file,
                                 const char** field_component_names,
                                 const char* point_mesh_name,
                                 real_t* field_data,
                                 int num_components);

// Reads a named multicomponent point field from the Silo file.
real_t* silo_file_read_point_field(silo_file_t* file,
                                   const char** field_component_names,
                                   const char* point_mesh_name,
                                   int num_components);

#endif
