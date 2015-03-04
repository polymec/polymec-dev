// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_SILO_FILE_H
#define POLYMEC_SILO_FILE_H

#include "core/polymec.h"
#include "core/mesh.h"
#include "core/point_cloud.h"
#include "core/slist.h"

// A Silo file can store various geometries (meshes) and data, using 
// "Poor Man's Parallel I/O" (PMPIO) to achieve scalable throughput.
typedef struct silo_file_t silo_file_t;

// Queries the given directory for a file or set of files matching the given 
// prefix, storing the number of files and the number of MPI processes used to 
// write them in the variables num_files and num_mpi_processes. If the cycles 
// linked list is non-NULL, it is filled with the cycle numbers available in 
// the files. This function returns true if the file/files are valid Polymec 
// Silo files (that is, if they were written using the silo_file mechanism)
// and false if they are non-Polymec Silo files or non-Silo files.
bool silo_file_query(const char* file_prefix,
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

// Reads a named scalar cell-centered field from the Silo file, returning a newly-
// allocated array of field data.
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

// Reads a named multicomponent cell-centered field from the Silo file, returning a 
// newly-allocated array of field data.
real_t* silo_file_read_cell_field(silo_file_t* file,
                                  const char** field_component_names,
                                  const char* mesh_name,
                                  int num_components);

// Writes a named scalar face-centered field, understood to exist on 
// the mesh with the given name, to the given Silo file.
void silo_file_write_scalar_face_field(silo_file_t* file,
                                       const char* field_name,
                                       const char* mesh_name,
                                       real_t* field_data);

// Reads a named scalar face-centered field from the Silo file, returning a newly-
// allocated array of field data.
real_t* silo_file_read_scalar_face_field(silo_file_t* file,
                                         const char* field_name,
                                         const char* mesh_name);

// Writes a named multicomponent face-centered field, understood to exist on 
// the mesh with the given name, to the given Silo file. The field data is 
// interpreted to be in component-minor order.
void silo_file_write_face_field(silo_file_t* file,
                                const char** field_component_names,
                                const char* mesh_name,
                                real_t* field_data,
                                int num_components);

// Reads a named multicomponent cell-centered field from the Silo file, returning a 
// newly-allocated array of field data.
real_t* silo_file_read_face_field(silo_file_t* file,
                                  const char** field_component_names,
                                  const char* mesh_name,
                                  int num_components);

// Adds a point mesh to the given Silo file.
void silo_file_write_point_cloud(silo_file_t* file,
                                 const char* cloud_name,
                                 point_cloud_t* cloud);

// Reads a point mesh from the Silo file.
point_cloud_t* silo_file_read_point_cloud(silo_file_t* file,
                                          const char* cloud_name);

// Writes a named scalar point field to the given Silo file, associated with 
// the given point cloud.
void silo_file_write_scalar_point_field(silo_file_t* file,
                                        const char* field_name,
                                        const char* cloud_name,
                                        real_t* field_data);

// Reads a named scalar point field from the Silo file.
real_t* silo_file_read_scalar_point_field(silo_file_t* file,
                                          const char* field_name,
                                          const char* cloud_name);

// Writes a named multicomponent point field to the given Silo file, 
// associated with the given point cloud. The field data is interpreted to 
// be in component-minor order.
void silo_file_write_point_field(silo_file_t* file,
                                 const char** field_component_names,
                                 const char* cloud_name,
                                 real_t* field_data,
                                 int num_components);

// Reads a named multicomponent point field from the Silo file.
real_t* silo_file_read_point_field(silo_file_t* file,
                                   const char** field_component_names,
                                   const char* cloud_name,
                                   int num_components);

// Adds a scalar expression (a definition of a new scalar variable in terms 
// of existing variables) to this Silo file. See the Silo manual for
// expression (Defvar) syntax.
void silo_file_write_scalar_expression(silo_file_t* file,
                                       const char* expression_name,
                                       const char* definition);

// Adds a vector expression (a definition of a new vector variable in terms 
// of existing variables) to this Silo file. See the Silo manual for 
// expression (Defvar) syntax.
void silo_file_write_vector_expression(silo_file_t* file,
                                       const char* expression_name,
                                       const char* definition);

// Adds a (rank 2) tensor expression (a definition of a new tensor variable 
// in terms of existing variables) to this Silo file. See the Silo manual for 
// expression (Defvar) syntax.
void silo_file_write_tensor_expression(silo_file_t* file,
                                       const char* expression_name,
                                       const char* definition);

// Writes a named string to the silo file.
void silo_file_write_string(silo_file_t* file,
                            const char* string_name,
                            char* string_data);

// Reads a named string from the silo file, returning a newly-allocated copy.
char* silo_file_read_string(silo_file_t* file,
                            const char* string_name);

// Writes an array of real numbers of the given size to the silo file.
// If the array has zero size, an empty entry will be written.
void silo_file_write_real_array(silo_file_t* file,
                                const char* array_name,
                                real_t* array_data,
                                int array_size);

// Reads an array of real numbers from the silo file, allocating it with 
// polymec_malloc, and filling in the given size. If the array exists but 
// has zero size, NULL will be returned.
real_t* silo_file_read_real_array(silo_file_t* file,
                                  const char* array_name,
                                  int* array_size);

// Writes an array of integers of the given size to the silo file. If 
// the array has zero size, an empty entry will be written.
void silo_file_write_int_array(silo_file_t* file,
                               const char* array_name,
                               int* array_data,
                               int array_size);

// Reads an array of integers from the silo file, allocating it with 
// polymec_malloc, and filling in the given size. If the array exists, but 
// has zero size, NULL will be returned.
int* silo_file_read_int_array(silo_file_t* file,
                              const char* array_name,
                              int* array_size);

// Reads and returns a newly allocated exchanger from the silo file, assigning
// it to the given communicator.
exchanger_t* silo_file_read_exchanger(silo_file_t* file, const char* exchanger_name, MPI_Comm comm);

// Writes an exchanger object with the given name.
void silo_file_write_exchanger(silo_file_t* file, const char* exchanger_name, exchanger_t* ex);

#endif
