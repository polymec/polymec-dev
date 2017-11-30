// Copyright (c) 2012-2017, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_SILO_FILE_H
#define POLYMEC_SILO_FILE_H

#include "core/polymec.h"
#include "core/point_cloud.h"
#include "core/slist.h"
#include "geometry/polymesh.h"
#include "model/neighbor_pairing.h"

// Enables GZIP compression at the given level for Silo files. This is 
// invoked globally and effects all file writes until it is set to a different 
// level. The level should be an integer from 0 to 9, where 0 is the fastest 
// compression possible, and 9 is the maximum (and slowest) level of 
// compression.
void silo_enable_compression(int level);

// This type represents a collection of metadata for field variables in 
// Silo files. Objects of this type are garbage-collected.
typedef struct 
{
  char* label;                // A visualization label. Owned by metadata.
  char* units;                // Units of measure. Owned by metadata.
  bool conserved;             // True if the field is conserved, false if not.
  bool extensive;             // True if the field is extensive, false if not.
  int vector_component;       // The index of the vector component that the 
                              // field represents, or -1 if the field is not 
                              // the component of a vector.
} silo_field_metadata_t;

// Creates a new empty object for storing field metadata.
silo_field_metadata_t* silo_field_metadata_new(void);

// A Silo file can store various geometries (meshes) and data, using 
// "Poor Man's Parallel I/O" (PMPIO) to achieve scalable throughput.
typedef struct silo_file_t silo_file_t;

// Queries the given directory for a file or set of files matching the given 
// prefix, storing the number of files and the number of MPI processes used to 
// write them in the variables num_files and num_mpi_processes. If the steps 
// linked list is non-NULL, it is filled with the step numbers available in 
// the files. This function returns true if the file/files are valid Polymec 
// Silo files (that is, if they were written using the silo_file mechanism)
// and false if they are non-Polymec Silo files or non-Silo files, or if they 
// don't exist.
bool silo_file_query(const char* file_prefix,
                     const char* directory,
                     int* num_files,
                     int* num_mpi_processes,
                     int_slist_t* steps);

// Creates and opens a new Silo file for writing simulation data, 
// returning the Silo file object. If step is non-negative, the file associates 
// itself with the given simulation step number, which is incorporated into 
// its filename. 
// * If directory is the blank string (""), a directory named 
//   <prefix>_<nprocs>procs is generated and used for parallel runs. For 
//   serial runs, the current working directory is used.
// * If the step is -1, the most recent step will be loaded, unless no files 
//   with step information can be found, in which case the single set of files 
//   containing no step information will be loaded. 
// * If the file cannot be created, this function returns NULL.
silo_file_t* silo_file_new(MPI_Comm comm,
                           const char* file_prefix,
                           const char* directory,
                           int num_files,
                           int step,
                           real_t time);

// Opens an existing Silo file for reading simulation data, returning the 
// Silo file object. 
// * If directory is the blank string (""), a directory named 
//   <prefix>_<nprocs>procs is generated and used for parallel runs. For 
//   serial runs, the current working directory is used.
// * If the step is -1, the most recent step will be loaded, unless no files 
//   with step information can be found, in which case the single set of files 
//   containing no step information will be loaded. 
// * If time is not NULL, it will store the time found in the file (or 0.0 if 
//   it does not exist in the file).
// * If the file does not exist or fails to load, this function returns NULL.
silo_file_t* silo_file_open(MPI_Comm comm,
                            const char* file_prefix,
                            const char* directory,
                            int step, 
                            real_t* time);

// Opens an existing Silo file for reading simulation data, masquerading as 
// the given MPI rank in the given communicator. Returns the Silo file object 
// that stores data originally written by the given rank. 
// * If directory is the blank string (""), a directory named 
//   <prefix>_<nprocs>procs is generated and used for parallel runs. For 
//   serial runs, the current working directory is used.
// * If the step is -1, the most recent step will be loaded, unless no files 
//   with step information can be found, in which case the single set of files 
//   containing no step information will be loaded. 
// * If time is not NULL, it will store the time found in the file (or 0.0 if 
//   it does not exist in the file).
// * If the file does not exist or fails to load, this function returns NULL.
silo_file_t* silo_file_open_as_rank(MPI_Comm comm,
                                    int mpi_rank,
                                    const char* file_prefix,
                                    const char* directory,
                                    int step, 
                                    real_t* time);

// Closes and destroys the given Silo file, writing all its data to disk.
void silo_file_close(silo_file_t* file);

// Writes a named arbitrary polyhedral mesh to the given Silo file.
void silo_file_write_polymesh(silo_file_t* file,
                              const char* mesh_name,
                              polymesh_t* mesh);

// Reads a named arbitrary polyhedral mesh from the given Silo file, returning 
// a newly-allocated mesh object.
polymesh_t* silo_file_read_polymesh(silo_file_t* file,
                                    const char* mesh_name);

// Returns true if the given Silo file contains an arbitrary polyhedral mesh 
// with the given name, false otherwise.
bool silo_file_contains_polymesh(silo_file_t* file, const char* mesh_name);

// Writes a named scalar field with the given centering on a polymesh, 
// understood to exist on the polymesh with the given name within the file, 
// to the given Silo file. If a silo_field_metadata object is passed as the 
// last argument, the metadata for the field will also be written to the 
// file--otherwise it can be NULL.
void silo_file_write_scalar_polymesh_field(silo_file_t* file,
                                           const char* field_name,
                                           const char* mesh_name,
                                           real_t* field_data,
                                           polymesh_centering_t centering,
                                           silo_field_metadata_t* field_metadata);

// Writes a named multicomponent field of the given centering on a polymesh, 
// understood to exist on the polymesh with the given name, to the given Silo 
// file. The field data is interpreted to be in component-minor order. If an 
// array of metadata objects is passed as the last argument, the metadata for 
// the field will also be written to the file--otherwise it can be NULL.
void silo_file_write_polymesh_field(silo_file_t* file,
                                    const char** field_component_names,
                                    const char* mesh_name,
                                    real_t* field_data,
                                    int num_components,
                                    polymesh_centering_t centering,
                                    silo_field_metadata_t** field_metadata);

// Reads a named scalar field with the given polymesh centering from the Silo 
// file, returning a newly-allocated array of field data. If a 
// silo_field_metadata object is passed as the last argument, the metadata for 
// the field will be read into it--otherwise it can be NULL.
real_t* silo_file_read_scalar_polymesh_field(silo_file_t* file,
                                             const char* field_name,
                                             const char* mesh_name,
                                             polymesh_centering_t centering,
                                             silo_field_metadata_t* field_metadata);

// Reads a named multicomponent field with the given polymesh centering from 
// the Silo file, returning a newly-allocated array of field data. If an array 
// of metadata objects is passed as the last argument, the metadata for the 
// field will be read into it--otherwise it can be NULL.
real_t* silo_file_read_polymesh_field(silo_file_t* file,
                                      const char** field_component_names,
                                      const char* mesh_name,
                                      int num_components,
                                      polymesh_centering_t centering,
                                      silo_field_metadata_t** field_metadata);

// Returns true if the given Silo file contains a (scalar) field 
// with the given centering and name, and associated with the given polymesh, 
// false otherwise.
bool silo_file_contains_polymesh_field(silo_file_t* file, 
                                       const char* field_name,
                                       const char* mesh_name,
                                       polymesh_centering_t centering);

// Adds a point cloud to the given Silo file.
void silo_file_write_point_cloud(silo_file_t* file,
                                 const char* cloud_name,
                                 point_cloud_t* cloud);

// Reads a point cloud from the Silo file.
point_cloud_t* silo_file_read_point_cloud(silo_file_t* file,
                                          const char* cloud_name);

// Returns true if the given Silo file contains a point cloud 
// with the given name, false otherwise.
bool silo_file_contains_point_cloud(silo_file_t* file, const char* cloud_name);

// Writes a named scalar point field to the given Silo file, associated with 
// the given point cloud. If a silo_field_metadata object is passed as the 
// last argument, the metadata for the field will be read into it--otherwise 
// it can be NULL.
void silo_file_write_scalar_point_field(silo_file_t* file,
                                        const char* field_name,
                                        const char* cloud_name,
                                        real_t* field_data,
                                        silo_field_metadata_t* field_metadata);

// Reads a named scalar point field from the Silo file, 
// allocated array of field data. If a silo_field_metadata object is 
// passed as the last argument, the metadata for the field will be read into 
// it--otherwise it can be NULL.
real_t* silo_file_read_scalar_point_field(silo_file_t* file,
                                          const char* field_name,
                                          const char* cloud_name,
                                          silo_field_metadata_t* field_metadata);

// Writes a named multicomponent point field to the given Silo file, 
// associated with the given point cloud. The field data is interpreted to 
// be in component-minor order. If an array of metadata objects is passed as 
// the last argument, the metadata for the field will also be written to the 
// file--otherwise it can be NULL.
void silo_file_write_point_field(silo_file_t* file,
                                 const char** field_component_names,
                                 const char* cloud_name,
                                 real_t* field_data,
                                 int num_components,
                                 silo_field_metadata_t** field_metadata);

// Reads a named multicomponent point field from the Silo file.
// If an array of metadata objects is passed as the last argument, the 
// metadata for the field will be read into it--otherwise it can be NULL.
real_t* silo_file_read_point_field(silo_file_t* file,
                                   const char** field_component_names,
                                   const char* cloud_name,
                                   int num_components,
                                   silo_field_metadata_t** field_metadata);

// Returns true if the given Silo file contains a (scalar) point field 
// with the given name and associated with the given cloud, false otherwise.
bool silo_file_contains_point_field(silo_file_t* file, 
                                    const char* field_name,
                                    const char* cloud_name);

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

// Return true if the silo file contains a string with the given name, 
// false if not.
bool silo_file_contains_string(silo_file_t* file,
                               const char* string_name);

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
                                size_t array_size);

// Reads an array of real numbers from the silo file, allocating it with 
// polymec_malloc, and filling in the given size. If the array exists but 
// has zero size, NULL will be returned.
real_t* silo_file_read_real_array(silo_file_t* file,
                                  const char* array_name,
                                  size_t* array_size);

// Writes an array of integers of the given size to the silo file. If 
// the array has zero size, an empty entry will be written.
void silo_file_write_int_array(silo_file_t* file,
                               const char* array_name,
                               int* array_data,
                               size_t array_size);

// Reads an array of integers from the silo file, allocating it with 
// polymec_malloc, and filling in the given size. If the array exists, but 
// has zero size, NULL will be returned.
int* silo_file_read_int_array(silo_file_t* file,
                              const char* array_name,
                              size_t* array_size);

// Reads and returns a newly allocated exchanger from the silo file, assigning
// it to the given communicator.
exchanger_t* silo_file_read_exchanger(silo_file_t* file, const char* exchanger_name, MPI_Comm comm);

// Writes an exchanger object with the given name.
void silo_file_write_exchanger(silo_file_t* file, const char* exchanger_name, exchanger_t* ex);

// Returns true if the given silo_file contains a stencil with the given 
// name, false otherwise.
bool silo_file_contains_stencil(silo_file_t* file, const char* stencil_name);

// Writes a stencil object with the given name.
void silo_file_write_stencil(silo_file_t* file,
                             const char* stencil_name,
                             stencil_t* stencil);

// This function extends the silo_file type to allow it to read in and 
// return a newly-allocated stencil object from the entry in the 
// file with the given name. The exchanger for the stencil is assigned
// to the given MPI communicator.
stencil_t* silo_file_read_stencil(silo_file_t* file,
                                  const char* stencil_name,
                                  MPI_Comm comm);

// Returns true if the given silo file contains a neighbor_pairing object 
// with the given name, false otherwise.
bool silo_file_contains_neighbor_pairing(silo_file_t* file,
                                         const char* neighbors_name);

// Writes a neighbor_pairing object to an entry with the given name.
void silo_file_write_neighbor_pairing(silo_file_t* file,
                                      const char* neighbors_name,
                                      neighbor_pairing_t* neighbors);

// This function extends the silo_file type to allow it to read in and 
// return a newly-allocated neighbor_pairing object from the entry in the 
// file with the given name. The exchanger for the neighbor pairing is assigned
// to the given MPI communicator.
neighbor_pairing_t* silo_file_read_neighbor_pairing(silo_file_t* file,
                                                    const char* neighbors_name,
                                                    MPI_Comm comm);

#endif
