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

#ifndef POLYMEC_IO_H
#define POLYMEC_IO_H

#include "core/polymec.h"
#include "core/mesh.h"

// The "io" interface type is an opaque type for descriptors used to 
// read and write parallel data.
typedef struct io_interface_t io_interface_t;

// Modes for mesh descriptor operation.
typedef enum
{
  IO_CLOSED,
  IO_READ,
  IO_WRITE
} io_mode_t;

// This descriptor allows one to retrieve data for a given dataset.
typedef struct io_dataset_t io_dataset_t;

// A function pointer type for creating file descriptors.
typedef void* (*io_create_file_func)(void*, const char*, const char*);

// A function pointer type for opening file descriptors.
typedef void* (*io_open_file_func)(void*, const char* , const char*, io_mode_t);

// A function pointer type for closing file descriptors.
typedef void (*io_close_file_func)(void*, void*);

// A function pointer for determining the number of datasets in a file.
typedef int (*io_get_num_datasets_func)(void*, void*, int*);

// A function pointer for reading data from all datasets in a file.
typedef void (*io_read_datasets_func)(void*, void*, io_dataset_t**, int);

// A function pointer for dumping data to datasets in a file.
typedef void (*io_write_datasets_func)(void*, void*, io_dataset_t**, int, int, int);

// A function pointer for writing a master file if needed.
typedef void (*io_write_master_func)(void*, void*, const char*, io_dataset_t**, int, int, int);

// A destructor function for the context object (if any).
typedef void (*io_dtor)(void*);

// This virtual table must be implemented by any I/O interface used 
// within polymec.
typedef struct 
{
  io_create_file_func           create_file;
  io_open_file_func             open_file;
  io_close_file_func            close_file;
  io_get_num_datasets_func      get_num_datasets;
  io_read_datasets_func         read_datasets;
  io_write_datasets_func        write_datasets;
  io_write_master_func          write_master;
  io_dtor                       dtor;
} io_vtable;

// Construct an I/O interface (subclass) object from the given name and 
// vtable. To take advantage of "poor man's parallel I/O", one can provide
// an MPI communicator, a number of files to be written, and an MPI tag.
io_interface_t* io_interface_new(void* context, 
                                 const char* name, 
                                 const char* suffix,
                                 const char* master_suffix,
                                 io_vtable vtable,
                                 MPI_Comm comm,
                                 int num_files,
                                 int mpi_tag);

// Construct a serial I/O interface (subclass) object from the given name and 
// vtable. 
io_interface_t* io_interface_new_serial(void* context, 
                                        const char* name, 
                                        const char* suffix,
                                        io_vtable vtable);

// Frees the given I/O interface.
void io_free(io_interface_t* interface);

// Open an I/O descriptor with the given interface.
void io_open(io_interface_t* interface, 
             const char* prefix, 
             const char* directory,  
             io_mode_t mode);

// Close the given file descriptor.
void io_close(io_interface_t* interface);

// Returns the number of datasets available via the given interface.
// This should be used when reading from an interface.
int io_num_datasets(io_interface_t* interface);

// Sets the number of datasets available via the given interface.
// This should be used when writing to an interface.
void io_set_num_datasets(io_interface_t* interface, int num_datasets);

// Appends the given dataset to the end of the list of datasets in the 
// given interface. The io_interface assumes control over dataset.
void io_append_dataset(io_interface_t* interface, io_dataset_t* dataset);

// Returns the default dataset descriptor for the file. Use this to retrieve
// data from files that can hold only one dataset.
io_dataset_t* io_default_dataset(io_interface_t* interface);

// Returns a dataset descriptor that can be used to retrieve data from the
// open file descriptor. If there is no dataset by the given name, this 
// returns NULL.
io_dataset_t* io_dataset(io_interface_t* interface, const char* dataset);

// Creates a new dataset that can be written to a file.
io_dataset_t* io_dataset_new(const char* name);

// Returns the name of the given dataset.
const char* io_dataset_name(io_dataset_t* dataset);

// Frees the given dataset descriptor.
void io_dataset_free(io_dataset_t* dataset);

// Retrieves a mesh from the descriptor. The caller assumes responsibility for 
// the storage of the mesh.
mesh_t* io_dataset_get_mesh(io_dataset_t* dataset);

// Writes a mesh the descriptor.
void io_dataset_put_mesh(io_dataset_t* dataset, mesh_t* mesh);

// Gathers metadata about the given field in the dataset. The size of the 
// field is consistent with its associated mesh.
void io_dataset_query_field(io_dataset_t* dataset, const char* field_name, int* num_components, mesh_centering_t* centering);

// Accesses field data from the descriptor. The caller assumes responsibility 
// for the field storage after this call.
void io_dataset_get_field(io_dataset_t* dataset, const char* field_name, double** field);

// Copies or passes a field to the descriptor. If the field is passed to 
// the descriptor, it will be destroyed by it subsequently.
void io_dataset_put_field(io_dataset_t* dataset, const char* field_name, double* field_data, int num_components, mesh_centering_t centering, bool copy);

// Allows iteration over the fields in the given dataset.
bool io_dataset_next_field(io_dataset_t* dataset, int* pos, char** field_name, double** field, int* num_components, mesh_centering_t* centering);

// Returns the number of fields stored in the dataset.
int io_dataset_num_fields(io_dataset_t* dataset);

// Gets an internally-stored string from the dataset descriptor, or NULL
// if no such string exists. Must be copied by the caller in order to be used.
char* io_dataset_get_string(io_dataset_t* dataset, const char* string_name);

// Copies a named string to the dataset descriptor.
void io_dataset_put_string(io_dataset_t* dataset, const char* string_name, const char* string);

// Allows iteration over the strings in the given dataset.
bool io_dataset_next_string(io_dataset_t* dataset, int* pos, char** string_name, char** string);

// Returns the number of strings stored in the dataset.
int io_dataset_num_strings(io_dataset_t* dataset);

#endif

