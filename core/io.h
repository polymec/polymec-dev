#ifndef ARBI_IO_H
#define ARBI_IO_H

#include "core/arbi.h"
#include "core/mesh.h"
#include "core/lite_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

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

// A function pointer type for opening file descriptors.
typedef int (*io_open_func)(void*, const char*, const char*, io_mode_t, MPI_Comm, int, int);

// A function pointer type for closing file descriptors.
typedef int (*io_close_func)(void*);

// A function pointer type for reading a (heavy) mesh from an open descriptor.
typedef int (*io_read_mesh_func)(void*, const char*, mesh_t**);

// A function pointer type for reading a light mesh from an open descriptor.
typedef int (*io_read_lite_mesh_func)(void*, const char*, lite_mesh_t**);

// Function pointer types for querying data fields from a descriptor.
typedef int (*io_query_field_func)(void*, const char*, const char*, int*, mesh_centering_t*);

// A function pointer type for reading a data field from a descriptor.
typedef int (*io_read_field_func)(void*, const char*, const char*, double*);

// A function pointer for reading a block of source code from a descriptor.
typedef int (*io_query_source_code_func)(void*, const char*, const char*, int*);

// A function pointer for reading a block of source code from a descriptor.
typedef int (*io_read_source_code_func)(void*, const char*, const char*, char*);

// A function pointer type for writing a (heavy) mesh to an open descriptor.
typedef int (*io_write_mesh_func)(void*, const char*, mesh_t*);

// A function pointer type for writing a light mesh to an open descriptor.
typedef int (*io_write_lite_mesh_func)(void*, const char*, lite_mesh_t*);

// A function pointer type for writing a data field to a descriptor.
typedef int (*io_write_field_func)(void*, const char*, const char*, double*, int, mesh_centering_t);

// A function pointer for writing a block of source code from a descriptor.
typedef int (*io_write_source_code_func)(void*, const char*, const char*, const char*);

// A destructor function for the context object (if any).
typedef void (*io_dtor)(void*);

// This virtual table must be implemented by any I/O interface used 
// within arbi.
typedef struct 
{
  io_open_func                  open;
  io_close_func                 close;
  io_read_mesh_func             read_mesh;
  io_read_lite_mesh_func        read_lite_mesh;
  io_query_field_func           query_field;
  io_read_field_func            read_field;
  io_query_source_code_func     query_source_code;
  io_read_source_code_func      read_source_code;
  io_write_mesh_func            write_mesh;
  io_write_lite_mesh_func       write_lite_mesh;
  io_write_field_func           write_field;
  io_write_source_code_func     write_source_code;
  io_dtor                       dtor;
} io_vtable;

// Construct an I/O interface (subclass) object from the given name and 
// vtable.
io_interface_t* io_interface_new(void* context, const char* name, io_vtable vtable);

// Frees the given I/O interface.
void io_free(io_interface_t* interface);

// Open an I/O descriptor with the given interface.
void io_open(io_interface_t* interface, 
             const char* prefix, 
             const char* directory,  
             io_mode_t mode,
             MPI_Comm comm,
             int num_files,
             int mpi_tag);

// Close the given file descriptor.
void io_close(io_interface_t* interface);

// Returns a dataset descriptor that can be used to retrieve data from the
// open file descriptor. If there is no dataset by the given name, this 
// returns NULL.
io_dataset_t* io_dataset(io_interface_t* interface, const char* dataset);

// Frees the given dataset descriptor.
void io_dataset_free(io_dataset_t* dataset);

// Reads a heavy mesh from the descriptor.
void io_dataset_read_mesh(io_dataset_t* dataset, mesh_t** mesh);

// Writes a heavy mesh the descriptor.
void io_dataset_write_mesh(io_dataset_t* dataset, mesh_t* mesh);

// Reads a lite mesh from the descriptor.
void io_dataset_read_lite_mesh(io_dataset_t* dataset, lite_mesh_t** mesh);

// Writes a lite mesh to the descriptor.
void io_dataset_write_lite_mesh(io_dataset_t* dataset, lite_mesh_t* mesh);

// Gathers metadata about the given field in the dataset. The size of the 
// field is consistent with its associated mesh.
void io_dataset_query_field(io_dataset_t* dataset, const char* field_name, int* num_components, mesh_centering_t* centering);

// Reads field data from the descriptor.
void io_dataset_read_field(io_dataset_t* dataset, const char* field_name, double* field_data);

// Writes field data to the descriptor.
void io_dataset_write_field(io_dataset_t* dataset, const char* field_name, double* field_data, int num_components, mesh_centering_t centering);

// Queries the dataset descriptor for a named block of source code, retrieving its length.
void io_dataset_query_source_code(io_dataset_t* dataset, const char* code_name, int* len);

// Reads a named block of source code from the dataset descriptor.
void io_dataset_read_source_code(io_dataset_t* dataset, const char* code_name, char* source_code);

// Writes a named block of source code to the dataset descriptor.
void io_dataset_write_source_code(io_dataset_t* dataset, const char* code_name, const char* source_code);

#ifdef __cplusplus
}
#endif

#endif

