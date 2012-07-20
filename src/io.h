#ifndef ARBI_IO_H
#define ARBI_IO_H

#include "arbi.h"
#include "mesh.h"
#include "lite_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

// The "io" interface type is an opaque type for descriptors used to 
// read and write parallel data.
typedef struct io_interface_t io_interface_t;

// Modes for descriptor operation.
typedef enum
{
  IO_CLOSED,
  IO_READ,
  IO_WRITE
} io_mode_t;

// A function pointer type for opening file descriptors.
typedef int (*io_open_func)(void*, const char*, const char*, io_mode_t, int, MPI_Comm, int, int);

// A function pointer type for closing file descriptors.
typedef int (*io_close_func)(void*);

// Function pointer types for querying meshes from a descriptor.
typedef int (*io_query_meshes_func)(void*, char**, int*);
typedef int (*io_query_lite_meshes_func)(void*, char**, int*);

// A function pointer type for reading a (heavy) mesh from an open descriptor.
typedef int (*io_read_mesh_func)(void*, const char*, mesh_t*);

// A function pointer type for reading a light mesh from an open descriptor.
typedef int (*io_read_lite_mesh_func)(void*, const char*, lite_mesh_t*);

// Function pointer types for querying data fields from a descriptor.
typedef int (*io_query_fields_func)(void*, char**, int*);
typedef int (*io_query_field_func)(void*, const char*, int*, int*);

// A function pointer type for reading a data field from a descriptor.
typedef int (*io_read_field_func)(void*, const char*, double*);

// A function pointer type for writing a (heavy) mesh to an open descriptor.
typedef int (*io_write_mesh_func)(void*, const char*, mesh_t*);

// A function pointer type for writing a light mesh to an open descriptor.
typedef int (*io_write_lite_mesh_func)(void*, const char*, lite_mesh_t*);

// A function pointer type for writing a data field to a descriptor.
typedef int (*io_write_field_func)(void*, const char*, double*, int, int);

// A destructor function for the context object (if any).
typedef void (*io_dtor)(void*);

// This virtual table must be implemented by any I/O interface used 
// within arbi.
typedef struct 
{
  io_open_func              open;
  io_close_func             close;
  io_query_meshes_func      query_meshes;
  io_query_lite_meshes_func query_lite_meshes;
  io_read_mesh_func         read_mesh;
  io_read_lite_mesh_func    read_lite_mesh;
  io_query_fields_func      query_fields;
  io_query_field_func       query_field;
  io_read_field_func        read_field;
  io_write_mesh_func        write_mesh;
  io_write_lite_mesh_func   write_lite_mesh;
  io_write_field_func       write_field;
  io_dtor                   dtor;
} io_vtable;

// Construct an I/O interface (subclass) object from the given name and 
// vtable.
io_interface_t* io_interface_new(void* context, const char* name, io_vtable vtable);

// Frees the given I/O interface.
void io_free(io_interface_t* interface);

// Open a file descriptor with the given interface.
void io_open(io_interface_t* interface, 
             const char* prefix, 
             const char* directory,  
             io_mode_t mode,
             int cycle, 
             MPI_Comm comm,
             int num_files,
             int mpi_tag);

// Close the given file descriptor.
void io_close(io_interface_t* interface);

// Points meshes to an array of names of heavy meshes accessible by the descriptor.
void io_query_meshes(io_interface_t* interface, char** mesh_names, int* num_meshes);

// Reads the heavy mesh with the given name from the descriptor.
void io_read_mesh(io_interface_t* interface, const char* mesh_name, mesh_t* mesh);

// Points lite_meshes to an array of names of heavy meshes accessible by the descriptor.
void io_query_lite_meshes(io_interface_t* interface, char** lite_mesh_names, int* num_lite_meshes);

// A function pointer type for reading a light mesh from an open descriptor.
void io_read_lite_mesh(io_interface_t* interface, const char* lite_mesh_name, lite_mesh_t* lite_mesh);

// Points field_names at an array containing the names of the fields within 
// the file.
void io_query_fields(io_interface_t* interface, char** field_names, int* num_fields);

// Fills num_data and num_components with the values for the given field.
void io_query_field(io_interface_t* interface, const char* field_name, int* num_data, int* num_components);

// Fills field_data with the data from the field.
void io_read_field(io_interface_t* interface, const char* field_name, double* field_data);

// Writes a (heavy) mesh to an open descriptor.
void io_write_mesh(io_interface_t* interface, const char* mesh_name, mesh_t* mesh);

// Writes a light mesh to an open descriptor.
void io_write_lite_mesh(io_interface_t* interface, const char* lite_mesh_name, lite_mesh_t* lite_mesh);

// Writes a field to an open descriptor.
void io_write_field(io_interface_t* interface, const char* field_name, double* field_data, int num_data, int num_components);

#ifdef __cplusplus
}
#endif

#endif

