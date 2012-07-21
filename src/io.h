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

// Modes for mesh descriptor operation.
typedef enum
{
  IO_CLOSED,
  IO_READ,
  IO_WRITE
} io_mode_t;

// This iterator allows one to traverse the data sets on heavy meshes 
// in a descriptor.
typedef struct io_mesh_iter io_mesh_iter;

// This iterator allows one to traverse the data sets on light meshes 
// in a descriptor.
typedef struct io_lite_mesh_iter io_lite_mesh_iter;

// This iterator allows one to traverse fields on a mesh.
typedef struct io_field_iter io_field_iter;

// Field centerings.
typedef enum
{
  IO_NODE,
  IO_EDGE,
  IO_FACE,
  IO_CELL
} io_field_centering_t;

// A function pointer type for opening file descriptors.
typedef int (*io_open_func)(void*, const char*, const char*, io_mode_t, int, MPI_Comm, int, int);

// A function pointer type for closing file descriptors.
typedef int (*io_close_func)(void*);

// Function pointer types for querying meshes from a descriptor.
typedef int (*io_query_num_meshes_func)(void*, int*);
typedef int (*io_query_mesh_names_func)(void*, char**);
typedef int (*io_query_num_lite_meshes_func)(void*, int*);
typedef int (*io_query_lite_meshe_names_func)(void*, char**);

// A function pointer type for reading a (heavy) mesh from an open descriptor.
typedef int (*io_read_mesh_func)(void*, const char*, mesh_t*);

// A function pointer type for reading a light mesh from an open descriptor.
typedef int (*io_read_lite_mesh_func)(void*, const char*, lite_mesh_t*);

// Function pointer types for querying data fields from a descriptor.
typedef int (*io_query_num_fields_func)(void*, const char*, int*);
typedef int (*io_query_field_names_func)(void*, const char*, char**);
typedef int (*io_query_field_func)(void*, const char*, const char*, int*, io_field_centering_t*);

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
  io_open_func                  open;
  io_close_func                 close;
  io_query_num_meshes_func      query_num_meshes;
  io_query_mesh_names_func      query_mesh_names;
  io_query_mesh_func            query_mesh;
  io_query_num_lite_meshes_func query_num_lite_meshes;
  io_query_lite_mesh_names_func query_lite_mesh_names;
  io_query_lite_mesh_func       query_lite_mesh;
  io_read_mesh_func             read_mesh;
  io_read_lite_mesh_func        read_lite_mesh;
  io_query_num_fields_func      query_field_names;
  io_query_field_names_func     query_field_names;
  io_query_field_func           query_field;
  io_read_field_func            read_field;
  io_write_mesh_func            write_mesh;
  io_write_lite_mesh_func       write_lite_mesh;
  io_write_field_func           write_field;
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

// Returns an iterator that can be used to traverse the heavy mesh 
// datasets available to this descriptor.
io_mesh_iter* io_meshes(io_interface_t* interface);

// Returns true if the given mesh iterator points to a mesh/dataset, 
// false if it does not.
bool io_mesh_iter_ok(io_mesh_iter* iter);

// Reads the heavy mesh from the descriptor, setting mesh_name to the name 
// of the mesh/dataset.
void io_mesh_iter_read(io_mesh_iter* iter, char** mesh_name, mesh_t* mesh);

// Writes a (heavy) mesh to an open descriptor via the given iterator.
void io_mesh_iter_write(io_mesh_iter* iter, const char* mesh_name, mesh_t* mesh);

// Increments the iterator.
void io_mesh_iter_next(io_mesh_iter* iter);

// Returns an iterator that can traverse the fields for the dataset.
io_field_iter* io_mesh_iter_fields(io_mesh_iter* iter);

// Returns an iterator that can be used to traverse the lite mesh
// datasets available to this descriptor.
io_lite_mesh_iter* io_lite_meshes(io_interface_t* interface);

// Returns true if the given mesh iterator points to a lite mesh/dataset, 
// false if it does not.
bool io_lite_mesh_iter_ok(io_lite_mesh_iter* iter);

// Reads the lite mesh from the descriptor, setting mesh_name to the name 
// of the mesh/dataset.
void io_lite_mesh_iter_read(io_lite_mesh_iter* interface, char** mesh_name, lite_mesh_t* lite_mesh);

// Writes a lite mesh to an open descriptor via the given iterator.
void io_lite_mesh_iter_write(io_lite_mesh_iter* iter, const char* mesh_name, mesh_t* mesh);

// Increments the iterator.
void io_lite_mesh_iter_next(io_lite_mesh_iter* iter);

// Returns an iterator that can traverse the fields for the dataset.
io_field_iter* io_lite_mesh_iter_fields(io_lite_mesh_iter* iter);

// Returns true if the given field iterator points to a field
// false if it does not.
bool io_field_iter_ok(io_field_iter* iter);

// Increments the iterator.
void io_field_iter_next(io_field_iter* iter);

// Gathers metadata about the current field. The size of the field is 
// consistent with its associated mesh.
void io_query_field(io_field_iter* iter, const char** field_name, int* num_components, io_field_centering_t* centering);

// Reads the field data from the iterator.
void io_field_iter_read(io_field_iter* iter, double* field_data);

// Writes a field to an open descriptor via the iterator.
void io_field_iter_write(io_field_iter* iter, const char* field_name, double* field_data, int num_data, int num_components);

#ifdef __cplusplus
}
#endif

#endif

