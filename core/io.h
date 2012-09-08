#ifndef ARBI_IO_H
#define ARBI_IO_H

#include "core/arbi.h"
#include "core/mesh.h"
#include "core/lite_mesh.h"
#include "tommytrie.h"

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

// A struct holding buffered data that has been read from or is to 
// be dumped to a file descriptor.
typedef struct
{
  tommy_trie* meshes;     
  tommy_trie* light_meshes; 
  tommy_trie* fields;       
  tommy_trie* source_code;       
} io_buffered_data_t;

// A datum in a tommy_trie storing an I/O dataset to be read or written 
// by the I/O interface.
typedef struct 
{
  tommy_node node;
  io_buffered_data_t* data;
} io_trie_dataset_t;

// This descriptor allows one to retrieve data for a given dataset.
typedef struct io_dataset_t io_dataset_t;

// A function pointer type for creating file descriptors.
typedef void* (*io_create_file_func)(void*, const char*, const char*);

// A function pointer type for opening file descriptors.
typedef void* (*io_open_file_func)(void*, const char* , const char*, io_mode_t);

// A function pointer type for closing file descriptors.
typedef void (*io_close_file_func)(void*, void*);

// A function pointer for reading data from all datasets in a file.
typedef void (*io_read_data_func)(void*, void*, io_buffered_data_t* data);

// A function pointer for dumping data to datasets in a file.
typedef void (*io_write_data_func)(void*, void*, io_buffered_data_t* data);

// A destructor function for the context object (if any).
typedef void (*io_dtor)(void*);

// This virtual table must be implemented by any I/O interface used 
// within arbi.
typedef struct 
{
  io_create_file_func           create_file;
  io_open_file_func             open_file;
  io_close_file_func            close_file;
  io_read_data_func             read_data;
  io_write_data_func            write_data;
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

// Returns the default dataset descriptor for the file. Use this to retrieve
// data from files that can hold only one dataset.
io_dataset_t* io_default_dataset(io_interface_t* interface);

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

