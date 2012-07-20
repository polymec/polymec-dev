#include <stdlib.h>
#include <string.h>
#include "io.h"

#ifdef __cplusplus
extern "C" {
#endif

// The "io" interface type is an opaque type for descriptors used to 
// read and write parallel data.
struct io_interface_t 
{
  void* context;
  char* name;
  io_mode_t mode;
  io_vtable vtable;
};

io_interface_t* io_interface_new(void* context, const char* name, io_vtable vtable)
{
  // Check the vtable.
  ASSERT(vtable.open != NULL);
  ASSERT(vtable.close != NULL);
  ASSERT((vtable.read_mesh != NULL) || (vtable.read_lite_mesh != NULL));
  ASSERT(((vtable.read_mesh != NULL) && (vtable.query_meshes != NULL)) ||
         ((vtable.read_lite_mesh != NULL) && (vtable.query_lite_meshes != NULL)));
  ASSERT(vtable.query_fields != NULL);
  ASSERT(vtable.query_field != NULL);
  ASSERT(vtable.read_field != NULL);
  ASSERT((vtable.write_mesh != NULL) || (vtable.write_lite_mesh != NULL));
  ASSERT(((vtable.read_mesh != NULL) && (vtable.write_mesh != NULL)) ||
         ((vtable.read_lite_mesh != NULL) && (vtable.write_lite_mesh != NULL)));
  ASSERT(vtable.write_field != NULL);

  // Allocate the interface.
  io_interface_t* i = malloc(sizeof(io_interface_t));
  i->context = context;
  i->name = strdup(name);
  i->mode = IO_CLOSED;
  i->vtable = vtable;
  return i;
}

void io_free(io_interface_t* interface)
{
  free(interface->name);
  if (interface->vtable.dtor != NULL)
    interface->vtable.dtor(interface->context);
  free(interface);
}

void io_open(io_interface_t* interface, 
             const char* prefix, 
             const char* directory,  
             io_mode_t mode,
             int cycle, 
             MPI_Comm comm,
             int num_files,
             int mpi_tag)
{
  ASSERT(interface->mode == IO_CLOSED);
  ASSERT(mode != IO_CLOSED);
  if (interface->vtable.open(interface->context, prefix, directory, mode, cycle, comm, num_files, mpi_tag) != ARBI_SUCCESS)
  {
    char err[1024];
    snprintf(err, 1024, "Could not open file descriptor for %s\n", prefix);
    arbi_error(err);
  }
  else
  {
    interface->mode = mode;
  }
}

void io_close(io_interface_t* interface)
{
  ASSERT(interface->mode != IO_CLOSED);
  if (interface->vtable.close(interface->context) != ARBI_SUCCESS)
  {
    arbi_error("Could not close file descriptor.\n");
  }
  else
  {
    interface->mode = IO_CLOSED;
  }
}

void io_query_meshes(io_interface_t* interface, char** mesh_names, int* num_meshes)
{
  ASSERT(interface->mode == IO_READ);
  if (interface->vtable.read_mesh == NULL)
  {
    char err[1024];
    snprintf(err, 1024, "The '%s' I/O interface does not support reading heavy meshes.\n", interface->name);
    arbi_error(err);
  }
  if (interface->vtable.query_meshes(interface->context, mesh_names, num_meshes) != ARBI_SUCCESS)
    arbi_error("Could not query meshes from file descriptor.\n");
}

void io_read_mesh(io_interface_t* interface, const char* mesh_name, mesh_t* mesh)
{
  ASSERT(interface->mode == IO_READ);
  if (interface->vtable.read_mesh == NULL)
  {
    char err[1024];
    snprintf(err, 1024, "The '%s' I/O interface does not support reading heavy meshes.\n", interface->name);
    arbi_error(err);
  }
  if (interface->vtable.read_mesh(interface->context, mesh_name, mesh) != ARBI_SUCCESS)
  {
    char err[1024];
    snprintf(err, 1024, "Could not read mesh '%s' from file descriptor.\n", mesh_name);
    arbi_error(err);
  }
}

void io_query_lite_meshes(io_interface_t* interface, char** lite_mesh_names, int* num_lite_meshes)
{
  ASSERT(interface->mode == IO_READ);
  if (interface->vtable.read_lite_mesh == NULL)
  {
    char err[1024];
    snprintf(err, 1024, "The '%s' I/O interface does not support reading lite meshes.\n", interface->name);
    arbi_error(err);
  }
  if (interface->vtable.query_lite_meshes(interface->context, lite_mesh_names, num_lite_meshes) != ARBI_SUCCESS)
    arbi_error("Could not query lite meshes from file descriptor.\n");
}

void io_read_lite_mesh(io_interface_t* interface, const char* lite_mesh_name, lite_mesh_t* lite_mesh)
{
  ASSERT(interface->mode == IO_READ);
  if (interface->vtable.read_lite_mesh == NULL)
  {
    char err[1024];
    snprintf(err, 1024, "The '%s' I/O interface does not support reading lite meshes.\n", interface->name);
    arbi_error(err);
  }
  if (interface->vtable.read_lite_mesh(interface->context, lite_mesh_name, lite_mesh) != ARBI_SUCCESS)
  {
    char err[1024];
    snprintf(err, 1024, "Could not read lite mesh '%s' from file descriptor.\n", lite_mesh_name);
    arbi_error(err);
  }
}

void io_query_fields(io_interface_t* interface, char** field_names, int* num_fields)
{
  ASSERT(interface->mode == IO_READ);
  if (interface->vtable.query_fields(interface->context, field_names, num_fields) != ARBI_SUCCESS)
    arbi_error("Could not query fields from file descriptor.\n");
}

void io_query_field(io_interface_t* interface, const char* field_name, int* num_data, int* num_components)
{
  ASSERT(interface->mode == IO_READ);
  if (interface->vtable.query_field(interface->context, field_name, num_data, num_components) != ARBI_SUCCESS)
  {
    char err[1024];
    snprintf(err, 1024, "Could not query field '%s' from file descriptor.\n", field_name);
    arbi_error(err);
  }
}

void io_read_field(io_interface_t* interface, const char* field_name, double* field_data)
{
  ASSERT(interface->mode == IO_READ);
  if (interface->vtable.read_field(interface->context, field_name, field_data) != ARBI_SUCCESS)
  {
    char err[1024];
    snprintf(err, 1024, "Could not read field '%s' from file descriptor.\n", field_name);
    arbi_error(err);
  }
}

void io_write_mesh(io_interface_t* interface, const char* mesh_name, mesh_t* mesh)
{
  ASSERT(interface->mode == IO_WRITE);
  ASSERT(mesh != NULL);
  if (interface->vtable.write_mesh == NULL)
  {
    char err[1024];
    snprintf(err, 1024, "The '%s' I/O interface does not support writing heavy meshes.\n", interface->name);
    arbi_error(err);
  }
  if (interface->vtable.write_mesh(interface->context, mesh_name, mesh) != ARBI_SUCCESS)
  {
    char err[1024];
    snprintf(err, 1024, "Could not write mesh '%s' to file descriptor.\n", mesh_name);
    arbi_error(err);
  }
}

void io_write_lite_mesh(io_interface_t* interface, const char* lite_mesh_name, lite_mesh_t* lite_mesh)
{
  ASSERT(interface->mode == IO_WRITE);
  ASSERT(lite_mesh != NULL);
  if (interface->vtable.write_lite_mesh == NULL)
  {
    char err[1024];
    snprintf(err, 1024, "The '%s' I/O interface does not support writing lite meshes.\n", interface->name);
    arbi_error(err);
  }
  if (interface->vtable.write_lite_mesh(interface->context, lite_mesh_name, lite_mesh) != ARBI_SUCCESS)
  {
    char err[1024];
    snprintf(err, 1024, "Could not write lite mesh '%s' to file descriptor.\n", lite_mesh_name);
    arbi_error(err);
  }
}

void io_write_field(io_interface_t* interface, const char* field_name, double* field_data, int num_data, int num_components)
{
  ASSERT(interface->mode == IO_WRITE);
  ASSERT(field_data != NULL);
  ASSERT(num_data > 0);
  ASSERT(num_components > 0);
  if (interface->vtable.write_field(interface->context, field_name, field_data, num_data, num_components) != ARBI_SUCCESS)
  {
    char err[1024];
    snprintf(err, 1024, "Could not write field '%s' to file descriptor.\n", field_name);
    arbi_error(err);
  }
}

#ifdef __cplusplus
}
#endif

