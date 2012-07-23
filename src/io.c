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
             MPI_Comm comm,
             int num_files,
             int mpi_tag)
{
  ASSERT(interface->mode == IO_CLOSED);
  ASSERT(mode != IO_CLOSED);
  if (interface->vtable.open(interface->context, prefix, directory, mode, comm, num_files, mpi_tag) != ARBI_SUCCESS)
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

struct io_dataset_t
{
  io_interface_t* interface;
  char* name;
};

io_dataset_t* io_dataset(io_interface_t* interface, const char* dataset)
{
  io_dataset_t* d = malloc(sizeof(io_dataset_t));
  d->interface = interface;
  d->name = strdup(dataset);
  return d;
}

void io_dataset_free(io_dataset_t* dataset)
{
  free(dataset->name);
  free(dataset);
}

void io_dataset_read_mesh(io_dataset_t* dataset, mesh_t* mesh)
{
  ASSERT(dataset->interface->mode == IO_READ);
  if (dataset->interface->vtable.read_mesh == NULL)
  {
    char err[1024];
    snprintf(err, 1024, "The '%s' I/O interface does not support reading heavy meshes.\n", dataset->interface->name);
    arbi_error(err);
  }
  if (dataset->interface->vtable.read_mesh(dataset->interface->context, dataset->name, mesh) != ARBI_SUCCESS)
  {
    char err[1024];
    snprintf(err, 1024, "Could not read mesh from dataset '%s'.\n", dataset->name);
    arbi_error(err);
  }
}

void io_dataset_write_mesh(io_dataset_t* dataset, mesh_t* mesh)
{
  ASSERT(dataset->interface->mode == IO_WRITE);
  ASSERT(mesh != NULL);
  if (dataset->interface->vtable.write_mesh == NULL)
  {
    char err[1024];
    snprintf(err, 1024, "The '%s' I/O interface does not support writing heavy meshes.\n", dataset->interface->name);
    arbi_error(err);
  }
  if (dataset->interface->vtable.write_mesh(dataset->interface->context, dataset->name, mesh) != ARBI_SUCCESS)
  {
    char err[1024];
    snprintf(err, 1024, "Could not write mesh to dataset '%s'.\n", dataset->name);
    arbi_error(err);
  }
}

void io_dataset_read_lite_mesh(io_dataset_t* dataset, lite_mesh_t* mesh)
{
  ASSERT(dataset->interface->mode == IO_READ);
  if (dataset->interface->vtable.read_lite_mesh == NULL)
  {
    char err[1024];
    snprintf(err, 1024, "The '%s' I/O interface does not support reading lite meshes.\n", dataset->interface->name);
    arbi_error(err);
  }
  if (dataset->interface->vtable.read_lite_mesh(dataset->interface->context, dataset->name, mesh) != ARBI_SUCCESS)
  {
    char err[1024];
    snprintf(err, 1024, "Could not read lite mesh from dataset '%s'.\n", dataset->name);
    arbi_error(err);
  }
}

void io_dataset_write_lite_mesh(io_dataset_t* dataset, lite_mesh_t* mesh)
{
  ASSERT(dataset->interface->mode == IO_WRITE);
  ASSERT(mesh != NULL);
  if (dataset->interface->vtable.write_lite_mesh == NULL)
  {
    char err[1024];
    snprintf(err, 1024, "The '%s' I/O interface does not support writing lite meshes.\n", dataset->interface->name);
    arbi_error(err);
  }
  if (dataset->interface->vtable.write_lite_mesh(dataset->interface->context, dataset->name, mesh) != ARBI_SUCCESS)
  {
    char err[1024];
    snprintf(err, 1024, "Could not write lite mesh to dataset '%s'.\n", dataset->name);
    arbi_error(err);
  }
}

void io_dataset_query_field(io_dataset_t* dataset, const char* field_name, int* num_components, io_field_centering_t* centering)
{
  ASSERT(dataset->interface->mode == IO_READ);
  if (dataset->interface->vtable.query_field(dataset->interface->context, dataset->name, field_name, num_components, centering) != ARBI_SUCCESS)
  {
    char err[1024];
    snprintf(err, 1024, "Could not query field '%s' from dataset '%s'.\n", field_name, dataset->name);
    arbi_error(err);
  }
}

void io_dataset_read_field(io_dataset_t* dataset, const char* field_name, double* field_data)
{
  ASSERT(dataset->interface->mode == IO_READ);
  if (dataset->interface->vtable.read_field(dataset->interface->context, field_name, field_data) != ARBI_SUCCESS)
  {
    char err[1024];
    snprintf(err, 1024, "Could not read field '%s' from dataset '%s'.\n", field_name, dataset->name);
    arbi_error(err);
  }
}

void io_dataset_write_field(io_dataset_t* dataset, const char* field_name, double* field_data, int num_components, io_field_centering_t centering)
{
  ASSERT(dataset->interface->mode == IO_WRITE);
  ASSERT(field_data != NULL);
  ASSERT(num_data > 0);
  ASSERT(num_components > 0);
  if (dataset->interface->vtable.write_field(dataset->interface->context, dataset->name, field_name, field_data, num_components, centering) != ARBI_SUCCESS)
  {
    char err[1024];
    snprintf(err, 1024, "Could not write field '%s' to dataset '%s'.\n", field_name, dataset->name);
    arbi_error(err);
  }
}

#ifdef __cplusplus
}
#endif

