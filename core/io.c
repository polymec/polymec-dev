#include <stdlib.h>
#include <string.h>
#include <dirent.h>
#include <sys/stat.h>
#include "core/io.h"

#ifdef USE_MPI
#include <mpi.h>
#include "pmpio.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

// The "io" interface type is an opaque type for descriptors used to 
// read and write parallel data.
struct io_interface_t 
{
  // Type information.
  void* context;
  char* name;
  io_vtable vtable;

  // File stuff.
  io_mode_t mode;
  void* file;
  char filename[1024];
#ifdef USE_MPI
  PMPIO_baton_t* baton;
  char dir[1024];
#endif

  // Intermediate storage.
  io_buffered_data_t* buffered_data;
};

io_interface_t* io_interface_new(void* context, const char* name, io_vtable vtable)
{
  // Check the vtable.
  ASSERT(vtable.create_file != NULL);
  ASSERT(vtable.open_file != NULL);
  ASSERT(vtable.close_file != NULL);
  ASSERT(vtable.read_data != NULL);
  ASSERT(vtable.write_data != NULL);

  // Allocate the interface.
  io_interface_t* i = malloc(sizeof(io_interface_t));
  i->context = context;
  i->name = strdup(name);
  i->vtable = vtable;
  i->mode = IO_CLOSED;
  i->file = NULL;
  i->buffered_data = NULL;

#if USE_MPI
  i->baton = NULL;
#endif

  return i;
}

void io_free(io_interface_t* interface)
{
  ASSERT(interface->mode == IO_CLOSED);
  ASSERT(interface->file == NULL);
  ASSERT(interface->buffered_data == NULL);

  free(interface->name);
  if (interface->vtable.dtor != NULL)
    interface->vtable.dtor(interface->context);
  free(interface);
}

#if USE_MPI
static void* pmpio_create_file(const char* filename, 
                               const char* dirname,
                               void* userData)
{
  io_interface_t* io = (io_interface_t*)userData;
  return io->vtable.create_file(io->context, filename, dirname);
}

static void* pmpio_open_file(const char* filename, 
                             const char* dirname,
                             PMPIO_iomode_t io_mode, 
                             void* userData)
{
  io_interface_t* io = (io_interface_t*)userData;
  io_mode_t mode = (pmpio_mode == PMPIO_READ) ? IO_READ : IO_WRITE;
  io->file = io->vtable.open_file(io->context, filename, dirname, mode);
  io->mode = mode;
  return io->file;
}

static void* pmpio_close_file(void* file,
                              void* userData)
{
  io_interface_t* io = (io_interface_t*)userData;
  io->vtable.close_file(io->context, io->file);
  io->file = NULL;
  io->mode = IO_CLOSED;
}
#endif

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
  ASSERT(interface->buffered_data == NULL);

  // If we're reading a file, make sure the master directory exists.
  // If we're writing a file, create the directory if it doesn't exist.
#if USE_MPI
  if (rank == 0)
  {
#endif
    DIR* dir = opendir(directory);
    if (dir == NULL)
    {
      if (mode == IO_READ)
        arbi_error("io_open: directory %s does not exist.", directory);
      else 
        mkdir((char*)directory, S_IRWXU | S_IRWXG);
    }
    else
    {
      closedir(dir);
    }
#if USE_MPI
    MPI_Barrier(comm);
  }
  else
  {
    MPI_Barrier(comm);
  }

  int nproc = 1, rank = 0;
  MPI_Comm_size(comm, &nproc);
  MPI_Comm_rank(comm, &rank);
  if (num_files == -1)
    num_files = nproc;
  ASSERT(num_files <= nproc);

  // We put the entire data set into a directory named after the 
  // prefix, and every process gets its own subdirectory therein.

  // Initialize poor man's I/O and figure out group ranks.
  PMPIO_iomode_t pmpio_mode = (mode == IO_READ) ? PMPIO_READ : PMPIO_WRITE;
  silo->baton = PMPIO_Init(num_files, pmpio_mode, comm, mpi_tag, 
                           &pmpio_create_file,
                           &pmpio_open_file,
                           &pmpio_close_file, (void*)interface);
  int group_rank = PMPIO_group_rank(baton, rank);
  int rank_in_group = PMPIO_rank_in_group(baton, rank);

  // Create a subdirectory for each group.
  char group_dirname[1024];
  snprintf(group_dirname, 1024, "%s/%d", directory, group_rank);
  if ((rank_in_group == 0) && (mode == IO_WRITE))
  {
    DIR* group_dir = opendir(group_dirname);
    if (group_dir == 0)
      mkdir((char*)group_dirname, S_IRWXU | S_IRWXG);
    else
      closedir(group_dir);
    MPI_Barrier(comm);
  }
  else if (mode == IO_WRITE)
  {
    MPI_Barrier(comm);
  }

  // Determine a file name.
  snprintf(interface->filename, 1024, "%s/%s.silo", group_dirname, prefix);
  snprintf(interface->dir, 1024, "domain_%d", rank_in_group);
#else
  snprintf(interface->filename, 1024, "%s/%s.silo", directory, prefix);
#endif

  // If we're reading from the file, read the contents.
  if (mode == IO_READ)
  {
#if USE_MPI
    interface->file = PMPIO_WaitForBaton(interface->baton, interface->filename, interface->dir);
#else 
    interface->file = interface->vtable.create_file(interface->filename, "/", NULL);
#endif
    if (interface->file == NULL)
    {
      char err[1024];
      snprintf(err, 1024, "io_open: Could not open file descriptor for %s\n", prefix);
      arbi_error(err);
    }
    interface->buffered_data = malloc(sizeof(io_buffered_data_t));
    interface->vtable.read_data(interface->context, interface->file, interface->buffered_data);
  }
  interface->mode = mode;
}

void io_close(io_interface_t* interface)
{
  ASSERT(interface->mode != IO_CLOSED);

  // Flush any buffered data.
  if ((interface->mode == IO_WRITE) && (interface->buffered_data != NULL))
    interface->vtable.write_data(interface->context, interface->file, interface->buffered_data);

#if USE_MPI
  PMPIO_HandOffBaton(interface->baton, interface->file);
  PMPIO_Finish(interface->baton);
#else
  interface->vtable.close_file(interface->context, interface->file);
#endif
  interface->mode = IO_CLOSED;

  // Clean up any buffered data.
  if (interface->buffered_data != NULL)
  {
    free(interface->buffered_data);
    interface->buffered_data = NULL;
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

void io_dataset_read_mesh(io_dataset_t* dataset, mesh_t** mesh)
{
  ASSERT(dataset->interface->mode == IO_READ);
#if 0
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
#endif
}

void io_dataset_write_mesh(io_dataset_t* dataset, mesh_t* mesh)
{
  ASSERT(dataset->interface->mode == IO_WRITE);
  ASSERT(mesh != NULL);
#if 0
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
#endif
}

void io_dataset_read_lite_mesh(io_dataset_t* dataset, lite_mesh_t** mesh)
{
  ASSERT(dataset->interface->mode == IO_READ);
#if 0
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
#endif
}

void io_dataset_write_lite_mesh(io_dataset_t* dataset, lite_mesh_t* mesh)
{
  ASSERT(dataset->interface->mode == IO_WRITE);
  ASSERT(mesh != NULL);
#if 0
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
#endif
}

void io_dataset_query_field(io_dataset_t* dataset, const char* field_name, int* num_components, mesh_centering_t* centering)
{
  ASSERT(dataset->interface->mode == IO_READ);
#if 0
  if (dataset->interface->vtable.query_field(dataset->interface->context, dataset->name, field_name, num_components, centering) != ARBI_SUCCESS)
  {
    char err[1024];
    snprintf(err, 1024, "Could not query field '%s' from dataset '%s'.\n", field_name, dataset->name);
    arbi_error(err);
  }
#endif
}

void io_dataset_read_field(io_dataset_t* dataset, const char* field_name, double* field_data)
{
  ASSERT(dataset->interface->mode == IO_READ);
#if 0
  if (dataset->interface->vtable.read_field(dataset->interface->context, dataset->name, field_name, field_data) != ARBI_SUCCESS)
  {
    char err[1024];
    snprintf(err, 1024, "Could not read field '%s' from dataset '%s'.\n", field_name, dataset->name);
    arbi_error(err);
  }
#endif
}

void io_dataset_write_field(io_dataset_t* dataset, const char* field_name, double* field_data, int num_components, mesh_centering_t centering)
{
  ASSERT(dataset->interface->mode == IO_WRITE);
  ASSERT(field_data != NULL);
  ASSERT(num_data > 0);
  ASSERT(num_components > 0);
#if 0
  if (dataset->interface->vtable.write_field(dataset->interface->context, dataset->name, field_name, field_data, num_components, centering) != ARBI_SUCCESS)
  {
    char err[1024];
    snprintf(err, 1024, "Could not write field '%s' to dataset '%s'.\n", field_name, dataset->name);
    arbi_error(err);
  }
#endif
}

void io_dataset_query_source_code(io_dataset_t* dataset, const char* code_name, int* len)
{
  ASSERT(dataset->interface->mode == IO_READ);
#if 0
  if (dataset->interface->vtable.query_source_code(dataset->interface->context, dataset->name, code_name, len) != ARBI_SUCCESS)
  {
    char err[1024];
    snprintf(err, 1024, "Could not query source code '%s' from dataset '%s'.\n", code_name, dataset->name);
    arbi_error(err);
  }
#endif
}

void io_dataset_read_source_code(io_dataset_t* dataset, const char* code_name, char* source_code)
{
  ASSERT(dataset->interface->mode == IO_READ);
#if 0
  if (dataset->interface->vtable.read_source_code(dataset->interface->context, dataset->name, code_name, source_code) != ARBI_SUCCESS)
  {
    char err[1024];
    snprintf(err, 1024, "Could not read source code '%s' from dataset '%s'.\n", code_name, dataset->name);
    arbi_error(err);
  }
#endif
}

void io_dataset_write_source_code(io_dataset_t* dataset, const char* code_name, const char* source_code)
{
  ASSERT(dataset->interface->mode == IO_WRITE);
#if 0
  if (dataset->interface->vtable.write_source_code(dataset->interface->context, dataset->name, code_name, source_code) != ARBI_SUCCESS)
  {
    char err[1024];
    snprintf(err, 1024, "Could not write source code '%s' to dataset '%s'.\n", code_name, dataset->name);
    arbi_error(err);
  }
#endif
}

#ifdef __cplusplus
}
#endif

