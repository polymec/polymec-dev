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
  char directory[1024];
  char prefix[1024];
  char filename[1024];
#ifdef USE_MPI
  PMPIO_baton_t* baton;
  char group_dir[1024];
  MPI_Comm comm;
  int rank, nproc;
  int num_files;
  int mpi_tag;
#endif

  // Intermediate storage.
  io_dataset_t** datasets;
  int num_datasets;
};

io_interface_t* io_interface_new(void* context, 
                                 const char* name, 
                                 io_vtable vtable,
                                 MPI_Comm comm, 
                                 int num_files,
                                 int mpi_tag)
{
  // Check the vtable.
  ASSERT(vtable.create_file != NULL);
  ASSERT(vtable.open_file != NULL);
  ASSERT(vtable.close_file != NULL);
  ASSERT((vtable.read_datasets != NULL) || (vtable.write_datasets != NULL));

  // Allocate the interface.
  io_interface_t* i = malloc(sizeof(io_interface_t));
  i->context = context;
  i->name = strdup(name);
  i->vtable = vtable;
  i->mode = IO_CLOSED;
  i->file = NULL;

#if USE_MPI
  i->baton = NULL;
  i->comm = comm;
  i->num_files = num_files;
  i->mpi_tag = mpi_tag;
#endif

  // Dataset storage.
  i->datasets = NULL;
  i->num_datasets = 0;

  return i;
}

io_interface_t* io_interface_new_serial(void* context, 
                                        const char* name, 
                                        io_vtable vtable)
{
  return io_interface_new(context, name, vtable, MPI_COMM_WORLD, 1, 0);
}

void io_free(io_interface_t* interface)
{
  ASSERT(interface->mode == IO_CLOSED);
  ASSERT(interface->file == NULL);

  free(interface->name);
  if (interface->vtable.dtor != NULL)
    interface->vtable.dtor(interface->context);

  for (int i = 0; i < interface->num_datasets; ++i)
    io_dataset_free(interface->datasets[i]);
  free(interface->datasets);

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
             io_mode_t mode)
{
  ASSERT(interface->mode == IO_CLOSED);
  ASSERT(mode != IO_CLOSED);

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
    MPI_Barrier(interface->comm);
  }
  else
  {
    MPI_Barrier(interface->comm);
  }

  interface->nproc = 1, interface->rank = 0;
  MPI_Comm_size(interface->comm, &interface->nproc);
  MPI_Comm_rank(interface->comm, &interface->rank);
  if (interface->num_files == -1)
    interface->num_files = interface->nproc;
  ASSERT(interface->num_files <= interface->nproc);

  // We put the entire data set into a directory named after the 
  // prefix, and every process gets its own subdirectory therein.

  // Initialize poor man's I/O and figure out group ranks.
  PMPIO_iomode_t pmpio_mode = (mode == IO_READ) ? PMPIO_READ : PMPIO_WRITE;
  interface->baton = PMPIO_Init(interface->num_files, pmpio_mode, interface->comm, interface->mpi_tag, 
                           &pmpio_create_file,
                           &pmpio_open_file,
                           &pmpio_close_file, (void*)interface);
  int group_rank = PMPIO_group_rank(interface->baton, interface->rank);
  int rank_in_group = PMPIO_rank_in_group(interface->baton, interface->rank);

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
    MPI_Barrier(interface->comm);
  }
  else if (mode == IO_WRITE)
  {
    MPI_Barrier(interface->comm);
  }

  // Determine a file name.
  snprintf(interface->filename, 1024, "%s/%s.silo", group_dirname, prefix);
  snprintf(interface->group_dir, 1024, "domain_%d", rank_in_group);
#else
  snprintf(interface->filename, 1024, "%s/%s.silo", directory, prefix);
#endif
  strncpy(interface->prefix, prefix, 1024);
  strncpy(interface->directory, directory, 1024);

  // If we're reading from the file, read the contents.
  interface->mode = mode;
  if (mode == IO_READ)
  {
#if USE_MPI
    interface->file = PMPIO_WaitForBaton(interface->baton, interface->filename, interface->group_dir);
#else 
    interface->file = interface->vtable.create_file(interface->context, interface->filename, "/");
#endif
    if (interface->file == NULL)
    {
      char err[1024];
      snprintf(err, 1024, "io_open: Could not open file descriptor for %s\n", prefix);
      arbi_error(err);
    }
    if (interface->vtable.get_num_datasets != NULL)
      interface->vtable.get_num_datasets(interface->context, interface->file, &interface->num_datasets);
    else
      interface->num_datasets = 1;
    interface->datasets = malloc(interface->num_datasets*sizeof(io_dataset_t*));
    memset(interface->datasets, 0, interface->num_datasets*sizeof(io_dataset_t*));
    interface->vtable.read_datasets(interface->context, interface->file, interface->datasets, interface->num_datasets);
  }
  // Otherwise, set some defaults for writing.
  else if (mode == IO_WRITE)
  {
    io_set_num_datasets(interface, 1); // 1 dataset by default.
  }
}

void io_close(io_interface_t* interface)
{
  ASSERT(interface->mode != IO_CLOSED);

  // Flush any buffered data.
  if ((interface->mode == IO_WRITE) && (interface->datasets != NULL))
  {
#if USE_MPI
    int rank_in_group = PMPIO_rank_in_group(interface->baton, interface->rank);
    int num_procs_per_file = interface->nproc / num_files;
#else
    int rank_in_group = 0;
    int num_procs_per_file = 1;
#endif
    interface->vtable.write_datasets(interface->context, interface->file, interface->datasets, interface->num_datasets, rank_in_group, num_procs_per_file);

#if USE_MPI
    // Write a master file if needed.
    if (interface->rank == 0)
    {
      char master_filename[1024];
      snprintf(master_filename, 1024, "%s/%s.silo", interface->directory, interface->prefix);
      void* master = interface->vtable.create_file(interface->context, master_filename, "/");
      interface->vtable.write_master(interface->context, master, interface->prefix, interface->datasets, interface->num_datasets, interface->num_files, num_procs_per_file);
      interface->vtable.close_file(interface->context, master);
    }
#endif
  }

#if USE_MPI
  PMPIO_HandOffBaton(interface->baton, interface->file);
  PMPIO_Finish(interface->baton);
#else
  interface->vtable.close_file(interface->context, interface->file);
#endif
  interface->mode = IO_CLOSED;
}

int io_num_datasets(io_interface_t* interface)
{
  return interface->num_datasets;
}

const char* io_dataset_name(io_interface_t* interface, int index)
{
  ASSERT(index >= 0);
  ASSERT(index < io_num_datasets(interface));
  return interface->datasets[index]->name;
}

void io_set_num_datasets(io_interface_t* interface, int num_datasets)
{
  ASSERT(interface->mode == IO_WRITE);
  ASSERT(num_datasets >= 1);
  ASSERT(interface->datasets == NULL);
  interface->num_datasets = num_datasets;
  interface->datasets = malloc(sizeof(io_dataset_t*)*num_datasets);
  memset(interface->datasets, 0, num_datasets*sizeof(io_dataset_t*));
}

io_dataset_t* io_dataset(io_interface_t* interface, const char* dataset)
{
  // Find the dataset with this name within the interface.
  for (int i = 0; i < interface->num_datasets; ++i)
  {
    if (!strcmp(interface->datasets[i]->name, dataset))
      return interface->datasets[i];
  }
  return NULL;
}

io_dataset_t* io_default_dataset(io_interface_t* interface)
{
  if (interface->num_datasets == 0)
    return NULL;
  else
    return interface->datasets[0];
}

io_dataset_t* io_dataset_new(io_interface_t* interface, const char* name,
                             int num_fields, int num_sources)
{
  ASSERT(num_fields >= 0);
  ASSERT(num_sources >= 0);
  int index = 0;
  while ((index < interface->num_datasets) && (interface->datasets[index] != NULL))
    ++index;
  if (index == interface->num_datasets) // No room!
    return NULL;

  io_dataset_t* d = malloc(sizeof(io_dataset_t));
  d->interface = interface;
  d->name = strdup(name);
  d->mesh = NULL;
  d->lite_mesh = NULL;
  d->fields = malloc(num_fields*sizeof(double*));
  d->field_names = malloc(num_fields*sizeof(char*));
  d->num_fields = 0;
  d->sources = malloc(num_sources*sizeof(char*));
  d->source_names = malloc(num_sources*sizeof(char*));
  d->num_sources = 0;

  // Place the dataset in its proper place within the interface.
  interface->datasets[index] = d;

  return d;
}

void io_dataset_free(io_dataset_t* dataset)
{
  free(dataset->name);
  if (dataset->mesh != NULL)
    mesh_free(dataset->mesh);
  if (dataset->lite_mesh != NULL)
    lite_mesh_free(dataset->lite_mesh);

  for (int i = 0; i < dataset->num_fields; ++i)
  {
    if (dataset->fields[i] != NULL)
      free(dataset->fields[i]);
    free(dataset->field_names[i]);
  }
  free(dataset->fields);
  free(dataset->field_names);

  for (int i = 0; i < dataset->num_sources; ++i)
  {
    if (dataset->sources[i] != NULL)
      free(dataset->sources[i]);
    free(dataset->source_names[i]);
  }
  free(dataset->sources);
  free(dataset->source_names);

  free(dataset);
}

void io_dataset_read_mesh(io_dataset_t* dataset, mesh_t** mesh)
{
  ASSERT(dataset->interface->mode == IO_READ);
  *mesh = dataset->mesh;
  dataset->mesh = NULL;
}

void io_dataset_write_mesh(io_dataset_t* dataset, mesh_t* mesh)
{
  ASSERT(dataset->interface->mode == IO_WRITE);
  ASSERT(mesh != NULL);
  dataset->mesh = mesh;
}

void io_dataset_read_lite_mesh(io_dataset_t* dataset, lite_mesh_t** mesh)
{
  ASSERT(dataset->interface->mode == IO_READ);
  *mesh = dataset->lite_mesh;
  dataset->lite_mesh = NULL;
}

void io_dataset_write_lite_mesh(io_dataset_t* dataset, lite_mesh_t* mesh)
{
  ASSERT(dataset->interface->mode == IO_WRITE);
  ASSERT(mesh != NULL);
  dataset->lite_mesh = mesh;
}

void io_dataset_query_field(io_dataset_t* dataset, const char* field_name, int* num_components, mesh_centering_t* centering)
{
  ASSERT(dataset->interface->mode == IO_READ);
  *num_components = -1;
  for (int i = 0; i < dataset->num_fields; ++i)
  {
    if (!strcmp(dataset->field_names[i], field_name))
    {
      *num_components = dataset->field_num_comps[i];
      *centering = dataset->field_centerings[i];
    }
  }
}

void io_dataset_read_field(io_dataset_t* dataset, const char* field_name, double** field)
{
  ASSERT(dataset->interface->mode == IO_READ);
  *field = NULL;
  for (int i = 0; i < dataset->num_fields; ++i)
  {
    if (!strcmp(dataset->field_names[i], field_name))
    {
      *field = dataset->fields[i];
      dataset->fields[i] = NULL;
    }
  }
}

void io_dataset_write_field(io_dataset_t* dataset, const char* field_name, double* field_data, int num_components, mesh_centering_t centering)
{
  ASSERT(dataset->interface->mode == IO_WRITE);
  ASSERT(field_data != NULL);
  ASSERT(num_components > 0);
  for (int i = 0; i < dataset->num_fields; ++i)
  {
    if (!strcmp(dataset->field_names[i], field_name))
    {
      dataset->fields[i] = field_data;
      dataset->field_num_comps[i] = num_components;
      dataset->field_centerings[i] = centering;
    }
  }
}

void io_dataset_query_source_code(io_dataset_t* dataset, const char* code_name, int* len)
{
  ASSERT(dataset->interface->mode == IO_READ);
  *len = -1;
  for (int i = 0; i < dataset->num_sources; ++i)
  {
    if (!strcmp(dataset->source_names[i], code_name))
    {
      *len = dataset->source_lengths[i];
    }
  }
}

void io_dataset_read_source_code(io_dataset_t* dataset, const char* code_name, char** source_code)
{
  ASSERT(dataset->interface->mode == IO_READ);
  for (int i = 0; i < dataset->num_sources; ++i)
  {
    if (!strcmp(dataset->source_names[i], code_name))
    {
      *source_code = dataset->sources[i];
      dataset->sources[i] = NULL;
    }
  }
}

void io_dataset_write_source_code(io_dataset_t* dataset, const char* code_name, const char* source_code)
{
  ASSERT(dataset->interface->mode == IO_WRITE);
  ASSERT(source_code != NULL);
  for (int i = 0; i < dataset->num_fields; ++i)
  {
    if (!strcmp(dataset->source_names[i], code_name))
    {
      dataset->sources[i] = (char*)source_code;
      dataset->source_lengths[i] = strlen(source_code);
    }
  }
}

#ifdef __cplusplus
}
#endif

