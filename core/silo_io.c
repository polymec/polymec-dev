#include <stdlib.h>
#include <sys/stat.h>
#include <dirent.h>
#include "silo.h"
#include "core/silo_io.h"

#ifdef USE_MPI
#include <mpi.h>
#include "pmpio.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

#ifdef USE_MPI

void* pmpio_create_file(const char* filename,
                        const char* dirname,
                        void* userData)
{
  int driver = DB_HDF5;
  DBfile* file = DBCreate(filename, 0, DB_LOCAL, 0, driver);
  DBMkDir(file, dirname);
  DBSetDir(file, dirname);
  return (void*)file;
}

void* pmpio_open_file(const char* filename, 
                      const char* dirname,
                      PMPIO_iomode_t iomode, 
                      void* userData)
{
  int driver = DB_HDF5;
  DBfile* file;
  if (iomode == PMPIO_WRITE)
  { 
    file = DBCreate(filename, 0, DB_LOCAL, 0, driver);
    DBMkDir(file, dirname);
    DBSetDir(file, dirname);
  }
  else
  {
    file = DBOpen(filename, driver, DB_READ);
    DBSetDir(file, dirname);
  }
  return (void*)file;
}

void pmpio_close_file(void* file, void* userData)
{
  DBClose((DBfile*)file);
}

#endif

typedef struct
{
  char filename[1024];
  char dir[1024];
#ifdef USE_MPI
  PMPIO_baton_t* baton;
#endif
  DBfile* file;
} silo_context;

static int silo_open(void* context, const char* prefix, const char* directory, io_mode_t mode, MPI_Comm comm, int num_files, int mpi_tag)
{
  silo_context* silo = (silo_context*)context;

  // Open a file in Silo/HDF5 format.
  char filename[1024];
#ifdef USE_MPI
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
                           &pmpio_createFile, &pmpio_openFile, 
                           &pmpio_closeFile, 0);
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
  snprintf(silo->filename, 1024, "%s/%s.silo", group_dirname, prefix);
  snprintf(silo->dir, 1024, "domain_%d", rank_in_group);
#else

  snprintf(silo->filename, 1024, "%s/%s.silo", directory, prefix);

  int driver = DB_HDF5;
  silo->file = DBCreate(silo->filename, 0, DB_LOCAL, 0, driver);
  DBSetDir(silo->file, "/");
#endif

  return 0;
}

static int silo_close(void* context)
{
  //silo->file = (DBfile*)PMPIO_WaitForBaton(silo->baton, silo->filename, silo->dir);
  //PMPIO_HandOffBaton(silo->baton, (void*)silo->file); // FIXME
  //PMPIO_Finish(silo->baton);
  //DBClose(file);
  return 0;
}

static int silo_read_mesh(void* context, const char* dataset, mesh_t** mesh)
{
  return 0;
}

static int silo_write_mesh(void* context, const char* dataset, mesh_t* mesh)
{
  return 0;
}

static int silo_query_field(void* context, const char* dataset, const char* field_name, int* num_components, mesh_centering_t* centering)
{
  return 0;
}

static int silo_read_field(void* context, const char* dataset, const char* field_name, double* data)
{
  return 0;
}

static int silo_write_field(void* context, const char* dataset, const char* field_name, double* data, int num_components, mesh_centering_t centering)
{
  return 0;
}

static int silo_query_source_code(void* context, const char* dataset, const char* code_name, int* length)
{
  return 0;
}

static int silo_read_source_code(void* context, const char* dataset, const char* code_name, char* source_code)
{
  return 0;
}

static int silo_write_source_code(void* context, const char* dataset, const char* code_name, const char* source_code)
{
  return 0;
}

static void silo_dtor(void* context)
{
  free(context);
}

io_interface_t* silo_io_new()
{
  silo_context* context = malloc(sizeof(silo_context));
  io_vtable vtable = {.open = &silo_open, 
                      .close = &silo_close,
                      .read_mesh = &silo_read_mesh,
                      .query_field = &silo_query_field,
                      .read_field = &silo_read_field,
                      .query_source_code = &silo_query_source_code,
                      .read_source_code = &silo_read_source_code,
                      .write_mesh = &silo_write_mesh,
                      .write_field = &silo_write_field,
                      .write_source_code = &silo_write_source_code,
                      .dtor = &silo_dtor};
  return io_interface_new(context, "SILO", vtable);
}

#ifdef __cplusplus
}
#endif

