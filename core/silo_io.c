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

typedef struct
{
} silo_context;

static void* silo_create_file(const char* filename,
                              const char* dirname,
                              void* userData)
{
  int driver = DB_HDF5;
  DBfile* file = DBCreate(filename, 0, DB_LOCAL, 0, driver);
  DBMkDir(file, dirname);
  DBSetDir(file, dirname);
  return (void*)file;
}

static void* silo_open_file(const char* filename, 
                            const char* dirname,
                            io_mode_t mode, 
                            void* userData)
{
  int driver = DB_HDF5;
  DBfile* file;
  if (mode == IO_WRITE)
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

static void silo_close_file(void* context, void* file)
{
  // In this method, we do all the writing.
  silo_context* silo = (silo_context*)context;

  DBClose(file);
}

static void silo_write_data(void* context, void* file, io_buffered_data_t* data)
{
}

static void silo_read_data(void* context, void* file, io_buffered_data_t* data)
{
}

static void silo_dtor(void* context)
{
  free(context);
}

io_interface_t* silo_io_new()
{
  silo_context* context = malloc(sizeof(silo_context));
  io_vtable vtable = {.create_file = &silo_create_file,
                      .open_file = &silo_open_file, 
                      .close_file = &silo_close_file,
                      .read_data = &silo_read_data,
                      .write_data = &silo_write_data,
                      .dtor = &silo_dtor};
  return io_interface_new(context, "SILO", vtable);
}

#ifdef __cplusplus
}
#endif

