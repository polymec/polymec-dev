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

void* silo_create_file(const char* filename,
                       const char* dirname,
                       void* userData)
{
  int driver = DB_HDF5;
  DBfile* file = DBCreate(filename, 0, DB_LOCAL, 0, driver);
  DBMkDir(file, dirname);
  DBSetDir(file, dirname);
  return (void*)file;
}

void* silo_open_file(const char* filename, 
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

void silo_close_file(void* context, void* file)
{
  // In this method, we do all the writing.
  silo_context* silo = (silo_context*)context;

  DBClose(file);
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
  io_vtable vtable = {.create_file = &silo_create_file,
                      .open_file = &silo_open_file, 
                      .close_file = &silo_close_file,
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

