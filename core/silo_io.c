#include <stdlib.h>
#include "core/silo_io.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct
{
} silo_context;

static int silo_open(void* context, const char* prefix, const char* directory, io_mode_t mode, MPI_Comm comm, int num_files, int mpi_tag)
{
  return 0;
}

static int silo_close(void* context)
{
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

