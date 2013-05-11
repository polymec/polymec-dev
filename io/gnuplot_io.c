#include <stdarg.h>
#include <string.h>
#include <math.h>
#include "io/gnuplot_io.h"
#include "core/point.h"
#include "core/slist.h"

#ifdef USE_MPI
#include <mpi.h>
#endif

static void* gnuplot_create_file(void* context, 
                                 const char* filename,
                                 const char* dirname)
{
  // Open a file for writing.
  FILE* fd = fopen(filename, "w");
  return (void*)fd;
}

static void* gnuplot_open_file(void* context, 
                               const char* filename, 
                               const char* dirname,
                               io_mode_t mode)
{
  FILE* fd = fopen(filename, "w");
  return (void*)fd;
}

static void gnuplot_close_file(void* context, void* file)
{
  // Close the file, writing any buffers.
  fclose((FILE*)file);
}

static int gnuplot_get_num_datasets(void* context, void* file, int* num_datasets)
{
  return 1;
}

// ASCI version
static void gnuplot_plot_write_datasets(void* context, void* f, io_dataset_t** datasets, int num_datasets, int rank_in_group, int procs_per_file)
{
  ASSERT(procs_per_file == 1);
  FILE* fd = (FILE*)f;

  // Only cell-centered fields are supported.
  io_dataset_t* dataset = datasets[0];
  int pos = 0, num_comps;
  char *field_name;
  double *field;
  mesh_centering_t centering;
  int num_fields = io_dataset_num_fields(dataset), i = 0;
  double* fields[num_fields];
  int field_comps[num_fields];
  while (io_dataset_next_field(dataset, &pos, &field_name, &field, &num_comps, &centering))
  {
    if (centering != MESH_CELL)
      polymec_error("Non-cell-centered field found: %s\nThe Gnuplot I/O interface only supports cell-centered fields.", field_name);
    fields[i] = field;
    field_comps[i] = num_comps;
    ++i;
  }

  // Write a header.
  fprintf(fd, "# cell x y z volume ");
  pos = 0;
  while (io_dataset_next_field(dataset, &pos, &field_name, &field, &num_comps, &centering))
  {
    if (num_comps == 1) // Scalar fields only for now!
      fprintf(fd, "%s ", field_name);
  }
  fprintf(fd, "\n");

  // Now write cell values.
  mesh_t* mesh = io_dataset_get_mesh(dataset);
  for (int c = 0; c < mesh->num_cells; ++c)
  {
    cell_t* cell = &mesh->cells[c];
    fprintf(fd, "%d %g %g %g %g ", c, cell->center.x, cell->center.y, cell->center.z, cell->volume);

    for (int i = 0; i < num_fields; ++i)
    {
      if (field_comps[i] == 1)
        fprintf(fd, "%g ", fields[i][c]);
    }
    fprintf(fd, "\n");
  }
}

io_interface_t* gnuplot_io_new()
{
  int nproc;
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  if (nproc > 1)
    polymec_error("The Gnuplot I/O interface is only intended for serial runs.");

  io_vtable vtable = {.create_file = gnuplot_create_file,
                      .open_file = gnuplot_open_file, 
                      .close_file = gnuplot_close_file,
                      .get_num_datasets = gnuplot_get_num_datasets,
                      .write_datasets = gnuplot_plot_write_datasets};
  return io_interface_new_serial(NULL, "Gnuplot", "gnuplot", vtable);
}

