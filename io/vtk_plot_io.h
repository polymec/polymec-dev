#ifndef POLYMEC_VTK_PLOT_IO_H
#define POLYMEC_VTK_PLOT_IO_H

#include "core/io.h"

#ifdef __cplusplus
extern "C" {
#endif

// Creates a VTK I/O interface designed for dumping (XML) plot files.
// This particular plotter dumps one file per process.
// If binary is set to true, the binary VTK XML format is used.
io_interface_t* vtk_plot_io_new(MPI_Comm comm, int mpi_tag, bool binary);

#ifdef __cplusplus
}
#endif

#endif

