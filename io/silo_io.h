#ifndef POLYMEC_SILO_IO_H
#define POLYMEC_SILO_IO_H

#include "core/io.h"

#ifdef __cplusplus
extern "C" {
#endif

// Creates a Silo I/O interface designed for dumping simulation results 
// and reading them in to restart simulations.
io_interface_t* silo_io_new(MPI_Comm comm, int num_files, int mpi_tag);

// Creates a Silo I/O interface designed for dumping plot files.
io_interface_t* silo_plot_io_new(MPI_Comm comm, int num_files, int mpi_tag);

#ifdef __cplusplus
}
#endif

#endif

