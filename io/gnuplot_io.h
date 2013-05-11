#ifndef POLYMEC_GNUPLOT_IO_H
#define POLYMEC_GNUPLOT_IO_H

#include "core/io.h"

// Creates a Gnuplot I/O interface designed for dumping text data that 
// is easily viewed with Gnuplot. This may only be used for serial runs.
io_interface_t* gnuplot_io_new();

#endif

