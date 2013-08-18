// Copyright 2012-2013 Jeffrey Johnson.
// 
// This file is part of Polymec, and is licensed under the Apache License, 
// Version 2.0 (the "License"); you may not use this file except in 
// compliance with the License. You may may find the text of the license in 
// the LICENSE file at the top-level source directory, or obtain a copy of 
// it at
// 
// http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef POLYMEC_SILO_IO_H
#define POLYMEC_SILO_IO_H

#include "core/io.h"

// Creates a Silo I/O interface designed for dumping simulation results 
// and reading them in to restart simulations.
io_interface_t* silo_io_new(MPI_Comm comm, int num_files, int mpi_tag);

// Creates a Silo I/O interface designed for dumping plot files.
io_interface_t* silo_plot_io_new(MPI_Comm comm, int num_files, int mpi_tag);

#endif

