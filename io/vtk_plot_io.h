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

#ifndef POLYMEC_VTK_PLOT_IO_H
#define POLYMEC_VTK_PLOT_IO_H

#include "core/io.h"

// Creates a VTK I/O interface designed for dumping (XML) plot files.
// This particular plotter dumps one file per process.
// If binary is set to true, the binary VTK XML format is used.
io_interface_t* vtk_plot_io_new(MPI_Comm comm, int mpi_tag, bool binary);

#endif

