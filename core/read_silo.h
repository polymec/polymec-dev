// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_READ_SILO_H
#define POLYMEC_READ_SILO_H

#include "core/polymec.h"
#include "core/mesh.h"
#include "core/unordered_map.h"

// Reads a silo mesh file, placing the mesh, cell-centered fields, and the 
// simulation time into their respective variables, and also reading any 
// stored mesh tags. Mesh is allocated by this function, whereas fields must 
// be a pre-allocated unordered map.
void read_silo_mesh(MPI_Comm comm,
                    const char* file_prefix,
                    const char* directory,
                    int cycle,
                    int num_files,
                    int mpi_tag,
                    mesh_t** mesh,
                    string_ptr_unordered_map_t* fields,
                    real_t* time);

#endif
