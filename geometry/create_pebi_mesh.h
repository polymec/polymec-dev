// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_PEBI_H
#define POLYMEC_PEBI_H

#include "core/mesh.h"

// Perpendicular-bisector mesh features.
const char* PEBI;       // Is a PEBI mesh (has no edges or nodes)

// Creates a PEBI mesh given a set of cell centers and adjacency (face) 
// information. The faces array contains 2*num_faces entries, where 
// faces[2*i] and faces[2*i+1] contain the indices of the cells connected by 
// face i. If a face is only connected to one cell, faces[2*i+1] == -1.
// Meanwhile, face_areas contains num_faces entries, with face_areas[i] holding
// the area of face i. No edge or node information is stored.
// NOTE: If face_centers is NULL, each face center is assumed to lie at the 
// midpoint between its two cells. 
mesh_t* create_pebi_mesh(MPI_Comm comm, 
                         point_t* cell_centers, real_t* cell_volumes, int num_cells,
                         int* faces, real_t* face_areas, point_t* face_centers,
                         int num_faces);
 
// Creates a PEBI mesh from the given unstructured mesh. Features are 
// copied--tags are not.
mesh_t* create_pebi_mesh_from_unstructured_mesh(mesh_t* mesh);

#endif

