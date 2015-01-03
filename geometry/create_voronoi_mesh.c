// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "geometry/create_voronoi_mesh.h"

mesh_t* create_voronoi_mesh(MPI_Comm comm, point_t* generators, 
                            int num_generators, bbox_t* bounding_box)
{
  // Check to see whether tetgen is installed.
  int status = system("tetgen");
  if (status != 0)
    polymec_error("create_voronoi_mesh: tetgen must be installed in your PATH.");

  // Generate a temporary directory for working with tetgen.
  char template[FILENAME_MAX], dir_path[FILENAME_MAX];
  sprintf(template, "voronoi-XXXXXX");
  bool made_dir = make_temp_directory(template, dir_path);
  if (!made_dir)
    polymec_error("create_voronoi_mesh: Could not open temporary directory.");

  // Make the input .node file containing the generators.
  char input_file[FILENAME_MAX];
  join_paths(dir_path, "voronoi.nodes", input_file);

  polymec_not_implemented("create_voronoi_mesh");
}

