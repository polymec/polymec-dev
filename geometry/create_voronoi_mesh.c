// Copyright (c) 2012-2014, Jeffrey N. Johnson
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without 
// modification, are permitted provided that the following conditions are met:
// 
// 1. Redistributions of source code must retain the above copyright notice, this 
// list of conditions and the following disclaimer.
// 
// 2. Redistributions in binary form must reproduce the above copyright notice, 
// this list of conditions and the following disclaimer in the documentation 
// and/or other materials provided with the distribution.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE 
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

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
  bool made_dir = make_temp_dir(template, dir_path);
  if (!made_dir)
    polymec_error("create_voronoi_mesh: Could not open temporary directory.");

  // Make the input .node file containing the generators.
  char input_file[FILENAME_MAX];
  join_paths(dir_path, "voronoi.nodes", input_file);

  polymec_not_implemented("create_voronoi_mesh");
}

