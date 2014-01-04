// Copyright (c) 2012-2013, Jeffrey N. Johnson
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

#include "core/point.h"

#ifndef POLYMEC_GENERATE_OCTAVE_SCRIPT_FOR_SURFACE_H
#define POLYMEC_GENERATE_OCTAVE_SCRIPT_FOR_SURFACE_H

// This function generates an Octave visualization script for the 
// given surface, sampling its values according to the number of samples 
// (the same in x, y, and z) within the given bounding box.
void generate_octave_script_for_surface(sp_func_t* surface, 
                                        int num_samples, 
                                        bbox_t* bounding_box,
                                        const char* script_name)
{
  // If you're a process whose rank != 0, bug off.
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank != 0) return;

  // Create a text file containing an Octave script that can be run to 
  // visualize this plot.
  FILE* fd = fopen(script_name, "w");
  fprintf(fd, "%% %s - A script for visualizing a sphere.\n", script_name);
  fprintf(fd, "%% Run with octave --persist %s\n\n", script_name);

  // Write out X, Y, and Z arrays.
  real_t hx = (bounding_box->x2 - bounding_box->x1) / num_samples;
  real_t hy = (bounding_box->y2 - bounding_box->y1) / num_samples;
  real_t hz = (bounding_box->z2 - bounding_box->z1) / num_samples;
  fprintf(fd, "X = [");
  for (int i = 0; i < num_samples; ++i)
    fprintf(fd, "%g ", bounding_box->x1 + (i+0.5)*hx);
  fprintf(fd, "];\n");

  fprintf(fd, "Y = [");
  for (int i = 0; i < num_samples; ++i)
    fprintf(fd, "%g ", bounding_box->y1 + (i+0.5)*hy);
  fprintf(fd, "];\n");

  fprintf(fd, "Z = [");
  for (int i = 0; i < num_samples; ++i)
    fprintf(fd, "%g ", bounding_box->z1 + (i+0.5)*hz);
  fprintf(fd, "];\n");

  // Write out the F array.
  point_t x;
  fprintf(fd, "F = zeros(%d, %d, %d);\n", num_samples, num_samples, num_samples);
  for (int k = 0; k < num_samples; ++k)
  {
    x.z = -1.0 + (k+0.5)*hz;
    for (int j = 0; j < num_samples; ++j)
    {
      x.y = -1.0 + (j+0.5)*hy;
      for (int i = 0; i < num_samples; ++i)
      {
        x.x = -1.0 + (i+0.5)*hx;
        real_t F;
        sp_func_eval(surface, &x, &F);
        fprintf(fd, "F(%d, %d, %d) = %g;\n", i+1, j+1, k+1, F);
      }
    }
  }

  // Dump out the rest of the script for Octave.
  static const char* instructions = 
    "[XX, YY, ZZ] = meshgrid(X, Y, Z);\n"
    "isosurface(XX, YY, ZZ, F, 0.0); \n"
//    "set(gca, 'PlotBoxAspectRatioMode', 'manual', 'PlotBoxAspectRatio', [1 1 1]);\n"
    "xlabel('x');\n"
    "ylabel('y');\n"
    "zlabel('z');\n";
  fprintf(fd, "%s", instructions);
  fclose(fd);
}

#endif

