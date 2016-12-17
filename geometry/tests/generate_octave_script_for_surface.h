// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/point.h"

#ifndef POLYMEC_GENERATE_OCTAVE_SCRIPT_FOR_SURFACE_H
#define POLYMEC_GENERATE_OCTAVE_SCRIPT_FOR_SURFACE_H

// This function generates an Octave visualization script for the 
// given surface, sampling its values according to the number of samples 
// (the same in x, y, and z) within the given bounding box.
static void generate_octave_script_for_surface(sd_func_t* surface, 
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
        F = sd_func_value(surface, &x);
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

