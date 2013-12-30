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

#include "core/symmetry_1d.h"
#include "core/create_rectilinear_mesh.h"
#include "core/create_uniform_mesh.h"
#include "core/cubic_lattice.h"

// This helper function computes the logarithmic spacing for N points spanning
// the interval [x1, x2] with the given log factor, placing those points into 
// the array points.
static void compute_log_spacing(double x1, double x2, double log_factor, int N, double* points)
{
  // Figure out the innermost cell length.
  double sum = 0.0;
  for (int i = 0; i < N; ++i)
    sum += pow(log_factor, 1.0*i);
  double dx0 = (x2 - x1) / sum;

  // Now compute the rest of the spacings.
  double last_dx = dx0;
  for (int i = 0; i < N-1; ++i)
  {
    double dx = (i < N-1) ? pow(log_factor, 1.0*i) * dx0
      : last_dx;
    points[i] = (i == 0) ? dx0 : points[i-1] + dx;
    last_dx = dx;
  }
}

mesh_t* create_uniform_cartesian_1d_mesh(MPI_Comm comm, double x1, double x2, int N)
{
  ASSERT(x2 > x1);
  ASSERT(N > 0);

  // This is really just a nonuniform Cartesian mesh with, um, uniform spacing.
  double* xs = malloc(sizeof(double) * N);
  xs[0] = x1;
  double dx = (x2 - x1) / N;
  for (int i = 1; i < N; ++i)
    xs[i] = xs[i-1] + dx;
  mesh_t* mesh = create_nonuniform_cartesian_1d_mesh(comm, xs, N);
  free(xs);

  return mesh;
}

mesh_t* create_logarithmic_cartesian_1d_mesh(MPI_Comm comm, double x1, double x2, double log_factor, int N)
{
  // This is really just a nonuniform Cartesian mesh with logarithmic spacings.
  double* xs = malloc(sizeof(double) * N);
  compute_log_spacing(x1, x2, log_factor, N, xs);
  mesh_t* mesh = create_nonuniform_cartesian_1d_mesh(comm, xs, N);
  free(xs);

  return mesh;
}

mesh_t* create_nonuniform_cartesian_1d_mesh(MPI_Comm comm, double* xs, int N)
{
  // This is really just a 1D rectilinear mesh with the given spacings.
  double ys[2] = {-0.5, 0.5};
  double zs[2] = {-0.5, 0.5};
  mesh_t* mesh = create_rectilinear_mesh(comm, xs, N, ys, 2, zs, 2);

  // Tag the boundary faces.
  tag_rectilinear_mesh_faces(mesh, N, 1, 1, "x1", "x2", "y1", "y2", "z1", "z2");

  // Add some symmetry-related features.
  mesh_add_feature(mesh, SYMMETRIC);
  mesh_add_feature(mesh, ONE_DIMENSIONAL);
  mesh_add_feature(mesh, CARTESIAN);

  return mesh;
}

mesh_t* create_uniform_cylindrical_1d_mesh(MPI_Comm comm, double r1, double r2, int N)
{
  ASSERT(r1 > 0.0);
  ASSERT(r2 > r1);
  ASSERT(N > 0);

  // This is really just a nonuniform cylindrical mesh with, um, uniform spacing.
  double* rs = malloc(sizeof(double) * N);
  mesh_t* mesh = create_nonuniform_cylindrical_1d_mesh(comm, rs, N);
  free(rs);

  return mesh;
}

mesh_t* create_logarithmic_cylindrical_1d_mesh(MPI_Comm comm, double r1, double r2, double log_factor, int N)
{
  ASSERT(r1 > 0.0);
  ASSERT(r2 > r1);
  ASSERT(N > 0);

  // This is really just a nonuniform cylindrical mesh with logarithmic spacings.
  double* rs = malloc(sizeof(double) * N);
  compute_log_spacing(r1, r2, log_factor, N, rs);
  mesh_t* mesh = create_nonuniform_cylindrical_1d_mesh(comm, rs, N);
  free(rs);

  return mesh;
}

mesh_t* create_nonuniform_cylindrical_1d_mesh(MPI_Comm comm, double* rs, int N)
{
  ASSERT(N > 0);

#ifndef NDEBUG
  // Make sure the radial spacings make sense.
  ASSERT(rs[0] >= 0.0);
  for (int i = 1; i < N; ++i)
  {
    ASSERT(rs[i] > rs[i-1]);
  }
#endif

  // Now we make the mesh, starting from a uniform mesh and performing 
  // any needed surgery. This is simply the easiest way to get the topology 
  // correct.
  bbox_t bbox = {.x1 = rs[0], .x2 = rs[N-1], 
                 .y1 = -1.0, .y2 = 1.0, 
                 .z1 = -1.0, .z2 = 1.0}; // Wrong, but it doesn't matter.
  mesh_t* mesh = create_uniform_mesh(comm, N, 1, 1, &bbox);

  // Retrieve the (topological) cubic lattice for the mesh.
  cubic_lattice_t* lattice = mesh_property(mesh, "lattice");
  ASSERT(lattice != NULL);

  if (rs[0] == 0.0)
  {
    // If the inner radius (rs[0]) is zero, the first cell will be a wedge 
    // whose edge lies at the origin. We have to collapse the leftmost face to 
    // an edge in this case.

    // Unhook the first face from the first cell. I happen to know that the 
    // -x face in the first cell is the first face in the mesh, so we can 
    // accomplish this unhooking by simply starting the list of faces for 
    // the first cell at 1 instead of 0.
    mesh->cell_face_offsets[0] = 1;

    // Now we effectively merge the vertical edges of that first face into 
    // one by changing the index of the edge for the -y face. Since edges 
    // are constructed in a "black box" fashion from pairs of nodes, we have 
    // to be more careful about this.
    // FIXME
  }
  else
  {
    // Otherwise, the mesh will topologically be rectilinear and we only need 
    // to set the geometry.
    // FIXME
  }

  // Compute the geometry of the mesh.
  mesh_compute_geometry(mesh);

  // Tag the boundary faces.
  tag_rectilinear_mesh_faces(mesh, N, 1, 1, "r1", "r2", "theta1", "theta2", "z1", "z2");

  // Add some symmetry-related features.
  mesh_add_feature(mesh, SYMMETRIC);
  mesh_add_feature(mesh, ONE_DIMENSIONAL);
  mesh_add_feature(mesh, CYLINDRICAL);

  return mesh;
}

mesh_t* create_uniform_spherical_1d_mesh(MPI_Comm comm, double r1, double r2, int N)
{
  ASSERT(r1 > 0.0);
  ASSERT(r2 > r1);
  ASSERT(N > 0);

  // This is really just a nonuniform spherical mesh with, um, uniform spacing.
  double* rs = malloc(sizeof(double) * N);
  mesh_t* mesh = create_nonuniform_spherical_1d_mesh(comm, rs, N);
  free(rs);

  return mesh;
}

mesh_t* create_logarithmic_spherical_1d_mesh(MPI_Comm comm, double r1, double r2, double log_factor, int N)
{
  ASSERT(r1 > 0.0);
  ASSERT(r2 > r1);
  ASSERT(N > 0);

  // This is really just a nonuniform spherical mesh with logarithmic spacings.
  double* rs = malloc(sizeof(double) * N);
  compute_log_spacing(r1, r2, log_factor, N, rs);
  mesh_t* mesh = create_nonuniform_spherical_1d_mesh(comm, rs, N);
  free(rs);

  return mesh;
}

mesh_t* create_nonuniform_spherical_1d_mesh(MPI_Comm comm, double* rs, int N)
{
  mesh_t* mesh = NULL;
  // FIXME

  // Add some symmetry-related features.
  mesh_add_feature(mesh, SYMMETRIC);
  mesh_add_feature(mesh, ONE_DIMENSIONAL);
  mesh_add_feature(mesh, SPHERICAL);

  return mesh;
}

