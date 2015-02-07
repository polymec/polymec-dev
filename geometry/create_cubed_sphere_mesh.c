// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/partition_mesh.h"
#include "geometry/cubic_lattice.h"
#include "geometry/create_uniform_mesh.h"
#include "geometry/create_cubed_sphere_mesh.h"
#include "geometry/create_welded_block_mesh.h"

// This transformation maps a point x within the cube 
// [-L/2, -L/2, -L/2] x [L/2, L/2, L/2] to its corresponding point y in the 
// sphere of radius L/2.
static void map_to_sphere(real_t L, point_t* x, point_t* y)
{
  // First, map x to a point (X, Y, Z) on the unit square.
  real_t X = 2.0*x->x/L, Y = 2.0*x->y/L, Z = 2.0*x->z/L;

  // Now map (X, Y, Z) to (X1, Y1, Z1) in the unit circle.
  // See http://mathproofs.blogspot.com/2005/07/mapping-cube-to-sphere.html.
  real_t one_third = 1.0/3.0;
  real_t X1 = X * sqrt(1.0 - 0.5*Y*Y - 0.5*Z*Z - one_third*Y*Y*Z*Z);
  real_t Y1 = Y * sqrt(1.0 - 0.5*Z*Z - 0.5*X*X - one_third*Z*Z*X*X);
  real_t Z1 = Z * sqrt(1.0 - 0.5*X*X - 0.5*Y*Y - one_third*X*X*Y*Y);

  // Finally, map to the circle with radius L/2.
  real_t theta = acos(Z1);
  real_t phi = atan2(Y1, X1);

  y->x = 0.5*L * cos(phi) * sin(theta);
  y->y = 0.5*L * sin(phi) * sin(theta);
  y->z = 0.5*L * cos(theta);
}

mesh_t* create_cubed_sphere_mesh(MPI_Comm comm,
                                 int n, real_t R, real_t l, 
                                 bool curved_center_block,
                                 const char* R_tag)
{
  return NULL; // FIXME
}

mesh_t* create_cubed_spherical_shell_mesh(MPI_Comm comm,
                                          int ns, int nr,
                                          real_t r, real_t R,
                                          const char* r_tag,
                                          const char* R_tag)
{
  return NULL; // FIXME
}

mesh_t* create_cubed_sphere_panel(MPI_Comm comm,
                                  int ns, int nr,
                                  real_t R, real_t l,
                                  bool curved_bottom,
                                  cubed_sphere_panel_t which_panel)
{
  ASSERT(ns > 0);
  ASSERT(nr > 0);
  ASSERT(R >= 0.0);
  ASSERT(l >= 0.0);
  ASSERT(l < R);

  // Generic bounding box -- doesn't really matter.
  bbox_t bbox = {.x1 = -0.5, .x2 = 0.5,
                 .y1 = -0.5, .y2 = 0.5,
                 .z1 = -0.5, .z2 = 0.5};

  // Indexing mechanism.
  cubic_lattice_t* lattice = cubic_lattice_new(ns, ns, nr);

#if 0
  // Constant spacings.
  real_t dx = l / nx, dy = l / nx, dz = ;
  real_t dz = L / nz;
  real_t dtheta = 0.5*M_PI / nx;
#endif

  mesh_t* panel = create_uniform_mesh(MPI_COMM_SELF, ns, ns, nr, &bbox);

  // Determine a prefix for the panel's tags.
  char tag_prefix[1024];
  if (which_panel == X1_PANEL)
    strcpy(tag_prefix, "x1_panel_");
  else if (which_panel == X2_PANEL)
    strcpy(tag_prefix, "x2_panel_");
  else if (which_panel == Y1_PANEL)
    strcpy(tag_prefix, "y1_panel_");
  else if (which_panel == Y2_PANEL)
    strcpy(tag_prefix, "y2_panel_");
  else if (which_panel == Z1_PANEL)
    strcpy(tag_prefix, "z1_panel_");
  else 
  {
    ASSERT(which_panel == Z2_PANEL);
    strcpy(tag_prefix, "z2_panel_");
  }

  // Tag the panel's faces accordingly.
  char tags[6][1024];
  for (int i = 0; i < 6; ++i)
    sprintf(tags[i], "%s_%d", tag_prefix, i);
  tag_rectilinear_mesh_faces(panel, tags[0], tags[1], tags[2], tags[3],
                             tags[4], tags[5]);
#if 0
  for (int k = 0; k <= nz; ++k)
  {
    real_t zk = -0.5*L + k*dz;
    for (int j = 0; j <= nx; ++j)
    {
      // Compute the radial spacing for this j index.
      real_t theta = 1.25*M_PI - j*dtheta;
      real_t cos_theta = cos(theta), sin_theta = sin(theta);

      // Find xc, the point of nearest approach on the center block surface.
      point_t xc = {.x = -0.5*l, .y = -0.5*l + j*dy, .z = 0.0}; 
      if (curved_center_block)
      {
        point_t x = xc;
        map_to_circle(l, &x, &xc);
      }

      // Find xR, the point on the outside of the cylinder for this j.
      point_t xR = {.x = R*cos_theta, .y = R*sin_theta, .z = 0.0};

      // Now find dR, the increment of the vector that connects xc to xR.
      vector_t dR;
      point_displacement(&xc, &xR, &dR);
      dR.x /= nx;
      dR.y /= nx;

      // Compute the node positions, proceeding from lower left to upper right.
      for (int i = 0; i <= nx; ++i)
      {
        int n = (int)cubic_lattice_node(lattice, i, j, k);
        point_t xn = {.x = xR.x - i*dR.x, .y = xR.y - i*dR.y, .z = zk};
        panel->nodes[n] = xn;
      }
    }
  }
#endif
  mesh_compute_geometry(panel);

  return panel;
}

