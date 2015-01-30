// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/partition_mesh.h"
#include "core/silo_file.h"
#include "geometry/cubic_lattice.h"
#include "geometry/create_uniform_mesh.h"
#include "geometry/create_cubed_cylinder_mesh.h"
#include "geometry/welded_block_mesh.h"

static void create_radial_blocks(int nx, int nz, 
                                 real_t R, real_t L,
                                 real_t l, real_t k,
                                 mesh_t** blocks)
{
  ASSERT(k == 0.0); // For now.
  // Generic bounding box -- doesn't really matter.
  bbox_t bbox = {.x1 = -0.5, .x2 = 0.5,
                 .y1 = -0.5, .y2 = 0.5,
                 .z1 = -0.5, .z2 = 0.5};

  // Indexing mechanism.
  cubic_lattice_t* lattice = cubic_lattice_new(nx, nx, nz);

  // Constant spacings.
  real_t dx = l / nx, dy = l / nx;
  real_t dz = L / nz;
  real_t dtheta = 0.5*M_PI / nx;

  // -x block
  blocks[0] = create_uniform_mesh(MPI_COMM_SELF, nx, nx, nz, &bbox);
  tag_rectilinear_mesh_faces(blocks[0], 
                             "west_outer", "west_seam",
                             "southwest_seam", "northwest_seam",
                             "west_bottom", "west_top");
  for (int k = 0; k <= nz; ++k)
  {
    real_t zk = -0.5*L + k*dz;
    for (int j = 0; j <= nx; ++j)
    {
      // Compute the radial spacing for this j index.
      real_t theta = 1.25*M_PI - j*dtheta;
      real_t cos_theta = cos(theta), sin_theta = sin(theta);

      // Find cx, the point of nearest approach on the center block surface.
      point_t xc = {.x = -0.5*l, .y = -0.5*l + j*dy, .z = 0.0}; 

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
        blocks[0]->nodes[n] = xn;
      }
    }
  }
  mesh_compute_geometry(blocks[0]);

  // -y block
  blocks[1] = create_uniform_mesh(MPI_COMM_SELF, nx, nx, nz, &bbox);
  tag_rectilinear_mesh_faces(blocks[1], 
                             "southwest_seam", "southeast_seam",
                             "south_outer", "south_seam",
                             "south_bottom", "south_top");
  for (int k = 0; k <= nz; ++k)
  {
    real_t zk = -0.5*L + k*dz;
    for (int i = 0; i <= nx; ++i)
    {
      // Compute the radial spacing for this j index.
      real_t theta = 1.25*M_PI + i*dtheta;
      real_t cos_theta = cos(theta), sin_theta = sin(theta);

      // Find cx, the point of nearest approach on the center block surface.
      point_t xc = {.x = -0.5*l + i*dx, .y = -0.5*l, .z = 0.0}; 

      // Find xR, the point on the outside of the cylinder for this j.
      point_t xR = {.x = R*cos_theta, .y = R*sin_theta, .z = 0.0};

      // Now find dR, the increment of the vector that connects xc to xR.
      vector_t dR;
      point_displacement(&xc, &xR, &dR);
      dR.x /= nx;
      dR.y /= nx;

      // Compute the node positions, proceeding from lower left to upper right.
      for (int j = 0; j <= nx; ++j)
      {
        int n = (int)cubic_lattice_node(lattice, i, j, k);
        point_t xn = {.x = xR.x - j*dR.x, .y = xR.y - j*dR.y, .z = zk};
        blocks[1]->nodes[n] = xn;
      }
    }
  }
  mesh_compute_geometry(blocks[1]);

  // +x block
  blocks[2] = create_uniform_mesh(MPI_COMM_SELF, nx, nx, nz, &bbox);
  tag_rectilinear_mesh_faces(blocks[2], 
                             "east_seam", "east_outer",
                             "southeast_seam", "northeast_seam",
                             "east_bottom", "east_top");
  for (int k = 0; k <= nz; ++k)
  {
    real_t zk = -0.5*L + k*dz;
    for (int j = 0; j <= nx; ++j)
    {
      // Compute the radial spacing for this j index.
      real_t theta = 1.75*M_PI + j*dtheta;
      real_t cos_theta = cos(theta), sin_theta = sin(theta);

      // Find cx, the point of nearest approach on the center block surface.
      point_t xc = {.x = 0.5*l, .y = -0.5*l + j*dy, .z = 0.0}; 

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
        point_t xn = {.x = xc.x + i*dR.x, .y = xc.y + i*dR.y, .z = zk};
        blocks[2]->nodes[n] = xn;
      }
    }
  }
  mesh_compute_geometry(blocks[2]);

  // +y block
  blocks[3] = create_uniform_mesh(MPI_COMM_SELF, nx, nx, nz, &bbox);
  tag_rectilinear_mesh_faces(blocks[3], 
                             "northwest_seam", "northeast_seam",
                             "north_outer", "north_seam",
                             "north_bottom", "north_top");
  for (int k = 0; k <= nz; ++k)
  {
    real_t zk = -0.5*L + k*dz;
    for (int i = 0; i <= nx; ++i)
    {
      // Compute the radial spacing for this j index.
      real_t theta = 0.75*M_PI - i*dtheta;
      real_t cos_theta = cos(theta), sin_theta = sin(theta);

      // Find cx, the point of nearest approach on the center block surface.
      point_t xc = {.x = -0.5*l + i*dx, .y = 0.5*l, .z = 0.0}; 

      // Find xR, the point on the outside of the cylinder for this j.
      point_t xR = {.x = R*cos_theta, .y = R*sin_theta, .z = 0.0};

      // Now find dR, the increment of the vector that connects xc to xR.
      vector_t dR;
      point_displacement(&xc, &xR, &dR);
      dR.x /= nx;
      dR.y /= nx;

      // Compute the node positions, proceeding from lower left to upper right.
      for (int j = 0; j <= nx; ++j)
      {
        int n = (int)cubic_lattice_node(lattice, i, j, k);
        point_t xn = {.x = xc.x + j*dR.x, .y = xc.y + j*dR.y, .z = zk};
        blocks[3]->nodes[n] = xn;
      }
    }
  }
  mesh_compute_geometry(blocks[3]);
}

mesh_t* create_cubed_cylinder_mesh(MPI_Comm comm, 
                                   int nx, int nz,
                                   real_t R, real_t L,
                                   real_t l, real_t k)
{
  ASSERT(nx > 0);
  ASSERT(nz > 0);
  ASSERT(R > 0.0);
  ASSERT(L > 0.0);
  ASSERT(l > 0.0);
  ASSERT(l < R);
  ASSERT(k >= 0.0);

  // Construct the center block.
  bbox_t bbox = {.x1 = -0.5*l, .x2 = 0.5*l,
                 .y1 = -0.5*l, .y2 = 0.5*l,
                 .z1 = -0.5*L, .z2 = 0.5*L};
  mesh_t* center_block = create_uniform_mesh(MPI_COMM_SELF, nx, nx, nz, &bbox);
  tag_rectilinear_mesh_faces(center_block, 
                             "west_seam", "east_seam",
                             "south_seam", "north_seam",
                             "center_bottom", "center_top");

  // If needed, deform its outer cells.
  if (k != 0.0)
  {
    // FIXME!
  }

  // Construct the radial blocks.
  mesh_t* radial_blocks[4];
  create_radial_blocks(nx, nz, R, L, l, k, radial_blocks);

  // Weld'em blocks.
  mesh_t* blocks[5] = {center_block, radial_blocks[0], radial_blocks[1], 
                       radial_blocks[2], radial_blocks[3]};
  mesh_t* mesh = welded_block_mesh(blocks, 5, 1e-10);

  // Clean up.
  for (int i = 0; i < 4; ++i)
    mesh_free(radial_blocks[i]);

  // Now partition the thing if we've been asked to.
  if (comm != MPI_COMM_SELF)
    partition_mesh(&mesh, comm, NULL, 0.0);

  return mesh;
}

mesh_t* create_cubed_cylindrical_shell_mesh(MPI_Comm comm, 
                                            int nx, int nz,
                                            real_t r, real_t R, real_t L)
{
  ASSERT(nx > 0);
  ASSERT(nz > 0);
  ASSERT(r > 1e-8);
  ASSERT(R > 1e-8);
  ASSERT(r < R);
  ASSERT(L > 0.0);

  // Construct the curvature and length of the inner interface.
  real_t k = 1.0/r;
  real_t l = 1.0/(k * sqrt(2.0));

  // Construct the radial blocks.
  mesh_t* radial_blocks[4];
  create_radial_blocks(nx, nz, R, L, l, k, radial_blocks);

  // Weld'em blocks.
  mesh_t* mesh = welded_block_mesh(radial_blocks, 4, 1e-10);

  // Clean up.
  for (int i = 0; i < 4; ++i)
    mesh_free(radial_blocks[i]);

  // Now partition the thing if we've been asked to.
  if (comm != MPI_COMM_SELF)
    partition_mesh(&mesh, comm, NULL, 0.0);

  return mesh;
}

