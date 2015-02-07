// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/partition_mesh.h"
#include "geometry/cubic_lattice.h"
#include "geometry/create_uniform_mesh.h"
#include "geometry/create_cubed_cylinder_mesh.h"
#include "geometry/create_welded_block_mesh.h"

// This transformation maps a point x within the square 
// [-L/2, -L/2] x [L/2, L/2] to its corresponding point y in the circle of 
// radius L/2. The z coordinate of y is set to that of x.
static void map_to_circle(real_t L, point_t* x, point_t* y)
{
  // First, map x to a point (X, Y) on the unit square.
  real_t X = 2.0*x->x/L, Y = 2.0*x->y/L;

  // Now map (X, Y) to (X1, Y1) in the unit circle.
  // See http://mathproofs.blogspot.com/2005/07/mapping-square-to-circle.html.
  real_t X1 = X * sqrt(1.0 - 0.5*Y*Y);
  real_t Y1 = Y * sqrt(1.0 - 0.5*X*X);

  // Finally, map to the circle with radius L/2.
  real_t theta = atan2(Y1, X1);
  real_t r = sqrt(X1*X1 + Y1*Y1);

  y->x = 0.5*L * r * cos(theta);
  y->y = 0.5*L * r * sin(theta);
  y->z = x->z;
}

static void create_radial_blocks(int nx, int nz, 
                                 real_t R, real_t L, real_t l, 
                                 bool curved_center_block,
                                 mesh_t** blocks)
{
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

      // Find xc, the point of nearest approach on the center block surface.
      point_t xc = {.x = -0.5*l + i*dx, .y = -0.5*l, .z = 0.0}; 
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

      // Find xc, the point of nearest approach on the center block surface.
      point_t xc = {.x = 0.5*l, .y = -0.5*l + j*dy, .z = 0.0}; 
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

      // Find xc, the point of nearest approach on the center block surface.
      point_t xc = {.x = -0.5*l + i*dx, .y = 0.5*l, .z = 0.0}; 
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
                                   real_t R, real_t L, real_t l, 
                                   bool curved_center_block,
                                   const char* R_tag,
                                   const char* bottom_tag,
                                   const char* top_tag)
{
  ASSERT(nx > 0);
  ASSERT(nz > 0);
  ASSERT(R > 0.0);
  ASSERT(L > 0.0);
  ASSERT(l > 0.0);
  ASSERT(l < R);

  // Construct the center block.
  bbox_t bbox = {.x1 = -0.5*l, .x2 = 0.5*l,
                 .y1 = -0.5*l, .y2 = 0.5*l,
                 .z1 = -0.5*L, .z2 = 0.5*L};
  mesh_t* center_block = create_uniform_mesh(MPI_COMM_SELF, nx, nx, nz, &bbox);
  tag_rectilinear_mesh_faces(center_block, 
                             "west_seam", "east_seam",
                             "south_seam", "north_seam",
                             "center_bottom", "center_top");

  // If needed, deform the center block.
  if (curved_center_block)
  {
    cubic_lattice_t* lattice = cubic_lattice_new(nx, nx, nz);
    for (int i = 0; i <= nx; ++i)
    {
      for (int j = 0; j <= nx; ++j)
      {

        for (int k = 0; k <= nz; ++k)
        {
          int n = (int)cubic_lattice_node(lattice, i, j, k);
          point_t x = center_block->nodes[n];
          map_to_circle(l, &x, &(center_block->nodes[n]));
        }
      }
    }
    mesh_compute_geometry(center_block);
  }

  // Construct the radial blocks.
  mesh_t* radial_blocks[4];
  create_radial_blocks(nx, nz, R, L, l, curved_center_block, radial_blocks);

  // Weld'em blocks.
  mesh_t* blocks[5] = {center_block, radial_blocks[0], radial_blocks[1], 
                       radial_blocks[2], radial_blocks[3]};
  mesh_t* mesh = create_welded_block_mesh(blocks, 5, 1e-10);

  // Replace the tags.
  {
    int N1, N2, N3, N4;
    int* outer1 = mesh_tag(radial_blocks[0]->face_tags, "west_outer", &N1);
    int* outer2 = mesh_tag(radial_blocks[1]->face_tags, "south_outer", &N2);
    int* outer3 = mesh_tag(radial_blocks[2]->face_tags, "east_outer", &N3);
    int* outer4 = mesh_tag(radial_blocks[3]->face_tags, "north_outer", &N4);
    int* rtag = mesh_create_tag(mesh->face_tags, R_tag, N1+N2+N3+N4);
    memcpy(&rtag[0], outer1, sizeof(int) * N1);
    memcpy(&rtag[N1], outer2, sizeof(int) * N2);
    memcpy(&rtag[N2], outer3, sizeof(int) * N3);
    memcpy(&rtag[N3], outer4, sizeof(int) * N4);
  }

  {
    int N1, N2, N3, N4, N5;
    int* bot1 = mesh_tag(radial_blocks[0]->face_tags, "west_bottom", &N1);
    int* bot2 = mesh_tag(radial_blocks[1]->face_tags, "south_bottom", &N2);
    int* bot3 = mesh_tag(radial_blocks[2]->face_tags, "east_bottom", &N3);
    int* bot4 = mesh_tag(radial_blocks[3]->face_tags, "north_bottom", &N4);
    int* bot5 = mesh_tag(center_block->face_tags, "center_bottom", &N5);
    int* btag = mesh_create_tag(mesh->face_tags, bottom_tag, N1+N2+N3+N4+N5);
    memcpy(&btag[0], bot1, sizeof(int) * N1);
    memcpy(&btag[N1], bot2, sizeof(int) * N2);
    memcpy(&btag[N2], bot3, sizeof(int) * N3);
    memcpy(&btag[N3], bot4, sizeof(int) * N4);
    memcpy(&btag[N4], bot5, sizeof(int) * N5);
  }

  {
    int N1, N2, N3, N4, N5;
    int* top1 = mesh_tag(radial_blocks[0]->face_tags, "west_top", &N1);
    int* top2 = mesh_tag(radial_blocks[1]->face_tags, "south_top", &N2);
    int* top3 = mesh_tag(radial_blocks[2]->face_tags, "east_top", &N3);
    int* top4 = mesh_tag(radial_blocks[3]->face_tags, "north_top", &N4);
    int* top5 = mesh_tag(center_block->face_tags, "center_top", &N5);
    int* ttag = mesh_create_tag(mesh->face_tags, top_tag, N1+N2+N3+N4+N5);
    memcpy(&ttag[0], top1, sizeof(int) * N1);
    memcpy(&ttag[N1], top2, sizeof(int) * N2);
    memcpy(&ttag[N2], top3, sizeof(int) * N3);
    memcpy(&ttag[N3], top4, sizeof(int) * N4);
    memcpy(&ttag[N4], top5, sizeof(int) * N5);
  }

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
                                            real_t r, real_t R, real_t L,
                                            const char* r_tag,
                                            const char* R_tag,
                                            const char* bottom_tag,
                                            const char* top_tag)
{
  ASSERT(nx > 0);
  ASSERT(nz > 0);
  ASSERT(r > 1e-8);
  ASSERT(R > 1e-8);
  ASSERT(r < R);
  ASSERT(L > 0.0);

  // Construct the curvature and length of the inner interface.
  real_t l = r;

  // Construct the radial blocks.
  mesh_t* radial_blocks[4];
  create_radial_blocks(nx, nz, R, L, l, true, radial_blocks);

  // Weld'em blocks.
  mesh_t* mesh = create_welded_block_mesh(radial_blocks, 4, 1e-10);

  // Replace the tags.
  {
    int N1, N2, N3, N4;
    int* outer1 = mesh_tag(radial_blocks[0]->face_tags, "west_outer", &N1);
    int* outer2 = mesh_tag(radial_blocks[1]->face_tags, "south_outer", &N2);
    int* outer3 = mesh_tag(radial_blocks[2]->face_tags, "east_outer", &N3);
    int* outer4 = mesh_tag(radial_blocks[3]->face_tags, "north_outer", &N4);
    int* rtag = mesh_create_tag(mesh->face_tags, R_tag, N1+N2+N3+N4);
    memcpy(&rtag[0], outer1, sizeof(int) * N1);
    memcpy(&rtag[N1], outer2, sizeof(int) * N2);
    memcpy(&rtag[N2], outer3, sizeof(int) * N3);
    memcpy(&rtag[N3], outer4, sizeof(int) * N4);
  }

  {
    int N1, N2, N3, N4;
    int* bot1 = mesh_tag(radial_blocks[0]->face_tags, "west_bottom", &N1);
    int* bot2 = mesh_tag(radial_blocks[1]->face_tags, "south_bottom", &N2);
    int* bot3 = mesh_tag(radial_blocks[2]->face_tags, "east_bottom", &N3);
    int* bot4 = mesh_tag(radial_blocks[3]->face_tags, "north_bottom", &N4);
    int* btag = mesh_create_tag(mesh->face_tags, bottom_tag, N1+N2+N3+N4);
    memcpy(&btag[0], bot1, sizeof(int) * N1);
    memcpy(&btag[N1], bot2, sizeof(int) * N2);
    memcpy(&btag[N2], bot3, sizeof(int) * N3);
    memcpy(&btag[N3], bot4, sizeof(int) * N4);
  }

  {
    int N1, N2, N3, N4;
    int* top1 = mesh_tag(radial_blocks[0]->face_tags, "west_top", &N1);
    int* top2 = mesh_tag(radial_blocks[1]->face_tags, "south_top", &N2);
    int* top3 = mesh_tag(radial_blocks[2]->face_tags, "east_top", &N3);
    int* top4 = mesh_tag(radial_blocks[3]->face_tags, "north_top", &N4);
    int* ttag = mesh_create_tag(mesh->face_tags, top_tag, N1+N2+N3+N4);
    memcpy(&ttag[0], top1, sizeof(int) * N1);
    memcpy(&ttag[N1], top2, sizeof(int) * N2);
    memcpy(&ttag[N2], top3, sizeof(int) * N3);
    memcpy(&ttag[N3], top4, sizeof(int) * N4);
  }

  // Clean up.
  for (int i = 0; i < 4; ++i)
    mesh_free(radial_blocks[i]);

  // Now partition the thing if we've been asked to.
  if (comm != MPI_COMM_SELF)
    partition_mesh(&mesh, comm, NULL, 0.0);

  return mesh;
}

