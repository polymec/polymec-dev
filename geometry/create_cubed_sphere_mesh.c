// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/partition_mesh.h"
#include "geometry/cubic_lattice.h"
#include "geometry/create_uniform_mesh.h"
#include "geometry/rotate_mesh.h"
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
  real_t X1 = X * sqrt(1.0 - 0.5*Y*Y - 0.5*Z*Z + one_third*Y*Y*Z*Z);
  real_t Y1 = Y * sqrt(1.0 - 0.5*Z*Z - 0.5*X*X + one_third*Z*Z*X*X);
  real_t Z1 = Z * sqrt(1.0 - 0.5*X*X - 0.5*Y*Y + one_third*X*X*Y*Y);

  // Finally, map to the circle with radius L/2.
  real_t r = sqrt(X1*X1 + Y1*Y1 + Z1*Z1);
  real_t theta = acos(Z1/(r + 1e-8));
  real_t phi = atan2(Y1, X1);

  y->x = 0.5*L * r * cos(phi) * sin(theta);
  y->y = 0.5*L * r * sin(phi) * sin(theta);
  y->z = 0.5*L * r * cos(theta);
}

// This renames the given boundary face tag on a block, also replacing any 
// such entry found in the rectilinear_boundary_tags property, which is 
// used by create_welded_block_mesh to identify seams.
static void rename_boundary_tag(mesh_t* block, 
                                const char* old_tag_name, 
                                const char* new_tag_name)
{
  // Rename the tag itself.
  mesh_rename_tag(block->face_tags, old_tag_name, new_tag_name);

  // Now find its entry in rectilinear_boundary_tags and replace it.
  string_array_t* btags = mesh_property(block, "rectilinear_boundary_tags");
  ASSERT(btags != NULL);
  for (int i = 0; i < btags->size; ++i)
  {
    if (strcmp(btags->data[i], old_tag_name) == 0)
    {
      // We simply replace the entry. The destructor for this entry doesn't 
      // change, so we allocate a new string with string_dup.
      string_free(btags->data[i]);
      btags->data[i] = string_dup(new_tag_name);
    }
  }
}

static void rename_surface_seam_tags(mesh_t** panels)
{
  rename_boundary_tag(panels[0], "panel_0_0", "seam_0_3");
  rename_boundary_tag(panels[0], "panel_0_1", "seam_0_1");
  rename_boundary_tag(panels[0], "panel_0_2", "seam_0_4");
  rename_boundary_tag(panels[0], "panel_0_3", "seam_0_5");

  rename_boundary_tag(panels[1], "panel_1_0", "seam_1_5");
  rename_boundary_tag(panels[1], "panel_1_1", "seam_1_4");
  rename_boundary_tag(panels[1], "panel_1_2", "seam_0_1");
  rename_boundary_tag(panels[1], "panel_1_3", "seam_1_2");

  rename_boundary_tag(panels[2], "panel_2_0", "seam_2_3");
  rename_boundary_tag(panels[2], "panel_2_1", "seam_1_2");
  rename_boundary_tag(panels[2], "panel_2_2", "seam_2_5");
  rename_boundary_tag(panels[2], "panel_2_3", "seam_2_4");

  rename_boundary_tag(panels[3], "panel_3_0", "seam_0_3");
  rename_boundary_tag(panels[3], "panel_3_1", "seam_2_3");
  rename_boundary_tag(panels[3], "panel_3_2", "seam_3_5");
  rename_boundary_tag(panels[3], "panel_3_3", "seam_3_4");

  rename_boundary_tag(panels[4], "panel_4_0", "seam_3_4");
  rename_boundary_tag(panels[4], "panel_4_1", "seam_1_4");
  rename_boundary_tag(panels[4], "panel_4_2", "seam_2_4");
  rename_boundary_tag(panels[4], "panel_4_3", "seam_0_4");

  rename_boundary_tag(panels[5], "panel_5_0", "seam_3_5");
  rename_boundary_tag(panels[5], "panel_5_1", "seam_1_5");
  rename_boundary_tag(panels[5], "panel_5_2", "seam_0_5");
  rename_boundary_tag(panels[5], "panel_5_3", "seam_3_5");
}

mesh_t* create_cubed_sphere_mesh(MPI_Comm comm,
                                 int ns, int nr, 
                                 real_t r, real_t R, 
                                 const char* R_tag)
{
  // Construct the center block.
  bbox_t bbox = {.x1 = -r, .x2 = r,
                 .y1 = -r, .y2 = r,
                 .z1 = -r, .z2 = r};
  mesh_t* center_block = create_uniform_mesh(MPI_COMM_SELF, ns, ns, ns, &bbox);
  tag_rectilinear_mesh_faces(center_block, 
                             "center_0", "center_1",
                             "center_2", "center_3",
                             "center_4", "center_5");

  // Deform the center block.
  cubic_lattice_t* lattice = cubic_lattice_new(ns, ns, ns);
  for (int i = 0; i <= ns; ++i)
  {
    for (int j = 0; j <= ns; ++j)
    {
      for (int k = 0; k <= ns; ++k)
      {
        int n = (int)cubic_lattice_node(lattice, i, j, k);
        point_t x = center_block->nodes[n];
        map_to_sphere(2.0*r, &x, &(center_block->nodes[n]));
      }
    }
  }
  mesh_compute_geometry(center_block);

  // Construct the panels. By default, they're all North-pole panels.
  mesh_t* panels[6];
  for (int i = 0; i < 6; ++i)
  {
    char tag_prefix[1024];
    snprintf(tag_prefix, 1024, "panel_%d", i);
    panels[i] = create_cubed_sphere_panel(MPI_COMM_SELF, ns, nr, r, R, tag_prefix);
  }

  // Rotate the panels so that they are properly oriented.
  // Panels 0-3 are equatorial panels, panel 4 is south, 5 is north.
  point_t x0 = {.x = 0.0, .y = 0.0, .z = 0.0};
  {
    vector_t omega = {.x = 1.0, .y = 0.0, .z = 0.0};
    rotate_mesh(panels[0], &x0, 0.5*M_PI, &omega);
    rotate_mesh(panels[2], &x0, -0.5*M_PI, &omega);
  }
  {
    vector_t omega = {.x = 0.0, .y = 1.0, .z = 0.0};
    rotate_mesh(panels[1], &x0, 0.5*M_PI, &omega);
    rotate_mesh(panels[3], &x0, -0.5*M_PI, &omega);
  }
  {
    vector_t omega = {.x = 1.0, .y = 0.0};
    rotate_mesh(panels[4], &x0, M_PI, &omega);
  }

  for (int i = 0; i < 6; ++i)
    mesh_compute_geometry(panels[i]);

  // Rename the surface panel seam tags so that they can be welded.
  rename_surface_seam_tags(panels);

  // Weld'em panels.
  mesh_t* blocks[7] = {center_block, panels[0], panels[1], panels[2],
                       panels[3], panels[4], panels[5]};
  mesh_t* mesh = create_welded_block_mesh(blocks, 7, 1e-10);

#if 0
  // Set up the inner / outer radius tags.
  {
    int N1, N2, N3, N4, N5, N6;
    int* R1 = mesh_tag(panels[0]->face_tags, "panel_0_6", &N1);
    int* R2 = mesh_tag(panels[1]->face_tags, "panel_1_6", &N2);
    int* R3 = mesh_tag(panels[2]->face_tags, "panel_2_6", &N3);
    int* R4 = mesh_tag(panels[3]->face_tags, "panel_3_6", &N4);
    int* R5 = mesh_tag(panels[4]->face_tags, "panel_4_6", &N5);
    int* R6 = mesh_tag(panels[5]->face_tags, "panel_5_6", &N6);
    int* rtag = mesh_create_tag(mesh->face_tags, R_tag, N1+N2+N3+N4+N5+N6);
    memcpy(&rtag[0], outer1, sizeof(int) * N1);
    memcpy(&rtag[N1], outer2, sizeof(int) * N2);
    memcpy(&rtag[N2], outer3, sizeof(int) * N3);
    memcpy(&rtag[N3], outer4, sizeof(int) * N4);
  }
#endif
  // Clean up.
  for (int i = 0; i < 6; ++i)
    mesh_free(panels[i]);

  // Now partition the thing if we've been asked to.
  if (comm != MPI_COMM_SELF)
    partition_mesh(&mesh, comm, NULL, 0.0);

  return mesh;
}

mesh_t* create_cubed_spherical_shell_mesh(MPI_Comm comm,
                                          int ns, int nr,
                                          real_t r, real_t R,
                                          const char* r_tag,
                                          const char* R_tag)
{
  // Construct the panels. By default, they're all North-pole panels.
  mesh_t* panels[6];
  for (int i = 0; i < 6; ++i)
  {
    char tag_prefix[1024];
    snprintf(tag_prefix, 1024, "panel_%d", i);
    panels[i] = create_cubed_sphere_panel(MPI_COMM_SELF, ns, nr, r, R, tag_prefix);
  }

  // Rotate the panels so that they are properly oriented.
  // Panels 0-3 are equatorial panels, panel 4 is south, 5 is north.
  point_t x0 = {.x = 0.0, .y = 0.0, .z = 0.0};
  {
    vector_t omega = {.x = 1.0, .y = 0.0, .z = 0.0};
    rotate_mesh(panels[0], &x0, 0.5*M_PI, &omega);
    rotate_mesh(panels[2], &x0, -0.5*M_PI, &omega);
  }
  {
    vector_t omega = {.x = 0.0, .y = 1.0, .z = 0.0};
    rotate_mesh(panels[1], &x0, 0.5*M_PI, &omega);
    rotate_mesh(panels[3], &x0, -0.5*M_PI, &omega);
  }
  {
    vector_t omega = {.x = 1.0, .y = 0.0};
    rotate_mesh(panels[4], &x0, M_PI, &omega);
  }

  for (int i = 0; i < 6; ++i)
    mesh_compute_geometry(panels[i]);

  // Rename the surface panel seam tags so that they can be welded.
  rename_surface_seam_tags(panels);

  // Weld'em panels.
  mesh_t* mesh = create_welded_block_mesh(panels, 6, 1e-10);

#if 0
  // Set up the inner / outer radius tags.
  {
    int N1, N2, N3, N4, N5, N6;
    int* R1 = mesh_tag(panels[0]->face_tags, "panel_0_6", &N1);
    int* R2 = mesh_tag(panels[1]->face_tags, "panel_1_6", &N2);
    int* R3 = mesh_tag(panels[2]->face_tags, "panel_2_6", &N3);
    int* R4 = mesh_tag(panels[3]->face_tags, "panel_3_6", &N4);
    int* R5 = mesh_tag(panels[4]->face_tags, "panel_4_6", &N5);
    int* R6 = mesh_tag(panels[5]->face_tags, "panel_5_6", &N6);
    int* rtag = mesh_create_tag(mesh->face_tags, R_tag, N1+N2+N3+N4+N5+N6);
    memcpy(&rtag[0], outer1, sizeof(int) * N1);
    memcpy(&rtag[N1], outer2, sizeof(int) * N2);
    memcpy(&rtag[N2], outer3, sizeof(int) * N3);
    memcpy(&rtag[N3], outer4, sizeof(int) * N4);
  }
#endif
  // Clean up.
  for (int i = 0; i < 6; ++i)
    mesh_free(panels[i]);

  // Now partition the thing if we've been asked to.
  if (comm != MPI_COMM_SELF)
    partition_mesh(&mesh, comm, NULL, 0.0);

  return mesh;
}

mesh_t* create_cubed_sphere_panel(MPI_Comm comm,
                                  int ns, int nr,
                                  real_t r, real_t R,
                                  const char* tag_prefix)
{
  ASSERT(ns >= 2);
  ASSERT(nr >= 2);
  ASSERT(r >= 0.0);
  ASSERT(R >= 0.0);
  ASSERT(r < R);

  // Generic bounding box -- doesn't really matter.
  bbox_t bbox = {.x1 = -0.5, .x2 = 0.5,
                 .y1 = -0.5, .y2 = 0.5,
                 .z1 = -0.5, .z2 = 0.5};

  // Indexing mechanism.
  cubic_lattice_t* lattice = cubic_lattice_new(ns, ns, nr);

  mesh_t* panel = create_uniform_mesh(MPI_COMM_SELF, ns, ns, nr, &bbox);

  // Tag the panel's faces accordingly.
  char tags[6][1024];
  for (int i = 0; i < 6; ++i)
    snprintf(tags[i], 1024, "%s_%d", tag_prefix, i);
  tag_rectilinear_mesh_faces(panel, tags[0], tags[1], tags[2], tags[3],
                             tags[4], tags[5]);

  // Grid spacings.
  real_t dr = (R - r) / (nr-1);
  real_t dx = 2.0 / ns;
  for (int k = 0; k <= nr; ++k)
  {
    // Use even radial spacing for now.
    real_t rk = r + k*dr;

    // Compute the latitude/longitude according to the Gnomonic cubed-sphere
    // projection for the North panel.
    for (int j = 0; j <= ns; ++j)
    {
      real_t yj = -1.0 + j*dx;
      for (int i = 0; i <= ns; ++i)
      {
        real_t xi = -1.0 + i*dx;
        real_t lon = atan2(xi, yj);
        real_t lat = atan2(1.0, sqrt(xi*xi + yj*yj));
        real_t theta = -lat + 0.5*M_PI, phi = lon; // spherical coordinates

        int n = (int)cubic_lattice_node(lattice, i, j, k);
        point_t xn = {.x = rk * cos(phi) * sin(theta),
                      .y = rk * sin(phi) * sin(theta),
                      .z = rk * cos(theta)};
        panel->nodes[n] = xn;
      }
    }
  }
  mesh_compute_geometry(panel);

  // Now partition the thing if we've been asked to.
  if (comm != MPI_COMM_SELF)
    partition_mesh(&panel, comm, NULL, 0.0);

  return panel;
}

