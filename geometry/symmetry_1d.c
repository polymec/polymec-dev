// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "geometry/symmetry.h"
#include "geometry/create_rectilinear_mesh.h"
#include "geometry/create_uniform_mesh.h"
#include "geometry/cubic_lattice.h"

// Symmetry-related features.
const char* SYMMETRIC = "symmetric";
const char* ONE_DIMENSIONAL = "one dimensional";
const char* TWO_DIMENSIONAL = "two dimensional";
const char* CARTESIAN_1D = "Cartesian 1D";
const char* CARTESIAN_2D = "Cartesian 2D";
const char* CYLINDRICAL_1D = "cylindrical 1D";
const char* CYLINDRICAL_RZ = "cylindrical r-z";
const char* SPHERICAL_1D = "spherical 1D";

// This helper function computes the logarithmic spacing for N points spanning
// the interval [x1, x2] with the given log factor, placing those points into 
// the array points.
static void compute_log_spacing(real_t x1, real_t x2, real_t log_factor, int N, real_t* points)
{
  // Figure out the innermost cell length.
  real_t sum = 0.0;
  for (int i = 0; i < N; ++i)
    sum += pow(log_factor, 1.0*i);
  real_t dx0 = (x2 - x1) / sum;

  // Now compute the rest of the spacings.
  real_t last_dx = dx0;
  for (int i = 0; i < N-1; ++i)
  {
    real_t dx = (i < N-1) ? pow(log_factor, 1.0*i) * dx0
      : last_dx;
    points[i] = (i == 0) ? dx0 : points[i-1] + dx;
    last_dx = dx;
  }
}

mesh_t* create_uniform_cartesian_1d_mesh(MPI_Comm comm, real_t x1, real_t x2, int N)
{
  ASSERT(x2 > x1);
  ASSERT(N > 0);

  // This is really just a nonuniform Cartesian mesh with, um, uniform spacing.
  real_t* xs = polymec_malloc(sizeof(real_t) * (N+1));
  xs[0] = x1;
  real_t dx = (x2 - x1) / N;
  for (int i = 1; i < N+1; ++i)
    xs[i] = xs[i-1] + dx;
  mesh_t* mesh = create_nonuniform_cartesian_1d_mesh(comm, xs, N);
  polymec_free(xs);

  return mesh;
}

mesh_t* create_logarithmic_cartesian_1d_mesh(MPI_Comm comm, real_t x1, real_t x2, real_t log_factor, int N)
{
  // This is really just a nonuniform Cartesian mesh with logarithmic spacings.
  real_t* xs = polymec_malloc(sizeof(real_t) * (N+1));
  compute_log_spacing(x1, x2, log_factor, N+1, xs);
  mesh_t* mesh = create_nonuniform_cartesian_1d_mesh(comm, xs, N);
  polymec_free(xs);

  return mesh;
}

mesh_t* create_nonuniform_cartesian_1d_mesh(MPI_Comm comm, real_t* xs, int N)
{
  // This is really just a 1D rectilinear mesh with the given spacings.
  real_t ys[2] = {-0.5, 0.5};
  real_t zs[2] = {-0.5, 0.5};
  mesh_t* mesh = create_rectilinear_mesh(comm, xs, N+1, ys, 2, zs, 2);

  // Tag the boundary faces.
  tag_rectilinear_mesh_faces(mesh, "x1", "x2", "y1", "y2", "z1", "z2");

  // Add some symmetry-related features.
  mesh_add_feature(mesh, SYMMETRIC);
  mesh_add_feature(mesh, ONE_DIMENSIONAL);
  mesh_add_feature(mesh, CARTESIAN_1D);

  return mesh;
}

mesh_t* create_uniform_cylindrical_1d_mesh(MPI_Comm comm, real_t r1, real_t r2, int N)
{
  ASSERT(r1 >= 0.0);
  ASSERT(r2 > r1);
  ASSERT(N > 0);

  // This is really just a nonuniform cylindrical mesh with, um, uniform spacing.
  real_t* rs = polymec_malloc(sizeof(real_t) * (N+1));
  rs[0] = r1;
  real_t dr = (r2 - r1) / N;
  for (int i = 1; i < N+1; ++i)
    rs[i] = rs[i-1] + dr;
  mesh_t* mesh = create_nonuniform_cylindrical_1d_mesh(comm, rs, N);
  polymec_free(rs);

  return mesh;
}

mesh_t* create_logarithmic_cylindrical_1d_mesh(MPI_Comm comm, real_t r1, real_t r2, real_t log_factor, int N)
{
  ASSERT(r1 >= 0.0);
  ASSERT(r2 > r1);
  ASSERT(N > 0);

  // This is really just a nonuniform cylindrical mesh with logarithmic spacings.
  real_t* rs = polymec_malloc(sizeof(real_t) * (N+1));
  compute_log_spacing(r1, r2, log_factor, N+1, rs);
  mesh_t* mesh = create_nonuniform_cylindrical_1d_mesh(comm, rs, N);
  polymec_free(rs);

  return mesh;
}

mesh_t* create_nonuniform_cylindrical_1d_mesh(MPI_Comm comm, real_t* rs, int N)
{
  ASSERT(N > 0);

#ifndef NDEBUG
  // Make sure the radial spacings make sense.
  ASSERT(rs[0] >= 0.0);
  for (int i = 1; i < N+1; ++i)
  {
    ASSERT(rs[i] > rs[i-1]);
  }
#endif

  mesh_t* mesh = NULL;
  cubic_lattice_t* lattice = NULL;

  // Compute the average radial spacing and use that as the axial spacing.
  real_t avg_dr = 0.0;
  for (int i = 1; i < N+1; ++i)
    avg_dr += rs[i] - rs[i-1];
  avg_dr /= N;
  real_t dz = avg_dr;

  if (reals_equal(rs[0], 0.0))
  {
    // If the inner radius (rs[0]) is zero, the innermost cell will be a wedge 
    // whose edge lies at the origin. We create a mesh of N-1 regular 
    // (hexahedral) cells to begin with.
    bbox_t bbox = {.x1 = rs[1], .x2 = rs[N-1], 
                   .y1 = -1.0, .y2 = 1.0, 
                   .z1 = -0.5*dz, .z2 = 0.5*dz}; // Wrong, but it doesn't matter.
    mesh = create_uniform_mesh(comm, N-1, 1, 1, &bbox);

    // Retrieve the (topological) cubic lattice for the mesh.
    lattice = mesh_property(mesh, "lattice");
    ASSERT(lattice != NULL);

    // Now we add the inner most cell, which is a wedge with 5 faces. The 
    // order of the faces is +x,-y,+y,-z,+z.
    mesh->num_cells += 1;
    mesh->cell_face_offsets = polymec_realloc(mesh->cell_face_offsets, 
                                              sizeof(int)*(mesh->num_cells+1));
    mesh->cell_face_offsets[mesh->num_cells] = mesh->cell_face_offsets[mesh->num_cells-1] + 5;
    int total_num_cell_faces = mesh->cell_face_offsets[mesh->num_cells];
    mesh->cell_faces = polymec_realloc(mesh->cell_faces, 
                                       sizeof(int)*total_num_cell_faces);
    mesh->cell_faces[total_num_cell_faces-5] = mesh->cell_faces[0]; // +x
    mesh->cell_faces[total_num_cell_faces-4] = mesh->num_faces;     // -y
    mesh->cell_faces[total_num_cell_faces-3] = mesh->num_faces + 1; // +y
    mesh->cell_faces[total_num_cell_faces-2] = mesh->num_faces + 2; // -z
    mesh->cell_faces[total_num_cell_faces-1] = mesh->num_faces + 3; // +z

    // Make sure we allocate the extra storage for cell volumes / centers.
    mesh->cell_volumes = polymec_realloc(mesh->cell_volumes, 
                                         sizeof(real_t)*mesh->num_cells);
    mesh->cell_centers = polymec_realloc(mesh->cell_centers, 
                                         sizeof(point_t)*mesh->num_cells);

    // Now add the 4 new faces, which introduce 2 new nodes and 5 new edges.
    // The -y/+y faces have 4 nodes/edges, and the -z/+z faces have 3.
    mesh->num_faces += 4;
    mesh->face_node_offsets = polymec_realloc(mesh->face_node_offsets, 
                                              sizeof(int)*(mesh->num_faces+1));
    mesh->face_node_offsets[mesh->num_faces-3] = mesh->face_node_offsets[mesh->num_faces-4] + 4; // -y
    mesh->face_node_offsets[mesh->num_faces-2] = mesh->face_node_offsets[mesh->num_faces-3] + 4; // +y
    mesh->face_node_offsets[mesh->num_faces-1] = mesh->face_node_offsets[mesh->num_faces-2] + 3; // -z
    mesh->face_node_offsets[mesh->num_faces  ] = mesh->face_node_offsets[mesh->num_faces-1] + 3; // +z
    int total_num_face_nodes = mesh->face_node_offsets[mesh->num_faces];
    mesh->face_nodes = polymec_realloc(mesh->face_nodes, 
                                       sizeof(int)*total_num_face_nodes);

    // We define the node orderings on these new faces such that they produce
    // outward normals when traversed according to the right hand rule.
    mesh->face_nodes[total_num_face_nodes-14] = mesh->num_nodes;                      // -y face node 0 (bottom axial node -- new)
    mesh->face_nodes[total_num_face_nodes-13] = (int)cubic_lattice_node(lattice, 0, 0, 0); // -y face node 1
    mesh->face_nodes[total_num_face_nodes-12] = (int)cubic_lattice_node(lattice, 0, 0, 1); // -y face node 2
    mesh->face_nodes[total_num_face_nodes-11] = mesh->num_nodes + 1;                  // -y face node 3 (top axial node -- new)

    mesh->face_nodes[total_num_face_nodes-10] = (int)cubic_lattice_node(lattice, 0, 1, 0); // +y face node 0
    mesh->face_nodes[total_num_face_nodes- 9] = mesh->num_nodes;                      // +y face node 1 (bottom axial node -- new)
    mesh->face_nodes[total_num_face_nodes- 8] = mesh->num_nodes + 1;                  // +y face node 2 (top axial node -- new)
    mesh->face_nodes[total_num_face_nodes- 7] = (int)cubic_lattice_node(lattice, 0, 1, 1); // +y face node 3

    mesh->face_nodes[total_num_face_nodes- 6] = mesh->num_nodes;                      // -z face node 0 (bottom axial node -- new)
    mesh->face_nodes[total_num_face_nodes- 5] = (int)cubic_lattice_node(lattice, 0, 1, 0); // -z face node 1 
    mesh->face_nodes[total_num_face_nodes- 4] = (int)cubic_lattice_node(lattice, 0, 0, 0); // -z face node 2

    mesh->face_nodes[total_num_face_nodes- 3] = mesh->num_nodes + 1;                  // +z face node 0 (top axial node -- new)
    mesh->face_nodes[total_num_face_nodes- 2] = (int)cubic_lattice_node(lattice, 0, 1, 1); // +z face node 1 
    mesh->face_nodes[total_num_face_nodes- 1] = (int)cubic_lattice_node(lattice, 0, 0, 1); // +z face node 2

    // Now for the face->edge connectivity.
    mesh->face_edge_offsets = polymec_realloc(mesh->face_edge_offsets, 
                                              sizeof(int)*(mesh->num_faces+1));
    mesh->face_edge_offsets[mesh->num_faces-3] = mesh->face_edge_offsets[mesh->num_faces-4] + 4; // -y
    mesh->face_edge_offsets[mesh->num_faces-2] = mesh->face_edge_offsets[mesh->num_faces-3] + 4; // +y
    mesh->face_edge_offsets[mesh->num_faces-1] = mesh->face_edge_offsets[mesh->num_faces-2] + 3; // -z
    mesh->face_edge_offsets[mesh->num_faces  ] = mesh->face_edge_offsets[mesh->num_faces-1] + 3; // +z
    int total_num_face_edges = mesh->face_edge_offsets[mesh->num_faces];
    mesh->face_edges = polymec_realloc(mesh->face_edges, 
                                       sizeof(int)*total_num_face_edges);

    mesh->face_edges[total_num_face_edges-14] = mesh->num_edges;                        // -y face edge 0 (axial edge -- new)
    mesh->face_edges[total_num_face_edges-13] = mesh->num_edges + 1;                    // -y face edge 1 (bottom -y new edge)
    mesh->face_edges[total_num_face_edges-12] = (int)cubic_lattice_z_edge(lattice, 0, 0, 0); // -y face edge 2 
    mesh->face_edges[total_num_face_edges-11] = mesh->num_edges + 2;                    // -y face edge 3 (top -y new edge)

    mesh->face_edges[total_num_face_edges-10] = (int)cubic_lattice_z_edge(lattice, 0, 1, 0); // +y face edge 0 
    mesh->face_edges[total_num_face_edges- 9] = mesh->num_edges + 3;                    // +y face edge 1 (bottom +y new edge)
    mesh->face_edges[total_num_face_edges- 8] = mesh->num_edges;                        // +y face edge 2 (axial edge -- new)
    mesh->face_edges[total_num_face_edges- 7] = mesh->num_edges + 4;                    // +y face edge 3 (top +y new edge)

    mesh->face_edges[total_num_face_edges- 6] = mesh->num_edges + 3;                    // -z face edge 1 (bottom +y new edge)
    mesh->face_edges[total_num_face_edges- 5] = (int)cubic_lattice_y_edge(lattice, 0, 0, 0); // -z face edge 2
    mesh->face_edges[total_num_face_edges- 4] = mesh->num_edges + 1;                    // -z face edge 3 (bottom -y new edge)

    mesh->face_edges[total_num_face_edges- 3] = mesh->num_edges + 2;                    // +z face edge 1 (top -y new edge)
    mesh->face_edges[total_num_face_edges- 2] = (int)cubic_lattice_y_edge(lattice, 0, 0, 1); // +z face edge 2
    mesh->face_edges[total_num_face_edges- 1] = mesh->num_edges + 4;                    // +z face edge 3 (top +y new edge)

    // Face->cell connectivity.
    mesh->face_cells = polymec_realloc(mesh->face_cells,
                                          sizeof(int)*2*mesh->num_faces);
    mesh->face_cells[1]                   = mesh->num_cells-1; // +x face
    mesh->face_cells[2*mesh->num_faces-8] = mesh->num_cells-1; // -y face
    mesh->face_cells[2*mesh->num_faces-7] = -1;
    mesh->face_cells[2*mesh->num_faces-6] = mesh->num_cells-1; // +y face
    mesh->face_cells[2*mesh->num_faces-5] = -1;
    mesh->face_cells[2*mesh->num_faces-4] = mesh->num_cells-1; // -z face
    mesh->face_cells[2*mesh->num_faces-3] = -1;
    mesh->face_cells[2*mesh->num_faces-2] = mesh->num_cells-1; // +z face
    mesh->face_cells[2*mesh->num_faces-1] = -1;

    // Make sure we allocate the extra storage for face centers/normals/areas.
    mesh->face_centers = polymec_realloc(mesh->face_centers, 
                                         sizeof(point_t)*mesh->num_faces);
    mesh->face_normals = polymec_realloc(mesh->face_normals, 
                                         sizeof(vector_t)*mesh->num_faces);
    mesh->face_areas = polymec_realloc(mesh->face_areas, 
                                       sizeof(real_t)*mesh->num_faces);

    // Add the two new nodes.
    mesh->num_nodes += 2;
    mesh->nodes = polymec_realloc(mesh->nodes, 
                                  sizeof(point_t)*mesh->num_nodes);

    // Set their coordinates.
    mesh->nodes[mesh->num_nodes-2].x = 0.0;
    mesh->nodes[mesh->num_nodes-2].y = 0.0;
    mesh->nodes[mesh->num_nodes-2].z = -0.5*dz;

    mesh->nodes[mesh->num_nodes-1].x = 0.0;
    mesh->nodes[mesh->num_nodes-1].y = 0.0;
    mesh->nodes[mesh->num_nodes-1].z = 0.5*dz;

    // Add the five new edges.
    mesh->num_edges += 5;
    mesh->edge_nodes = polymec_realloc(mesh->edge_nodes, 
                                       sizeof(int)*2*(mesh->num_edges));

    // Axial edge.
    mesh->edge_nodes[2*mesh->num_edges-10] = mesh->num_nodes - 2; 
    mesh->edge_nodes[2*mesh->num_edges- 9] = mesh->num_nodes - 1; 

    // Bottom -y new edge.
    mesh->edge_nodes[2*mesh->num_edges- 8] = (int)cubic_lattice_node(lattice, 0, 0, 0);
    mesh->edge_nodes[2*mesh->num_edges- 7] = mesh->num_nodes - 2; 

    // Top -y new edge.
    mesh->edge_nodes[2*mesh->num_edges- 6] = (int)cubic_lattice_node(lattice, 0, 0, 1);
    mesh->edge_nodes[2*mesh->num_edges- 5] = mesh->num_nodes - 1; 

    // Bottom +y new edge.
    mesh->edge_nodes[2*mesh->num_edges- 4] = (int)cubic_lattice_node(lattice, 0, 1, 0);
    mesh->edge_nodes[2*mesh->num_edges- 3] = mesh->num_nodes - 2; 

    // Top +y new edge.
    mesh->edge_nodes[2*mesh->num_edges- 2] = (int)cubic_lattice_node(lattice, 0, 1, 1);
    mesh->edge_nodes[2*mesh->num_edges- 1] = mesh->num_nodes - 1; 
  }
  else
  {
    // Otherwise, the mesh will topologically be rectilinear and we only need 
    // to set the geometry.
    bbox_t bbox = {.x1 = rs[0], .x2 = rs[N-1], 
                   .y1 = -1.0, .y2 = 1.0, 
                   .z1 = -0.5*dz, .z2 = 0.5*dz}; 
    mesh = create_uniform_mesh(comm, N, 1, 1, &bbox);

    // Retrieve the (topological) cubic lattice for the mesh.
    lattice = mesh_property(mesh, "lattice");
    ASSERT(lattice != NULL);
  }

  // Adjust the geometry of the mesh and compute volumes/areas.
  // For each of the non-inner cells, the x and y coordinates of 
  // each nodes should be set to its respective cylindrical position.

  // Choose the angle to make cells maximally isotropic. For isotropy,
  // sin(phi/2) = avg_dr/r.
  real_t dphi = 2.0*asin(avg_dr/rs[N-1]);
  real_t phi1 = -0.5*dphi, phi2 = 0.5*dphi;
  int Nrad = (reals_equal(rs[0], 0.0)) ? N : N+1;
  for (int i = 0; i < Nrad; ++i)
  {
    // The x coordinate already contains the radius.
    int n1 = (int)cubic_lattice_node(lattice, i, 0, 0);
    real_t r = rs[i];
    int n2 = (int)cubic_lattice_node(lattice, i, 1, 0);
    int n3 = (int)cubic_lattice_node(lattice, i, 0, 1);
    int n4 = (int)cubic_lattice_node(lattice, i, 1, 1);

    // Find the new x,y coordinates along the cylinder.
    mesh->nodes[n1].x = r * cos(phi1);
    mesh->nodes[n1].y = r * sin(phi1);

    mesh->nodes[n2].x = r * cos(phi2);
    mesh->nodes[n2].y = r * sin(phi2);

    mesh->nodes[n3].x = r * cos(phi1);
    mesh->nodes[n3].y = r * sin(phi1);

    mesh->nodes[n4].x = r * cos(phi2);
    mesh->nodes[n4].y = r * sin(phi2);
  }
  
  mesh_compute_geometry(mesh);

  // Add some symmetry-related features.
  mesh_add_feature(mesh, SYMMETRIC);
  mesh_add_feature(mesh, ONE_DIMENSIONAL);
  mesh_add_feature(mesh, CYLINDRICAL_1D);

  return mesh;
}

mesh_t* create_uniform_spherical_1d_mesh(MPI_Comm comm, real_t r1, real_t r2, int N)
{
  ASSERT(r1 >= 0.0);
  ASSERT(r2 > r1);
  ASSERT(N > 0);

  // This is really just a nonuniform spherical mesh with, um, uniform spacing.
  real_t* rs = polymec_malloc(sizeof(real_t) * (N+1));
  rs[0] = r1;
  real_t dr = (r2 - r1) / N;
  for (int i = 1; i < N+1; ++i)
    rs[i] = rs[i-1] + dr;
  mesh_t* mesh = create_nonuniform_spherical_1d_mesh(comm, rs, N);
  polymec_free(rs);

  return mesh;
}

mesh_t* create_logarithmic_spherical_1d_mesh(MPI_Comm comm, real_t r1, real_t r2, real_t log_factor, int N)
{
  ASSERT(r1 >= 0.0);
  ASSERT(r2 > r1);
  ASSERT(N > 0);

  // This is really just a nonuniform spherical mesh with logarithmic spacings.
  real_t* rs = polymec_malloc(sizeof(real_t) * (N+1));
  compute_log_spacing(r1, r2, log_factor, N+1, rs);
  mesh_t* mesh = create_nonuniform_spherical_1d_mesh(comm, rs, N);
  polymec_free(rs);

  return mesh;
}

mesh_t* create_nonuniform_spherical_1d_mesh(MPI_Comm comm, real_t* rs, int N)
{
  ASSERT(N > 0);

#ifndef NDEBUG
  // Make sure the radial spacings make sense.
  ASSERT(rs[0] >= 0.0);
  for (int i = 1; i < N+1; ++i)
  {
    ASSERT(rs[i] > rs[i-1]);
  }
#endif

  mesh_t* mesh = NULL;
  cubic_lattice_t* lattice = NULL;

  // Compute the average radial spacing and use that as the axial spacing.
  real_t avg_dr = 0.0;
  for (int i = 1; i < N+1; ++i)
    avg_dr += rs[i] - rs[i-1];
  avg_dr /= N;
  real_t dz = avg_dr;

  if (reals_equal(rs[0], 0.0))
  {
    // If the inner radius (rs[0]) is zero, the innermost cell will be a 
    // square pyramid whose tip lies at the origin. We create a mesh of N-1 
    // regular (hexahedral) cells to begin with.
    bbox_t bbox = {.x1 = rs[1], .x2 = rs[N-1], 
                   .y1 = -1.0, .y2 = 1.0, 
                   .z1 = -0.5*dz, .z2 = 0.5*dz}; // Wrong, but it doesn't matter.
    mesh = create_uniform_mesh(comm, N-1, 1, 1, &bbox);

    // Retrieve the (topological) cubic lattice for the mesh.
    lattice = mesh_property(mesh, "lattice");
    ASSERT(lattice != NULL);

    // Now we add the inner most cell, which is a pyramid with 5 faces. The 
    // order of the faces is +x,-y,+y,-z,+z.
    mesh->num_cells += 1;
    mesh->cell_face_offsets = polymec_realloc(mesh->cell_face_offsets, 
                                              sizeof(int)*(mesh->num_cells+1));
    mesh->cell_face_offsets[mesh->num_cells] = mesh->cell_face_offsets[mesh->num_cells-1] + 5;
    int total_num_cell_faces = mesh->cell_face_offsets[mesh->num_cells];
    mesh->cell_faces = polymec_realloc(mesh->cell_faces, 
                                       sizeof(int)*total_num_cell_faces);
    mesh->cell_faces[total_num_cell_faces-5] = mesh->cell_faces[0]; // +x
    mesh->cell_faces[total_num_cell_faces-4] = mesh->num_faces;     // -y
    mesh->cell_faces[total_num_cell_faces-3] = mesh->num_faces + 1; // +y
    mesh->cell_faces[total_num_cell_faces-2] = mesh->num_faces + 2; // -z
    mesh->cell_faces[total_num_cell_faces-1] = mesh->num_faces + 3; // +z

    // Make sure we allocate the extra storage for cell volumes / centers.
    mesh->cell_volumes = polymec_realloc(mesh->cell_volumes, 
                                         sizeof(real_t)*mesh->num_cells);
    mesh->cell_centers = polymec_realloc(mesh->cell_centers, 
                                         sizeof(point_t)*mesh->num_cells);

    // Now add the 4 new faces, which introduce 1 new node and 4 new edges.
    // The new faces all have 3 nodes/edges.
    mesh->num_faces += 4;
    mesh->face_node_offsets = polymec_realloc(mesh->face_node_offsets, 
                                              sizeof(int)*(mesh->num_faces+1));
    mesh->face_node_offsets[mesh->num_faces-3] = mesh->face_node_offsets[mesh->num_faces-4] + 3; // -y
    mesh->face_node_offsets[mesh->num_faces-2] = mesh->face_node_offsets[mesh->num_faces-3] + 3; // +y
    mesh->face_node_offsets[mesh->num_faces-1] = mesh->face_node_offsets[mesh->num_faces-2] + 3; // -z
    mesh->face_node_offsets[mesh->num_faces  ] = mesh->face_node_offsets[mesh->num_faces-1] + 3; // +z
    int total_num_face_nodes = mesh->face_node_offsets[mesh->num_faces];
    mesh->face_nodes = polymec_realloc(mesh->face_nodes, 
                                       sizeof(int)*total_num_face_nodes);

    // We define the node orderings on these new faces such that they produce
    // outward normals when traversed according to the right hand rule.
    mesh->face_nodes[total_num_face_nodes-12] = mesh->num_nodes;                      // -y face node 0 (origin node -- new)
    mesh->face_nodes[total_num_face_nodes-11] = (int)cubic_lattice_node(lattice, 0, 0, 0); // -y face node 1
    mesh->face_nodes[total_num_face_nodes-10] = (int)cubic_lattice_node(lattice, 0, 0, 1); // -y face node 2

    mesh->face_nodes[total_num_face_nodes- 9] = (int)cubic_lattice_node(lattice, 0, 1, 0); // +y face node 0
    mesh->face_nodes[total_num_face_nodes- 8] = mesh->num_nodes;                      // +y face node 1 (origin node -- new)
    mesh->face_nodes[total_num_face_nodes- 7] = (int)cubic_lattice_node(lattice, 0, 1, 1); // +y face node 3

    mesh->face_nodes[total_num_face_nodes- 6] = mesh->num_nodes;                      // -z face node 0 (origin node -- new)
    mesh->face_nodes[total_num_face_nodes- 5] = (int)cubic_lattice_node(lattice, 0, 1, 0); // -z face node 1 
    mesh->face_nodes[total_num_face_nodes- 4] = (int)cubic_lattice_node(lattice, 0, 0, 0); // -z face node 2

    mesh->face_nodes[total_num_face_nodes- 3] = mesh->num_nodes;                      // +z face node 0 (origin node -- new)
    mesh->face_nodes[total_num_face_nodes- 2] = (int)cubic_lattice_node(lattice, 0, 1, 1); // +z face node 1 
    mesh->face_nodes[total_num_face_nodes- 1] = (int)cubic_lattice_node(lattice, 0, 0, 1); // +z face node 2

    // Now for the face->edge connectivity.
    mesh->face_edge_offsets = polymec_realloc(mesh->face_edge_offsets, 
                                              sizeof(int)*(mesh->num_faces+1));
    mesh->face_edge_offsets[mesh->num_faces-3] = mesh->face_edge_offsets[mesh->num_faces-4] + 3; // -y
    mesh->face_edge_offsets[mesh->num_faces-2] = mesh->face_edge_offsets[mesh->num_faces-3] + 3; // +y
    mesh->face_edge_offsets[mesh->num_faces-1] = mesh->face_edge_offsets[mesh->num_faces-2] + 3; // -z
    mesh->face_edge_offsets[mesh->num_faces  ] = mesh->face_edge_offsets[mesh->num_faces-1] + 3; // +z
    int total_num_face_edges = mesh->face_edge_offsets[mesh->num_faces];
    mesh->face_edges = polymec_realloc(mesh->face_edges, 
                                       sizeof(int)*total_num_face_edges);

    mesh->face_edges[total_num_face_edges-12] = mesh->num_edges;                        // -y face edge 1 (bottom -y new edge)
    mesh->face_edges[total_num_face_edges-11] = (int)cubic_lattice_z_edge(lattice, 0, 0, 0); // -y face edge 2 
    mesh->face_edges[total_num_face_edges-10] = mesh->num_edges + 1;                    // -y face edge 3 (top -y new edge)

    mesh->face_edges[total_num_face_edges- 9] = (int)cubic_lattice_z_edge(lattice, 0, 1, 0); // +y face edge 0 
    mesh->face_edges[total_num_face_edges- 8] = mesh->num_edges + 2;                    // +y face edge 1 (bottom +y new edge)
    mesh->face_edges[total_num_face_edges- 7] = mesh->num_edges + 3;                    // +y face edge 3 (top +y new edge)

    mesh->face_edges[total_num_face_edges- 6] = mesh->num_edges + 2;                    // -z face edge 1 (bottom +y new edge)
    mesh->face_edges[total_num_face_edges- 5] = (int)cubic_lattice_y_edge(lattice, 0, 0, 0); // -z face edge 2
    mesh->face_edges[total_num_face_edges- 4] = mesh->num_edges;                        // -z face edge 3 (bottom -y new edge)

    mesh->face_edges[total_num_face_edges- 3] = mesh->num_edges + 1;                    // +z face edge 1 (top -y new edge)
    mesh->face_edges[total_num_face_edges- 2] = (int)cubic_lattice_y_edge(lattice, 0, 0, 1); // +z face edge 2
    mesh->face_edges[total_num_face_edges- 1] = mesh->num_edges + 3;                    // +z face edge 3 (top +y new edge)

    // Face->cell connectivity.
    mesh->face_cells = polymec_realloc(mesh->face_cells,
                                       sizeof(int)*2*mesh->num_faces);
    mesh->face_cells[1]                   = mesh->num_cells-1; // +x face
    mesh->face_cells[2*mesh->num_faces-8] = mesh->num_cells-1; // -y face
    mesh->face_cells[2*mesh->num_faces-7] = -1;
    mesh->face_cells[2*mesh->num_faces-6] = mesh->num_cells-1; // +y face
    mesh->face_cells[2*mesh->num_faces-5] = -1;
    mesh->face_cells[2*mesh->num_faces-4] = mesh->num_cells-1; // -z face
    mesh->face_cells[2*mesh->num_faces-3] = -1;
    mesh->face_cells[2*mesh->num_faces-2] = mesh->num_cells-1; // +z face
    mesh->face_cells[2*mesh->num_faces-1] = -1;

    // Make sure we allocate the extra storage for face centers/normals/areas.
    mesh->face_centers = polymec_realloc(mesh->face_centers, 
                                         sizeof(point_t)*mesh->num_faces);
    mesh->face_normals = polymec_realloc(mesh->face_normals, 
                                         sizeof(vector_t)*mesh->num_faces);
    mesh->face_areas = polymec_realloc(mesh->face_areas, 
                                       sizeof(real_t)*mesh->num_faces);

    // Add the new node at the origin.
    mesh->num_nodes += 1;
    mesh->nodes = polymec_realloc(mesh->nodes, sizeof(point_t)*mesh->num_nodes);

    mesh->nodes[mesh->num_nodes-2].x = 0.0;
    mesh->nodes[mesh->num_nodes-2].y = 0.0;
    mesh->nodes[mesh->num_nodes-2].z = 0.0;

    // Add the four new edges.
    mesh->num_edges += 4;
    mesh->edge_nodes = polymec_realloc(mesh->edge_nodes, sizeof(int)*2*(mesh->num_edges));

    // Bottom -y new edge.
    mesh->edge_nodes[2*mesh->num_edges- 8] = (int)cubic_lattice_node(lattice, 0, 0, 0);
    mesh->edge_nodes[2*mesh->num_edges- 7] = mesh->num_nodes - 2; 

    // Top -y new edge.
    mesh->edge_nodes[2*mesh->num_edges- 6] = (int)cubic_lattice_node(lattice, 0, 0, 1);
    mesh->edge_nodes[2*mesh->num_edges- 5] = mesh->num_nodes - 1; 

    // Bottom +y new edge.
    mesh->edge_nodes[2*mesh->num_edges- 4] = (int)cubic_lattice_node(lattice, 0, 1, 0);
    mesh->edge_nodes[2*mesh->num_edges- 3] = mesh->num_nodes - 2; 

    // Top +y new edge.
    mesh->edge_nodes[2*mesh->num_edges- 2] = (int)cubic_lattice_node(lattice, 0, 1, 1);
    mesh->edge_nodes[2*mesh->num_edges- 1] = mesh->num_nodes - 1; 
  }
  else
  {
    // Otherwise, the mesh will topologically be rectilinear and we only need 
    // to set the geometry.
    bbox_t bbox = {.x1 = rs[0], .x2 = rs[N-1], 
                   .y1 = -1.0, .y2 = 1.0, 
                   .z1 = -0.5*dz, .z2 = 0.5*dz}; 
    mesh = create_uniform_mesh(comm, N, 1, 1, &bbox);

    // Retrieve the (topological) cubic lattice for the mesh.
    lattice = mesh_property(mesh, "lattice");
    ASSERT(lattice != NULL);
  }

  // Adjust the geometry of the mesh and compute volumes/areas.
  // For each of the non-inner cells, the x and y coordinates of 
  // each nodes should be set to its respective cylindrical position.

  // Choose the angle to make cells maximally isotropic. For isotropy,
  // dx ~ dy ~ dz, so
  // sin(dphi/2) = avg_dr/r,  (dy ~ dx) and 
  // cos(theta + pi/2) - cos(-theta + pi/2) = avg_dr/r, or
  // sin(theta) - sin(-theta) = avg_dr/r, or
  // 2*sin(theta) = avg_dr/r.
  real_t dphi = 2.0*asin(avg_dr/rs[N-1]);
  real_t phi1 = -0.5*dphi, phi2 = 0.5*dphi;
  real_t theta = 2.0*asin(0.5*avg_dr/rs[N-1]);
  real_t theta1 = 0.5*M_PI - theta, theta2 = 0.5*M_PI + theta;
  int Nrad = (reals_equal(rs[0], 0.0)) ? N : N+1;
  for (int i = 0; i < Nrad; ++i)
  {
    // The x coordinate already contains the radius.
    int n1 = (int)cubic_lattice_node(lattice, i, 0, 0);
    real_t r = rs[i];
    int n2 = (int)cubic_lattice_node(lattice, i, 1, 0);
    int n3 = (int)cubic_lattice_node(lattice, i, 0, 1);
    int n4 = (int)cubic_lattice_node(lattice, i, 1, 1);

    // Find the new x,y coordinates along the sphere.
    mesh->nodes[n1].x = r * sin(theta1) * cos(phi1);
    mesh->nodes[n1].y = r * sin(theta1) * sin(phi1);
    mesh->nodes[n1].z = r * cos(theta1);

    mesh->nodes[n2].x = r * sin(theta1) * cos(phi2);
    mesh->nodes[n2].y = r * sin(theta1) * sin(phi2);
    mesh->nodes[n2].z = r * cos(theta1);

    mesh->nodes[n3].x = r * sin(theta2) * cos(phi1);
    mesh->nodes[n3].y = r * sin(theta2) * sin(phi1);
    mesh->nodes[n3].z = r * cos(theta2);

    mesh->nodes[n4].x = r * sin(theta2) * cos(phi2);
    mesh->nodes[n4].y = r * sin(theta2) * sin(phi2);
    mesh->nodes[n4].z = r * cos(theta2);
  }
  
  mesh_compute_geometry(mesh);

  // Add some symmetry-related features.
  mesh_add_feature(mesh, SYMMETRIC);
  mesh_add_feature(mesh, ONE_DIMENSIONAL);
  mesh_add_feature(mesh, SPHERICAL_1D);

  return mesh;
}

