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

  mesh_t* mesh = NULL;

  if (rs[0] == 0.0)
  {
    // If the inner radius (rs[0]) is zero, the innermost cell will be a wedge 
    // whose edge lies at the origin. We create a mesh of N-1 regular 
    // (hexahedral) cells to begin with.
    bbox_t bbox = {.x1 = rs[1], .x2 = rs[N-1], 
                   .y1 = -1.0, .y2 = 1.0, 
                   .z1 = -1.0, .z2 = 1.0}; // Wrong, but it doesn't matter.
    create_uniform_mesh(comm, N-1, 1, 1, &bbox);

    // Retrieve the (topological) cubic lattice for the mesh.
    cubic_lattice_t* lattice = mesh_property(mesh, "lattice");
    ASSERT(lattice != NULL);

    // Now we add the inner most cell, which is a wedge with 5 faces. The 
    // order of the faces is +r,-y,+y,-z,+z.
    mesh->num_cells += 1;
    mesh->cell_face_offsets = ARENA_REALLOC(mesh->arena, 
                                            mesh->cell_face_offsets, 
                                            sizeof(int)*(mesh->num_cells+1), 0);
    mesh->cell_face_offsets[mesh->num_cells] = mesh->cell_face_offsets[mesh->num_cells-1] + 5;
    int total_num_cell_faces = mesh->cell_face_offsets[mesh->num_cells];
    mesh->cell_faces = ARENA_REALLOC(mesh->arena, 
                                     mesh->cell_faces, 
                                     sizeof(int)*total_num_cell_faces, 0);
    mesh->cell_faces[total_num_cell_faces-5] = mesh->cell_faces[0]; // +r
    mesh->cell_faces[total_num_cell_faces-4] = mesh->num_faces;     // -y
    mesh->cell_faces[total_num_cell_faces-3] = mesh->num_faces + 1; // +y
    mesh->cell_faces[total_num_cell_faces-2] = mesh->num_faces + 2; // -z
    mesh->cell_faces[total_num_cell_faces-1] = mesh->num_faces + 3; // +z

    // Now add the 4 new faces, which introduce 2 new nodes and 5 new edges.
    // The -y/+y faces have 4 nodes/edges, and the -z/+z faces have 3.
    mesh->num_faces += 4;
    mesh->face_node_offsets = ARENA_REALLOC(mesh->arena, 
                                            mesh->face_node_offsets, 
                                            sizeof(int)*(mesh->num_faces+1), 0);
    mesh->face_node_offsets[mesh->num_faces+1-3] = mesh->face_node_offsets[mesh->num_faces-4] + 4; // -y
    mesh->face_node_offsets[mesh->num_faces+1-2] = mesh->face_node_offsets[mesh->num_faces-3] + 4; // +y
    mesh->face_node_offsets[mesh->num_faces+1-1] = mesh->face_node_offsets[mesh->num_faces-2] + 3; // -z
    mesh->face_node_offsets[mesh->num_faces+1  ] = mesh->face_node_offsets[mesh->num_faces-1] + 3; // +z
    int total_num_face_nodes = mesh->face_node_offsets[mesh->num_faces];
    mesh->face_nodes = ARENA_REALLOC(mesh->arena, 
                                     mesh->face_nodes, 
                                     sizeof(int)*total_num_face_nodes, 0);

    // We define the node orderings on these new faces such that they produce
    // outward normals when traversed according to the right hand rule.
    mesh->face_nodes[total_num_face_nodes-14] = mesh->num_nodes;                      // -y face node 0 (bottom axial node -- new)
    mesh->face_nodes[total_num_face_nodes-13] = cubic_lattice_node(lattice, 0, 0, 0); // -y face node 1
    mesh->face_nodes[total_num_face_nodes-12] = cubic_lattice_node(lattice, 0, 0, 1); // -y face node 2
    mesh->face_nodes[total_num_face_nodes-11] = mesh->num_nodes + 1;                  // -y face node 3 (top axial node -- new)

    mesh->face_nodes[total_num_face_nodes-10] = cubic_lattice_node(lattice, 0, 1, 0); // +y face node 0
    mesh->face_nodes[total_num_face_nodes- 9] = mesh->num_nodes;                      // +y face node 1 (bottom axial node -- new)
    mesh->face_nodes[total_num_face_nodes- 8] = mesh->num_nodes + 1;                  // +y face node 1 (bottom axial node -- new)
    mesh->face_nodes[total_num_face_nodes- 7] = cubic_lattice_node(lattice, 0, 1, 1); // -y face node 1

    mesh->face_nodes[total_num_face_nodes- 6] = mesh->num_nodes;                      // -z face node 0 (bottom axial node -- new)
    mesh->face_nodes[total_num_face_nodes- 5] = cubic_lattice_node(lattice, 0, 1, 0); // -z face node 1 
    mesh->face_nodes[total_num_face_nodes- 4] = cubic_lattice_node(lattice, 0, 0, 0); // -z face node 2

    mesh->face_nodes[total_num_face_nodes- 3] = mesh->num_nodes + 1;                  // +z face node 0 (top axial node -- new)
    mesh->face_nodes[total_num_face_nodes- 2] = cubic_lattice_node(lattice, 0, 1, 1); // +z face node 1 
    mesh->face_nodes[total_num_face_nodes- 1] = cubic_lattice_node(lattice, 0, 0, 1); // +z face node 2

    // Now for the face->edge connectivity.
    mesh->face_edge_offsets = ARENA_REALLOC(mesh->arena, 
                                            mesh->face_edge_offsets, 
                                            sizeof(int)*(mesh->num_faces+1), 0);
    mesh->face_edge_offsets[mesh->num_faces+1-3] = mesh->face_edge_offsets[mesh->num_faces-4] + 4; // -y
    mesh->face_edge_offsets[mesh->num_faces+1-2] = mesh->face_edge_offsets[mesh->num_faces-3] + 4; // +y
    mesh->face_edge_offsets[mesh->num_faces+1-1] = mesh->face_edge_offsets[mesh->num_faces-2] + 3; // -z
    mesh->face_edge_offsets[mesh->num_faces+1  ] = mesh->face_edge_offsets[mesh->num_faces-1] + 3; // +z
    int total_num_face_edges = mesh->face_edge_offsets[mesh->num_faces];
    mesh->face_edges = ARENA_REALLOC(mesh->arena, 
                                     mesh->face_edges, 
                                     sizeof(int)*total_num_face_edges, 0);

    mesh->face_edges[total_num_face_edges-14] = mesh->num_edges;                        // -y face edge 0 (axial edge -- new)
    mesh->face_edges[total_num_face_edges-13] = mesh->num_edges + 1;                    // -y face edge 1 (bottom -y new edge)
    mesh->face_edges[total_num_face_edges-12] = cubic_lattice_z_edge(lattice, 0, 0, 0); // -y face edge 2 
    mesh->face_edges[total_num_face_edges-11] = mesh->num_edges + 2;                    // -y face edge 3 (top -y new edge)

    mesh->face_edges[total_num_face_edges-10] = cubic_lattice_z_edge(lattice, 0, 1, 0); // +y face edge 0 
    mesh->face_edges[total_num_face_edges- 9] = mesh->num_edges + 3;                    // +y face edge 1 (bottom +y new edge)
    mesh->face_edges[total_num_face_edges- 8] = mesh->num_edges;                        // +y face edge 2 (axial edge -- new)
    mesh->face_edges[total_num_face_edges- 7] = mesh->num_edges + 4;                    // +y face edge 3 (top +y new edge)

    mesh->face_edges[total_num_face_edges- 6] = mesh->num_edges + 3;                    // -z face edge 1 (bottom +y new edge)
    mesh->face_edges[total_num_face_edges- 5] = cubic_lattice_y_edge(lattice, 0, 0, 0); // -z face edge 2
    mesh->face_edges[total_num_face_edges- 4] = mesh->num_edges + 1;                    // -z face edge 3 (bottom -y new edge)

    mesh->face_edges[total_num_face_edges- 3] = mesh->num_edges + 2;                    // +z face edge 1 (top -y new edge)
    mesh->face_edges[total_num_face_edges- 2] = cubic_lattice_y_edge(lattice, 0, 0, 1); // +z face edge 2
    mesh->face_edges[total_num_face_edges- 1] = mesh->num_edges + 4;                    // +z face edge 3 (top +y new edge)

    // Face->cell connectivity.
    mesh->face_cells = ARENA_REALLOC(mesh->arena, 
                                     mesh->face_cells,
                                     sizeof(int)*2*mesh->num_faces, 0);
    mesh->face_cells[1]                   = mesh->num_cells-1; // +r face
    mesh->face_cells[2*mesh->num_faces-4] = mesh->num_cells-1; // -y face
    mesh->face_cells[2*mesh->num_faces-3] = mesh->num_cells-1; // +y face
    mesh->face_cells[2*mesh->num_faces-2] = mesh->num_cells-1; // -z face
    mesh->face_cells[2*mesh->num_faces-1] = mesh->num_cells-1; // +z face

    // Add the two new nodes.
    mesh->num_nodes += 2;
    mesh->nodes = ARENA_REALLOC(mesh->arena, 
                                mesh->nodes, 
                                sizeof(point_t)*mesh->num_nodes, 0);
    mesh->nodes[mesh->num_nodes-2].x = rs[0];
    mesh->nodes[mesh->num_nodes-2].y = -0.5;
    mesh->nodes[mesh->num_nodes-2].z = -0.5;
    mesh->nodes[mesh->num_nodes-1].x = rs[0];
    mesh->nodes[mesh->num_nodes-1].y = +0.5;
    mesh->nodes[mesh->num_nodes-1].z = +0.5;

    // Add the five new edges.
    mesh->num_edges += 5;
    mesh->edge_nodes = ARENA_REALLOC(mesh->arena, mesh->edge_nodes, 
                                     sizeof(int)*2*(mesh->num_edges), 0);

    // Axial edge.
    mesh->edge_nodes[2*mesh->num_edges-10] = mesh->num_nodes - 2; 
    mesh->edge_nodes[2*mesh->num_edges- 9] = mesh->num_nodes - 1; 

    // Bottom -y new edge.
    mesh->edge_nodes[2*mesh->num_edges- 8] = cubic_lattice_node(lattice, 0, 0, 0);
    mesh->edge_nodes[2*mesh->num_edges- 7] = mesh->num_nodes - 2; 

    // Top -y new edge.
    mesh->edge_nodes[2*mesh->num_edges- 6] = cubic_lattice_node(lattice, 0, 0, 1);
    mesh->edge_nodes[2*mesh->num_edges- 5] = mesh->num_nodes - 1; 

    // Bottom +y new edge.
    mesh->edge_nodes[2*mesh->num_edges- 4] = cubic_lattice_node(lattice, 0, 1, 0);
    mesh->edge_nodes[2*mesh->num_edges- 3] = mesh->num_nodes - 2; 

    // Top +y new edge.
    mesh->edge_nodes[2*mesh->num_edges- 2] = cubic_lattice_node(lattice, 0, 1, 1);
    mesh->edge_nodes[2*mesh->num_edges- 1] = mesh->num_nodes - 1; 
  }
  else
  {
    // Otherwise, the mesh will topologically be rectilinear and we only need 
    // to set the geometry.
    bbox_t bbox = {.x1 = rs[0], .x2 = rs[N-1], 
                   .y1 = -1.0, .y2 = 1.0, 
                   .z1 = -1.0, .z2 = 1.0}; 
    create_uniform_mesh(comm, N, 1, 1, &bbox);
  }

  // Adjust the geometry of the mesh and compute volumes/areas.
  // FIXME
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

