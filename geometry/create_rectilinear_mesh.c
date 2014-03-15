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

#include "core/unordered_set.h"
#include "geometry/create_rectilinear_mesh.h"
#include "geometry/cubic_lattice.h"

mesh_t* create_rectilinear_mesh(MPI_Comm comm, 
                                real_t* xs, int nxs, 
                                real_t* ys, int nys, 
                                real_t* zs, int nzs)
{
  ASSERT(nxs > 1);
  ASSERT(nys > 1);
  ASSERT(nzs > 1);

#ifndef NDEBUG
  for (int i = 1; i < nxs; ++i)
  {
    ASSERT(xs[i]-xs[i-1] > 0.0);
  }
  for (int j = 1; j < nys; ++j)
  {
    ASSERT(ys[j]-ys[j-1] > 0.0);
  }
  for (int k = 1; k < nzs; ++k)
  {
    ASSERT(zs[k]-zs[k-1] > 0.0);
  }
#endif

  // Numbers of cells in each direction.
  int nx = nxs-1, ny = nys-1, nz = nzs-1;

  // Create a cubic lattice object for indexing.
  cubic_lattice_t* lattice = cubic_lattice_new(nx, ny, nz);

//  // Compute the number of ghost cells.
//  int num_ghost_cells = 2*num_ghost*(ny*nz + nz*nx + nx*ny);

  // Create the mesh.
  // FIXME: Not parallel safe.
  mesh_t* mesh = mesh_new(comm, 
                          cubic_lattice_num_cells(lattice), 0,
                          cubic_lattice_num_faces(lattice),
                          cubic_lattice_num_nodes(lattice));
  mesh->cell_faces = ARENA_REALLOC(mesh->arena, mesh->cell_faces, sizeof(int)*6*mesh->num_cells, 0);
  mesh->face_nodes = ARENA_REALLOC(mesh->arena, mesh->face_nodes, sizeof(int)*4*mesh->num_faces, 0);
  mesh->cell_face_offsets[mesh->num_cells] = 6*mesh->num_cells;
  mesh->face_node_offsets[mesh->num_faces] = 4*mesh->num_faces;

  int_unordered_set_t* processed_nodes = int_unordered_set_new();
  for (int k = 0; k < nz; ++k)
  {
    for (int j = 0; j < ny; ++j)
    {
      for (int i = 0; i < nx; ++i)
      {
        // Hook up the cell and faces.
        int cell = cubic_lattice_cell(lattice, i, j, k);
        mesh->cell_face_offsets[cell] = 6*cell;

        // Note that we define the nodes of the faces to be ordered such that 
        // the right-hand rule produces a normal vector that always points 
        // in the positive x, y, or z direction. Thus, the faces whose normal 
        // vectors point in the opposite directions get their ones complement.
        mesh->cell_faces[6*cell+0] = ~cubic_lattice_x_face(lattice, i, j, k);
        mesh->cell_faces[6*cell+1] = cubic_lattice_x_face(lattice, i+1, j, k);
        mesh->cell_faces[6*cell+2] = ~cubic_lattice_y_face(lattice, i, j, k);
        mesh->cell_faces[6*cell+3] = cubic_lattice_y_face(lattice, i, j+1, k);
        mesh->cell_faces[6*cell+4] = ~cubic_lattice_z_face(lattice, i, j, k);
        mesh->cell_faces[6*cell+5] = cubic_lattice_z_face(lattice, i, j, k+1);

        // Hook up each face to its nodes.
        int nodes[6][4];

        // We use the reference cell below, which is typical of 3D
        // finite element schemes:
        //             
        //     7o----6o      z^  y
        //     /|    /|       | /
        //   4o----5o |       |/   x
        //    |3o---|2o       +---->
        //    |/    |/       
        //   0o----1o      
        //
        // The faces are numbered 0-5, with 0-1 being the (-/+) x faces,
        // 2-3 the (-/+) y faces, and 4-5 the (-/+) z faces.
        int node_indices[8];
        node_indices[0] = cubic_lattice_node(lattice, i, j, k);
        node_indices[1] = cubic_lattice_node(lattice, i+1, j, k);
        node_indices[2] = cubic_lattice_node(lattice, i+1, j+1, k);
        node_indices[3] = cubic_lattice_node(lattice, i, j+1, k);
        node_indices[4] = cubic_lattice_node(lattice, i, j, k+1);
        node_indices[5] = cubic_lattice_node(lattice, i+1, j, k+1);
        node_indices[6] = cubic_lattice_node(lattice, i+1, j+1, k+1);
        node_indices[7] = cubic_lattice_node(lattice, i, j+1, k+1);

        // Face 0 (-x) -- backward traversal
        nodes[0][0] = node_indices[7];
        nodes[0][1] = node_indices[4];
        nodes[0][2] = node_indices[0];
        nodes[0][3] = node_indices[3];

        // Face 1 (+x) -- forward traversal
        nodes[1][0] = node_indices[1];
        nodes[1][1] = node_indices[2];
        nodes[1][2] = node_indices[6];
        nodes[1][3] = node_indices[5];

        // Face 2 (-y) -- backward traversal
        nodes[2][0] = node_indices[4];
        nodes[2][1] = node_indices[5];
        nodes[2][2] = node_indices[1];
        nodes[2][3] = node_indices[0];

        // Face 3 (+y) -- forward traversal
        nodes[3][0] = node_indices[2];
        nodes[3][1] = node_indices[3];
        nodes[3][2] = node_indices[7];
        nodes[3][3] = node_indices[6];

        // Face 4 (-z) -- backward traversal
        nodes[4][0] = node_indices[0];
        nodes[4][1] = node_indices[1];
        nodes[4][2] = node_indices[2];
        nodes[4][3] = node_indices[3];

        // Face 5 (+z) -- forward traversal
        nodes[5][0] = node_indices[4];
        nodes[5][1] = node_indices[5];
        nodes[5][2] = node_indices[6];
        nodes[5][3] = node_indices[7];

        // Hook everything up.
        for (int f = 0; f < 6; ++f)
        {
          int face = mesh->cell_faces[6*cell+f];
          if (face < 0) face = ~face;
          if (mesh->face_cells[2*face] == -1)
            mesh->face_cells[2*face] = cell;
          else if (mesh->face_cells[2*face+1] == -1)
            mesh->face_cells[2*face+1] = cell;

          mesh->face_node_offsets[face] = 4*face;
          for (int n = 0; n < 4; ++n)
            mesh->face_nodes[4*face+n] = nodes[f][n];
        }

        // Assign the node positions for the cell.
        static const int i_offsets[] = {0, 1, 1, 0, 0, 1, 1, 0};
        static const int j_offsets[] = {0, 0, 1, 1, 0, 0, 1, 1};
        static const int k_offsets[] = {0, 0, 0, 0, 1, 1, 1, 1};
        for (int n = 0; n < 8; ++n)
        {
          if (!int_unordered_set_contains(processed_nodes, node_indices[n]))
          {
            point_t* node = &mesh->nodes[node_indices[n]];
            node->x = xs[i + i_offsets[n]];
            node->y = ys[j + j_offsets[n]];
            node->z = zs[k + k_offsets[n]];
            int_unordered_set_insert(processed_nodes, node_indices[n]);
          }
        }
      }
    }
  }


#if 0
  // Set up ghost cells.
  int gindex = mesh->num_cells;
  if (num_ghost > 0)
  {
    // x faces.
    for (int j = 0; j < ny; ++j)
    {
      for (int k = 0; k < nz; ++k)
      {
        // Hook up the ghost cells.
        int low_cell = cubic_lattice_cell(lattice, 0, j, k);
        cell_t* low_ghost_cell = &mesh->cells[gindex++];
        face_t* low_face = mesh->cells[low_cell].faces[0];
        ASSERT(low_face->cell2 == NULL); 
        low_face->cell2 = low_ghost_cell;
        mesh_attach_face_to_cell(mesh, low_face, low_ghost_cell);

        int high_cell = cubic_lattice_cell(lattice, nx-1, j, k);
        cell_t* high_ghost_cell = &mesh->cells[gindex++];
        face_t* high_face = mesh->cells[high_cell].faces[1];
        ASSERT(high_face->cell2 == NULL); 
        high_face->cell2 = high_ghost_cell;
        mesh_attach_face_to_cell(mesh, high_face, high_ghost_cell);

        // FIXME: Do ghost cells need node/edge connectivity?
      }
    }

    // y faces.
    for (int k = 0; k < nz; ++k)
    {
      for (int i = 0; i < nx; ++i)
      {
        // Hook up the ghost cells.
        int low_cell = cubic_lattice_cell(lattice, i, 0, k);
        cell_t* low_ghost_cell = &mesh->cells[gindex++];
        face_t* low_face = mesh->cells[low_cell].faces[2];
        ASSERT(low_face->cell2 == NULL); 
        low_face->cell2 = low_ghost_cell;
        mesh_attach_face_to_cell(mesh, low_face, low_ghost_cell);

        int high_cell = cubic_lattice_cell(lattice, i, ny-1, k);
        cell_t* high_ghost_cell = &mesh->cells[gindex++];
        face_t* high_face = mesh->cells[high_cell].faces[3];
        ASSERT(high_face->cell2 == NULL); 
        high_face->cell2 = high_ghost_cell;
        mesh_attach_face_to_cell(mesh, high_face, high_ghost_cell);

        // FIXME: Do ghost cells need node/edge connectivity?
      }
    }

    // z faces.
    for (int i = 0; i < nx; ++i)
    {
      for (int j = 0; j < ny; ++j)
      {
        // Hook up the ghost cells.
        int low_cell = cubic_lattice_cell(lattice, i, j, 0);
        cell_t* low_ghost_cell = &mesh->cells[gindex++];
        face_t* low_face = mesh->cells[low_cell].faces[4];
        ASSERT(low_face->cell2 == NULL); 
        low_face->cell2 = low_ghost_cell;
        mesh_attach_face_to_cell(mesh, low_face, low_ghost_cell);

        int high_cell = cubic_lattice_cell(lattice, i, j, nz-1);
        cell_t* high_ghost_cell = &mesh->cells[gindex++];
        face_t* high_face = mesh->cells[high_cell].faces[5];
        ASSERT(high_face->cell2 == NULL); 
        high_face->cell2 = high_ghost_cell;
        mesh_attach_face_to_cell(mesh, high_face, high_ghost_cell);

        // FIXME: Do ghost cells need node/edge connectivity?
      }
    }
  }
  ASSERT(gindex == (mesh->num_cells + mesh->num_ghost_cells));
#endif

  // Construct edge information.
  mesh_construct_edges(mesh);

  // Compute mesh geometry.
  mesh_compute_geometry(mesh);

  // Stash the lattice in the "lattice" property.
  mesh_set_property(mesh, "lattice", (void*)lattice, NULL);

  // Clean up.
  int_unordered_set_free(processed_nodes);

  return mesh;
}

void tag_rectilinear_mesh_faces(mesh_t* mesh, 
                                  int nx, int ny, int nz,
                                  const char* x1_tag, 
                                  const char* x2_tag, 
                                  const char* y1_tag,
                                  const char* y2_tag,
                                  const char* z1_tag,
                                  const char* z2_tag)
{
  // Tag the boundaries of the mesh.
  cubic_lattice_t* lattice = mesh_property(mesh, "lattice");
  ASSERT(lattice != NULL);
  int* x1tag = mesh_create_tag(mesh->face_tags, x1_tag, ny*nz);
  int* x2tag = mesh_create_tag(mesh->face_tags, x2_tag, ny*nz);
  for (int j = 0; j < ny; ++j)
  {
    for (int k = 0; k < nz; ++k)
    {
      x1tag[nz*j + k] = cubic_lattice_x_face(lattice, 0, j, k);
      x2tag[nz*j + k] = cubic_lattice_x_face(lattice, nx, j, k);
    }
  }

  int* y1tag = mesh_create_tag(mesh->face_tags, y1_tag, nx*nz);
  int* y2tag = mesh_create_tag(mesh->face_tags, y2_tag, nx*nz);
  for (int i = 0; i < nx; ++i)
  {
    for (int k = 0; k < nz; ++k)
    {
      y1tag[nz*i + k] = cubic_lattice_y_face(lattice, i, 0, k);
      y2tag[nz*i + k] = cubic_lattice_y_face(lattice, i, ny, k);
    }
  }

  int* z1tag = mesh_create_tag(mesh->face_tags, z1_tag, nx*ny);
  int* z2tag = mesh_create_tag(mesh->face_tags, z2_tag, nx*ny);
  for (int i = 0; i < nx; ++i)
  {
    for (int j = 0; j < ny; ++j)
    {
      z1tag[ny*i + j] = cubic_lattice_z_face(lattice, i, j, 0);
      z2tag[ny*i + j] = cubic_lattice_z_face(lattice, i, j, nz);
    }
  }
  lattice = NULL;
}

