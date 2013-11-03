// Copyright 2012-2013 Jeffrey Johnson.
// 
// This file is part of Polymec, and is licensed under the Apache License, 
// Version 2.0 (the "License"); you may not use this file except in 
// compliance with the License. You may may find the text of the license in 
// the LICENSE file at the top-level source directory, or obtain a copy of 
// it at
// 
// http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "core/unordered_set.h"
#include "geometry/create_rectilinear_lattice_mesh.h"
#include "geometry/cubic_lattice.h"

mesh_t* create_rectilinear_lattice_mesh(double* xs, int nxs, 
                                        double* ys, int nys, 
                                        double* zs, int nzs)
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
  mesh_t* mesh = mesh_new(cubic_lattice_num_cells(lattice), 0,
                          cubic_lattice_num_faces(lattice),
                          cubic_lattice_num_edges(lattice),
                          cubic_lattice_num_nodes(lattice));
  mesh->cell_faces = ARENA_REALLOC(mesh->arena, mesh->cell_faces, sizeof(int)*6*mesh->num_cells, 0);
  mesh->face_nodes = ARENA_REALLOC(mesh->arena, mesh->face_nodes, sizeof(int)*4*mesh->num_faces, 0);
  mesh->face_edges = ARENA_REALLOC(mesh->arena, mesh->face_nodes, sizeof(int)*4*mesh->num_faces, 0);
  mesh->cell_face_offsets[mesh->num_cells] = 6*mesh->num_cells;
  mesh->face_node_offsets[mesh->num_faces] = 4*mesh->num_faces;
  mesh->face_edge_offsets[mesh->num_faces] = 4*mesh->num_faces;
  for (int f = 0; f < mesh->num_faces; ++f)
  {
    mesh->face_cells[2*f] = -1;
    mesh->face_cells[2*f+1] = -1;
  }

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

        // Hook up each face to its edges, and each edge to its nodes. 
        // edges[i][j] is the jth edge for the ith face.
        // nodes[i][j][k] is the kth node for the jth edge for the ith face.
        int nodes[6][4][2];
        int edges[6][4];

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

        // The edges 0-3 traverse the bottom from 0->1->2->3->0.
        int edge_indices[12];
        edge_indices[0] = cubic_lattice_x_edge(lattice, i, j, k);
        edge_indices[1] = cubic_lattice_y_edge(lattice, i+1, j, k);
        edge_indices[2] = cubic_lattice_x_edge(lattice, i, j+1, k);
        edge_indices[3] = cubic_lattice_y_edge(lattice, i, j, k);

        // The edges 4-7 scale the sides, connecting 0-3 to 4-7.
        edge_indices[4] = cubic_lattice_z_edge(lattice, i, j, k);
        edge_indices[5] = cubic_lattice_z_edge(lattice, i+1, j, k);
        edge_indices[6] = cubic_lattice_z_edge(lattice, i+1, j+1, k);
        edge_indices[7] = cubic_lattice_z_edge(lattice, i, j+1, k);

        // The edges 8-11 traverse the top from 4->5->6->7->4.
        edge_indices[8]  = cubic_lattice_x_edge(lattice, i, j, k+1);
        edge_indices[9]  = cubic_lattice_y_edge(lattice, i+1, j, k+1);
        edge_indices[10] = cubic_lattice_x_edge(lattice, i, j+1, k+1);
        edge_indices[11] = cubic_lattice_y_edge(lattice, i, j, k+1);

        // Face 0 (-x) -- backward traversal
        edges[0][0] = edge_indices[7];
        nodes[0][0][0] = node_indices[7];
        nodes[0][0][1] = node_indices[4];

        edges[0][1] = edge_indices[11];
        nodes[0][1][0] = node_indices[4];
        nodes[0][1][1] = node_indices[0];

        edges[0][2] = edge_indices[4];
        nodes[0][2][0] = node_indices[0];
        nodes[0][2][1] = node_indices[3];

        edges[0][3] = edge_indices[3];
        nodes[0][3][0] = node_indices[3];
        nodes[0][3][1] = node_indices[7];

        // Face 1 (+x) -- forward traversal
        edges[1][0] = edge_indices[1];
        nodes[1][0][0] = node_indices[1];
        nodes[1][0][1] = node_indices[2];

        edges[1][1] = edge_indices[6];
        nodes[1][1][0] = node_indices[2];
        nodes[1][1][1] = node_indices[6];

        edges[1][2] = edge_indices[9];
        nodes[1][2][0] = node_indices[6];
        nodes[1][2][1] = node_indices[5];

        edges[1][3] = edge_indices[5];
        nodes[1][3][0] = node_indices[5];
        nodes[1][3][1] = node_indices[1];

        // Face 2 (-y) -- backward traversal
        edges[2][0] = edge_indices[4];
        nodes[2][0][0] = node_indices[4];
        nodes[2][0][1] = node_indices[5];

        edges[2][1] = edge_indices[8];
        nodes[2][1][0] = node_indices[5];
        nodes[2][1][1] = node_indices[1];

        edges[2][2] = edge_indices[5];
        nodes[2][2][0] = node_indices[1];
        nodes[2][2][1] = node_indices[0];

        edges[2][3] = edge_indices[0];
        nodes[2][3][0] = node_indices[0];
        nodes[2][3][1] = node_indices[4];

        // Face 3 (+y) -- forward traversal
        edges[3][0] = edge_indices[2];
        nodes[3][0][0] = node_indices[2];
        nodes[3][0][1] = node_indices[3];

        edges[3][1] = edge_indices[7];
        nodes[3][1][0] = node_indices[3];
        nodes[3][1][1] = node_indices[7];

        edges[3][2] = edge_indices[10];
        nodes[3][2][0] = node_indices[7];
        nodes[3][2][1] = node_indices[6];

        edges[3][3] = edge_indices[6];
        nodes[3][3][0] = node_indices[6];
        nodes[3][3][1] = node_indices[2];

        // Face 4 (-z) -- backward traversal
        edges[4][0] = edge_indices[0];
        nodes[4][0][0] = node_indices[0];
        nodes[4][0][1] = node_indices[1];

        edges[4][1] = edge_indices[1];
        nodes[4][1][0] = node_indices[1];
        nodes[4][1][1] = node_indices[2];

        edges[4][2] = edge_indices[2];
        nodes[4][2][0] = node_indices[2];
        nodes[4][2][1] = node_indices[3];

        edges[4][3] = edge_indices[3];
        nodes[4][3][0] = node_indices[3];
        nodes[4][3][1] = node_indices[0];

        // Face 5 (+z) -- forward traversal
        edges[5][0] = edge_indices[8];
        nodes[5][0][0] = node_indices[4];
        nodes[5][0][1] = node_indices[5];

        edges[5][1] = edge_indices[9];
        nodes[5][1][0] = node_indices[5];
        nodes[5][1][1] = node_indices[6];

        edges[5][2] = edge_indices[10];
        nodes[5][2][0] = node_indices[6];
        nodes[5][2][1] = node_indices[7];

        edges[5][3] = edge_indices[11];
        nodes[5][3][0] = node_indices[7];
        nodes[5][3][1] = node_indices[4];

        // Hook everything up.
        for (int f = 0; f < 6; ++f)
        {
          int face = mesh->cell_faces[6*cell+f];
          if (face < 0) face = ~face;
          if (mesh->face_cells[2*face] == -1)
            mesh->face_cells[2*face] = cell;
          else if (mesh->face_cells[2*face] == -1)
            mesh->face_cells[2*face+1] = cell;

          mesh->face_edge_offsets[face] = 4*face;
          for (int e = 0; e < 4; ++e)
          {
            int edge = edges[f][e];
            mesh->face_edges[4*face+e] = edge;
            mesh->edge_nodes[2*edge] = nodes[f][e][0];
            mesh->edge_nodes[2*edge+1] = nodes[f][e][1];
          }

          mesh->face_node_offsets[face] = 4*face;
          for (int n = 0; n < 4; ++n)
            mesh->face_nodes[4*face+n] = nodes[f][n][0];
        }

        // Assign the node positions for a uniform grid spanning [0,1]x[0,1]x[0,1].
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

  // Compute mesh geometry.
  mesh_compute_geometry(mesh);

  // Stash the lattice in the "lattice" property.
  mesh_set_property(mesh, "lattice", (void*)lattice, NULL);

  // Clean up.
  int_unordered_set_free(processed_nodes);

  return mesh;
}

