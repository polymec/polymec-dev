#include "core/edit_mesh.h"
#include "core/unordered_set.h"
#include "geometry/create_cubic_lattice_mesh.h"
#include "geometry/cubic_lattice.h"

#ifdef __cplusplus
extern "C" {
#endif

mesh_t* create_cubic_lattice_mesh_with_bbox(int nx, int ny, int nz, bbox_t* bbox)
{
  ASSERT(nx > 0);
  ASSERT(ny > 0);
  ASSERT(nz > 0);

  ASSERT((bbox == NULL) || (bbox->x2 > bbox->x1));
  ASSERT((bbox == NULL) || (bbox->y2 > bbox->y1));
  ASSERT((bbox == NULL) || (bbox->z2 > bbox->z1));
//  ASSERT(num_ghost >= 0);

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

  // Grid spacings.
  double Lx = (bbox != NULL) ? bbox->x2 - bbox->x1 : 1.0;
  double Ly = (bbox != NULL) ? bbox->y2 - bbox->y1 : 1.0;
  double Lz = (bbox != NULL) ? bbox->z2 - bbox->z1 : 1.0;
  double x1 = (bbox != NULL) ? bbox->x1 : 0.0;
  double y1 = (bbox != NULL) ? bbox->y1 : 0.0;
  double z1 = (bbox != NULL) ? bbox->z1 : 0.0;
  double dx = Lx/nx, dy = Ly/ny, dz = Lz/nz;

  int_unordered_set_t* processed_faces = int_unordered_set_new();
  int_unordered_set_t* processed_nodes = int_unordered_set_new();
  for (int k = 0; k < nz; ++k)
  {
    for (int j = 0; j < ny; ++j)
    {
      for (int i = 0; i < nx; ++i)
      {
        // Hook up the cell and faces.
        int cell = cubic_lattice_cell(lattice, i, j, k);
        int faces[6];
        faces[0] = cubic_lattice_x_face(lattice, i, j, k);
        faces[1] = cubic_lattice_x_face(lattice, i+1, j, k);
        faces[2] = cubic_lattice_y_face(lattice, i, j, k);
        faces[3] = cubic_lattice_y_face(lattice, i, j+1, k);
        faces[4] = cubic_lattice_z_face(lattice, i, j, k);
        faces[5] = cubic_lattice_z_face(lattice, i, j, k+1);

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

        // Face 0 (-x)
        edges[0][0] = edge_indices[3];
        nodes[0][0][0] = node_indices[3];
        nodes[0][0][1] = node_indices[0];

        edges[0][1] = edge_indices[4];
        nodes[0][1][0] = node_indices[0];
        nodes[0][1][1] = node_indices[4];

        edges[0][2] = edge_indices[11];
        nodes[0][2][0] = node_indices[4];
        nodes[0][2][1] = node_indices[7];

        edges[0][3] = edge_indices[7];
        nodes[0][3][0] = node_indices[7];
        nodes[0][3][1] = node_indices[3];

        // Face 1 (+x)
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

        // Face 2 (-y)
        edges[2][0] = edge_indices[0];
        nodes[2][0][0] = node_indices[0];
        nodes[2][0][1] = node_indices[1];

        edges[2][1] = edge_indices[5];
        nodes[2][1][0] = node_indices[1];
        nodes[2][1][1] = node_indices[5];

        edges[2][2] = edge_indices[8];
        nodes[2][2][0] = node_indices[5];
        nodes[2][2][1] = node_indices[4];

        edges[2][3] = edge_indices[4];
        nodes[2][3][0] = node_indices[4];
        nodes[2][3][1] = node_indices[0];

        // Face 3 (+y)
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

        // Face 4 (-z)
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

        // Face 5 (+z)
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
          face_t* face = &mesh->faces[faces[f]];
          mesh_add_face_to_cell(mesh, face, &mesh->cells[cell]);

          if (!int_unordered_set_contains(processed_faces, faces[f]))
          {
            for (int e = 0; e < 4; ++e)
            {
              mesh_add_edge_to_face(mesh, &mesh->edges[edges[f][e]], face);
              mesh->edges[edges[f][e]].node1 = &mesh->nodes[nodes[f][e][0]];
              mesh->edges[edges[f][e]].node2 = &mesh->nodes[nodes[f][e][1]];
            }

            // Face geometry.
            switch (f)
            {
              case 0: // -x
                face->area = Ly*Lz / (ny*nz);
                face->center.x = x1 + i*dx;
                face->center.y = y1 + (j+0.5)*dy;
                face->center.z = z1 + (k+0.5)*dz;
                break;
              case 1: // +x
                face->area = Ly*Lz / (ny*nz);
                face->center.x = x1 + (i+1)*dx;
                face->center.y = y1 + (j+0.5)*dy;
                face->center.z = z1 + (k+0.5)*dz;
                break;
              case 2: // -y
                face->area = Lx*Lz / (nx*nz);
                face->center.x = x1 + (i+0.5)*dx;
                face->center.y = y1 + j*dy;
                face->center.z = z1 + (k+0.5)*dz;
                break;
              case 3: // +y
                face->area = Lx*Lz / (nx*nz);
                face->center.x = x1 + (i+0.5)*dx;
                face->center.y = y1 + (j+1)*dy;
                face->center.z = z1 + (k+0.5)*dz;
                break;
              case 4: // -z
                face->area = Lx*Ly / (nx*ny);
                face->center.x = x1 + (i+0.5)*dx;
                face->center.y = y1 + (j+0.5)*dy;
                face->center.z = z1 + k*dz;
                break;
              case 5: // +z
                face->area = Lx*Ly / (nx*ny);
                face->center.x = x1 + (i+0.5)*dx;
                face->center.y = y1 + (j+0.5)*dy;
                face->center.z = z1 + (k+1)*dz;
                break;
            }

            // We're done processing this face.
            int_unordered_set_insert(processed_faces, faces[f]);
          }
        }

        // Assign the node positions for a uniform grid spanning [0,1]x[0,1]x[0,1].
        static const int i_offsets[] = {0, 1, 1, 0, 0, 1, 1, 0};
        static const int j_offsets[] = {0, 0, 1, 1, 0, 0, 1, 1};
        static const int k_offsets[] = {0, 0, 0, 0, 1, 1, 1, 1};
        for (int n = 0; n < 8; ++n)
        {
          if (!int_unordered_set_contains(processed_nodes, node_indices[n]))
          {
            node_t* node = &mesh->nodes[node_indices[n]];
            node->x = x1 + (i + i_offsets[n]) * dx;
            node->y = y1 + (j + j_offsets[n]) * dy;
            node->z = z1 + (k + k_offsets[n]) * dz;
            int_unordered_set_insert(processed_nodes, node_indices[n]);
          }
        }

        // Cell geometry.
        mesh->cells[cell].volume = Lx*Ly*Lz / (nx*ny*nz);
        mesh->cells[cell].center.x = x1 + (i+0.5)*dx;
        mesh->cells[cell].center.y = y1 + (j+0.5)*dy;
        mesh->cells[cell].center.z = z1 + (k+0.5)*dz;
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
        mesh_add_face_to_cell(mesh, low_face, low_ghost_cell);

        int high_cell = cubic_lattice_cell(lattice, nx-1, j, k);
        cell_t* high_ghost_cell = &mesh->cells[gindex++];
        face_t* high_face = mesh->cells[high_cell].faces[1];
        ASSERT(high_face->cell2 == NULL); 
        high_face->cell2 = high_ghost_cell;
        mesh_add_face_to_cell(mesh, high_face, high_ghost_cell);

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
        mesh_add_face_to_cell(mesh, low_face, low_ghost_cell);

        int high_cell = cubic_lattice_cell(lattice, i, ny-1, k);
        cell_t* high_ghost_cell = &mesh->cells[gindex++];
        face_t* high_face = mesh->cells[high_cell].faces[3];
        ASSERT(high_face->cell2 == NULL); 
        high_face->cell2 = high_ghost_cell;
        mesh_add_face_to_cell(mesh, high_face, high_ghost_cell);

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
        mesh_add_face_to_cell(mesh, low_face, low_ghost_cell);

        int high_cell = cubic_lattice_cell(lattice, i, j, nz-1);
        cell_t* high_ghost_cell = &mesh->cells[gindex++];
        face_t* high_face = mesh->cells[high_cell].faces[5];
        ASSERT(high_face->cell2 == NULL); 
        high_face->cell2 = high_ghost_cell;
        mesh_add_face_to_cell(mesh, high_face, high_ghost_cell);

        // FIXME: Do ghost cells need node/edge connectivity?
      }
    }
  }
  ASSERT(gindex == (mesh->num_cells + mesh->num_ghost_cells));
#endif

  // Clean up.
  int_unordered_set_free(processed_faces);
  int_unordered_set_free(processed_nodes);

  return mesh;
}

mesh_t* create_cubic_lattice_mesh(int nx, int ny, int nz)
{
  bbox_t bbox = {.x1 = 0.0, .x2 = 1.0, .y1 = 0.0, .y2 = 1.0, .z1 = 0.0, .z2 = 1.0};
  return create_cubic_lattice_mesh_with_bbox(nx, ny, nz, &bbox);
}

#ifdef __cplusplus
}
#endif

