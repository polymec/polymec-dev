#include "core/edit_mesh.h"
#include "core/unordered_set.h"
#include "geometry/create_cubic_lattice_mesh.h"
#include "geometry/cubic_lattice.h"

#ifdef __cplusplus
extern "C" {
#endif

mesh_t* create_cubic_lattice_mesh(int nx, int ny, int nz)
{
  // Create a cubic lattice object for indexing.
  cubic_lattice_t* lattice = cubic_lattice_new(nx, ny, nz);

  // Create the mesh.
  // FIXME: Not parallel safe.
  mesh_t* mesh = mesh_new(cubic_lattice_num_cells(lattice), 0,
                          cubic_lattice_num_faces(lattice),
                          cubic_lattice_num_edges(lattice),
                          cubic_lattice_num_nodes(lattice));

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
          mesh_add_face_to_cell(mesh, &mesh->faces[faces[f]], &mesh->cells[cell]);

          if (!int_unordered_set_contains(processed_faces, faces[f]))
          {
            for (int e = 0; e < 4; ++e)
            {
              mesh_add_edge_to_face(mesh, &mesh->edges[edges[f][e]], &mesh->faces[faces[f]]);
              mesh->edges[edges[f][e]].node1 = &mesh->nodes[nodes[f][e][0]];
              mesh->edges[edges[f][e]].node2 = &mesh->nodes[nodes[f][e][1]];
            }
            int_unordered_set_insert(processed_faces, faces[f]);
          }
        }

        // Assign the node positions for a uniform grid spanning [0,1]x[0,1]x[0,1].
        double dx = 1.0/nx, dy = 1.0/ny, dz = 1.0/nz;
        static const int i_offsets[] = {0, 1, 1, 0, 0, 1, 1, 0};
        static const int j_offsets[] = {0, 0, 1, 1, 0, 0, 1, 1};
        static const int k_offsets[] = {0, 0, 0, 0, 1, 1, 1, 1};
        for (int n = 0; n < 8; ++n)
        {
          if (!int_unordered_set_contains(processed_nodes, node_indices[n]))
          {
            node_t* node = &mesh->nodes[node_indices[n]];
            node->x = (i + i_offsets[n]) * dx;
            node->y = (j + j_offsets[n]) * dy;
            node->z = (k + k_offsets[n]) * dz;
            int_unordered_set_insert(processed_nodes, node_indices[n]);
          }
        }
      }
    }
  }

  // Clean up.
  int_unordered_set_free(processed_faces);
  int_unordered_set_free(processed_nodes);

  return mesh;
}

#ifdef __cplusplus
}
#endif

