#include "core/edit_mesh.h"
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

        // Hook up each face to its edges. 
        // edges[i][j] is the jth edge for the ith face.
        int edges[6][4];

        edges[0][0] = cubic_lattice_y_edge(lattice, i, j, k);
        edges[0][1] = cubic_lattice_z_edge(lattice, i, j, k+1);
        edges[0][2] = cubic_lattice_y_edge(lattice, i, j+1, k);
        edges[0][2] = cubic_lattice_z_edge(lattice, i, j, k);

        edges[1][0] = cubic_lattice_y_edge(lattice, i+1, j, k);
        edges[1][1] = cubic_lattice_z_edge(lattice, i+1, j, k+1);
        edges[1][2] = cubic_lattice_y_edge(lattice, i+1, j+1, k);
        edges[1][2] = cubic_lattice_z_edge(lattice, i+1, j, k);

        // FIXME
        edges[2][0] = cubic_lattice_x_edge(lattice, i, j, k);
        edges[2][1] = cubic_lattice_z_edge(lattice, i, j, k+1);
        edges[2][2] = cubic_lattice_x_edge(lattice, i, j+1, k+1);
        edges[2][2] = cubic_lattice_z_edge(lattice, i, j, k);

        edges[3][0] = cubic_lattice_x_edge(lattice, i+1, j, k);
        edges[3][1] = cubic_lattice_z_edge(lattice, i+1, j, k+1);
        edges[3][2] = cubic_lattice_x_edge(lattice, i+1, j+1, k);
        edges[3][2] = cubic_lattice_z_edge(lattice, i+1, j, k);

        edges[4][0] = cubic_lattice_x_edge(lattice, i, j, k);
        edges[4][1] = cubic_lattice_y_edge(lattice, i, j, k+1);
        edges[4][2] = cubic_lattice_x_edge(lattice, i, j+1, k+1);
        edges[4][2] = cubic_lattice_y_edge(lattice, i, j, k);

        edges[5][0] = cubic_lattice_x_edge(lattice, i+1, j, k);
        edges[5][1] = cubic_lattice_y_edge(lattice, i+1, j, k+1);
        edges[5][2] = cubic_lattice_x_edge(lattice, i+1, j+1, k);
        edges[5][2] = cubic_lattice_y_edge(lattice, i+1, j, k);

        // Hook up the nodes to their edges.
        // edges[i][j][k] is the kth node for the jth edge for the ith face.
        int nodes[6][4][2];

        // FIXME

        // Hook everything up.
        for (int f = 0; f < 6; ++f)
        {
          mesh_add_face_to_cell(mesh, &mesh->faces[faces[f]], &mesh->cells[cell]);
          for (int e = 0; e < 4; ++e)
          {
            mesh_add_edge_to_face(mesh, &mesh->edges[edges[f][e]], &mesh->faces[faces[f]]);
            mesh->edges[edges[f][e]].node1 = &mesh->nodes[nodes[f][e][0]];
            mesh->edges[edges[f][e]].node2 = &mesh->nodes[nodes[f][e][1]];
          }
        }
      }
    }
  }

  return mesh;
}

#ifdef __cplusplus
}
#endif

