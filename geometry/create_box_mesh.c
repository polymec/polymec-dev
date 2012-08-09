#include "geometry/create_box_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

mesh_t* create_box_mesh(int N[], double low[], double high[])
{
  // FIXME: Only serial for now!
  return NULL;

#if 0
  // Allocate the basic elements.
  int num_cells = N[0]*N[1]*N[2];
  int num_ghost_cells = 0; // FIXME
  int num_faces = (N[0]+1)*N[1]*N[2] + 
                  N[0]*(N[1]+1)*N[2] + 
                  N[0]*N[1]*(N[2]+1);
  int num_edges = N[0]*(N[1]+1)*(N[2]+1) + 
                  (N[0]+1)*N[1]*(N[2]+1) + 
                  (N[0]+1)*(N[1]+1)*N[2];
  int num_nodes = (N[0]+1)*(N[1]+1)*(N[2]+1);
  mesh_t* mesh = mesh_new(num_cells, num_ghost_cells, num_faces, 
                          num_edges, num_nodes);

  // All cells are hexahedra.
  for (int i = 0; i < num_cells; ++i)
  {
    mesh->cells[i]->faces = calloc(6, sizeof(face_t*));
    mesh->cells[i]->num_faces = 6;
    mesh->cells[i]->volume = V;
  }

  // All faces are rectangles.
  for (int i = 0; i < num_faces; ++i)
  {
    mesh->faces[i]->edges = calloc(4, sizeof(edge_t*));
    mesh->faces[i]->num_edges = 4;
  }

  // Now let's hook everything together.

  // cells <-> faces
  int face_index = 0;
  for (int i = 0; i < num_cells; ++i)
  {
    // Determine lattice coordinates.
    int I = i % N[0];
    int J = i/N[0] % N[1];
    int K = i/(N[0]*N[1]) % N[2];
    cell_t* cell = mesh->cells[i];

    // Hook up x-faces.
    if (I > 0) // non-left-most cell
    {
      cell->faces[0] = mesh->cells[i-1]->faces[1];
      cell->faces[0]->cell2 = cell;
    }
    else // left-most cell
    {
      cell->faces[0] = mesh->faces[face_index++]; // ghost face
    }
    cell->faces[1] = mesh->faces[face_index++]; // right x-face.
    cell->faces[1]->cell1 = cell;

    // Hook up y-faces.
    if (J > 0) // non-bottom-most cell
    {
      cell->faces[2] = mesh->cells[i-N[0]]->faces[3];
      cell->faces[2]->cell2 = cell;
    }
    else // bottom-most cell
    {
      cell->faces[2] = mesh->faces[face_index++]; // ghost face
    }
    cell->faces[3] = mesh->faces[face_index++]; // top y-face.
    cell->faces[3]->cell1 = cell;

    // Hook up z-faces.
    if (K > 0) // non-inner-most cell
    {
      cell->faces[4] = mesh->cells[i-N[0]*N[1]]->faces[5];
      cell->faces[4]->cell2 = cell;
    }
    else // inner-most cell
    {
      cell->faces[4] = mesh->faces[face_index++]; // ghost face
    }
    cell->faces[5] = mesh->faces[face_index++]; // outer z-face.
    cell->faces[5]->cell1 = cell;
  }

  // Hook up ghost cells.
  for (int i = num_cells; i < num_cells + num_ghost_cells; ++i)
  {

  }

  // faces <-> edges
  int edge_index = 0;
  for (int i = 0; i < num_faces; ++i)
  {
  }

  // edges <-> nodes
  int node_index = 0;
  for (int i = 0; i < num_edges; ++i)
  {
  }

  // Place the nodes in space.
  double Lx = high[0] - low[0],
         Ly = high[1] - low[1],
         Lz = high[2] - low[2],
         dx = Lx/N[0], dy = Ly/N[1], dz = Lz/N[2],
         V = dx*dy*dz, Ax = dy*dz, Ay = dx*dz, Az = dx*dy;
  for (int c = 0; c < num_cells; ++c)
  {
    // Determine lattice coordinates.
    int I = i % N[0];
    int J = i/N[0] % N[1];
    int K = i/(N[0]*N[1]) % N[2];
    node_t pos[3] = {.x = I*dx, .y = J*dy, .z = K*dz};
    cell_t* cell = mesh->cells[c];

    // -x face
    face_t* face = cell->faces[0]; // -x face
    face->area = Ax;
    face->edges[0]->node1 = pos;
    pos[1]
    face->edges[0]->node1 = pos;


    }
    cell->volume = V;
  }

  return mesh;
#endif
}

#ifdef __cplusplus
}
#endif

