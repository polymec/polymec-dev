#include "geometry/create_block_mesh.h"
#include "core/index_space.h"
#include "core/unordered_set.h"
#include "core/edit_mesh.h"
#include "geometry/cubic_lattice.h"

#ifdef __cplusplus
extern "C" {
#endif

// This data structure represents a hexahedral cell, giving relationships 
// between indices of mesh elements.
typedef struct
{
  int faces[6], edges[12], nodes[8];
} hexahedral_cell_t;

static void create_block_grid(block_t* block, 
                              int low_cell, 
                              int high_cell,
                              hexahedral_cell_t* hexes, 
                              mesh_t* mesh)
{
  ASSERT(block != NULL);
  ASSERT(num_block_cells > 0);
  ASSERT(hexes != NULL);
  ASSERT(mesh != NULL);

  int_unordered_set_t* processed_faces = int_unordered_set_new();
  int_unordered_set_t* processed_nodes = int_unordered_set_new();
  for (int c = low_cell; c < high_cell; ++c)
  {
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

    // Face 0 (-x)
    edges[0][0] = hexes[c].edges[3];
    nodes[0][0][0] = hexes[c].nodes[3];
    nodes[0][0][1] = hexes[c].nodes[0];

    edges[0][1] = hexes[c].edges[4];
    nodes[0][1][0] = hexes[c].nodes[0];
    nodes[0][1][1] = hexes[c].nodes[4];

    edges[0][2] = hexes[c].edges[11];
    nodes[0][2][0] = hexes[c].nodes[4];
    nodes[0][2][1] = hexes[c].nodes[7];

    edges[0][3] = hexes[c].edges[7];
    nodes[0][3][0] = hexes[c].nodes[7];
    nodes[0][3][1] = hexes[c].nodes[3];

    // Face 1 (+x)
    edges[1][0] = hexes[c].edges[1];
    nodes[1][0][0] = hexes[c].nodes[1];
    nodes[1][0][1] = hexes[c].nodes[2];

    edges[1][1] = hexes[c].edges[6];
    nodes[1][1][0] = hexes[c].nodes[2];
    nodes[1][1][1] = hexes[c].nodes[6];

    edges[1][2] = hexes[c].edges[9];
    nodes[1][2][0] = hexes[c].nodes[6];
    nodes[1][2][1] = hexes[c].nodes[5];

    edges[1][3] = hexes[c].edges[5];
    nodes[1][3][0] = hexes[c].nodes[5];
    nodes[1][3][1] = hexes[c].nodes[1];

    // Face 2 (-y)
    edges[2][0] = hexes[c].edges[0];
    nodes[2][0][0] = hexes[c].nodes[0];
    nodes[2][0][1] = hexes[c].nodes[1];

    edges[2][1] = hexes[c].edges[5];
    nodes[2][1][0] = hexes[c].nodes[1];
    nodes[2][1][1] = hexes[c].nodes[5];

    edges[2][2] = hexes[c].edges[8];
    nodes[2][2][0] = hexes[c].nodes[5];
    nodes[2][2][1] = hexes[c].nodes[4];

    edges[2][3] = hexes[c].edges[4];
    nodes[2][3][0] = hexes[c].nodes[4];
    nodes[2][3][1] = hexes[c].nodes[0];

    // Face 3 (+y)
    edges[3][0] = hexes[c].edges[2];
    nodes[3][0][0] = hexes[c].nodes[2];
    nodes[3][0][1] = hexes[c].nodes[3];

    edges[3][1] = hexes[c].edges[7];
    nodes[3][1][0] = hexes[c].nodes[3];
    nodes[3][1][1] = hexes[c].nodes[7];

    edges[3][2] = hexes[c].edges[10];
    nodes[3][2][0] = hexes[c].nodes[7];
    nodes[3][2][1] = hexes[c].nodes[6];

    edges[3][3] = hexes[c].edges[6];
    nodes[3][3][0] = hexes[c].nodes[6];
    nodes[3][3][1] = hexes[c].nodes[2];

    // Face 4 (-z)
    edges[4][0] = hexes[c].edges[0];
    nodes[4][0][0] = hexes[c].nodes[0];
    nodes[4][0][1] = hexes[c].nodes[1];

    edges[4][1] = hexes[c].edges[1];
    nodes[4][1][0] = hexes[c].nodes[1];
    nodes[4][1][1] = hexes[c].nodes[2];

    edges[4][2] = hexes[c].edges[2];
    nodes[4][2][0] = hexes[c].nodes[2];
    nodes[4][2][1] = hexes[c].nodes[3];

    edges[4][3] = hexes[c].edges[3];
    nodes[4][3][0] = hexes[c].nodes[3];
    nodes[4][3][1] = hexes[c].nodes[0];

    // Face 5 (+z)
    edges[5][0] = hexes[c].edges[8];
    nodes[5][0][0] = hexes[c].nodes[4];
    nodes[5][0][1] = hexes[c].nodes[5];

    edges[5][1] = hexes[c].edges[9];
    nodes[5][1][0] = hexes[c].nodes[5];
    nodes[5][1][1] = hexes[c].nodes[6];

    edges[5][2] = hexes[c].edges[10];
    nodes[5][2][0] = hexes[c].nodes[6];
    nodes[5][2][1] = hexes[c].nodes[7];

    edges[5][3] = hexes[c].edges[11];
    nodes[5][3][0] = hexes[c].nodes[7];
    nodes[5][3][1] = hexes[c].nodes[4];

    // Hook everything up.
    for (int f = 0; f < 6; ++f)
    {
      face_t* face = &mesh->faces[hexes[c].faces[f]];
      mesh_add_face_to_cell(mesh, face, &mesh->cells[c]);

      if (!int_unordered_set_contains(processed_faces, hexes[c].faces[f]))
      {
        for (int e = 0; e < 4; ++e)
        {
          mesh_add_edge_to_face(mesh, &mesh->edges[edges[f][e]], face);
          mesh->edges[edges[f][e]].node1 = &mesh->nodes[nodes[f][e][0]];
          mesh->edges[edges[f][e]].node2 = &mesh->nodes[nodes[f][e][1]];
        }

#if 0
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
#endif

        // We're done processing this face.
        int_unordered_set_insert(processed_faces, hexes[c].faces[f]);
      }
    }

    // Assign the node positions for a uniform grid spanning [0,1]x[0,1]x[0,1].
    static const int i_offsets[] = {0, 1, 1, 0, 0, 1, 1, 0};
    static const int j_offsets[] = {0, 0, 1, 1, 0, 0, 1, 1};
    static const int k_offsets[] = {0, 0, 0, 0, 1, 1, 1, 1};
    for (int n = 0; n < 8; ++n)
    {
      if (!int_unordered_set_contains(processed_nodes, hexes[c].nodes[n]))
      {
#if 0
        node_t* node = &mesh->nodes[hexes.nodes[n]];
        node->x = x1 + (i + i_offsets[n]) * dx;
        node->y = y1 + (j + j_offsets[n]) * dy;
        node->z = z1 + (k + k_offsets[n]) * dz;
#endif
        int_unordered_set_insert(processed_nodes, hexes[c].nodes[n]);
      }
    }

#if 0
    // Cell geometry.
    mesh->cells[cell].volume = Lx*Ly*Lz / (nx*ny*nz);
    mesh->cells[cell].center.x = x1 + (i+0.5)*dx;
    mesh->cells[cell].center.y = y1 + (j+0.5)*dy;
    mesh->cells[cell].center.z = z1 + (k+0.5)*dz;
#endif
  }
}

mesh_t* create_block_mesh(MPI_Comm comm, 
                          block_assembly_t* assembly, 
                          block_grid_t* resolutions, 
                          int num_ghosts)
{
  // Compute the number of cells in the assembly.
  int num_blocks = block_assembly_num_blocks(assembly), 
      num_global_cells = 0, num_block_cells[num_blocks];
  cubic_lattice_t* lattices[num_blocks];
  for (int b = 0; b < num_blocks; ++b)
  {
    int N1 = resolutions[b].N1, N2 = resolutions[b].N2, N3 = resolutions[b].N3;
    lattices[b] = cubic_lattice_new(N1, N2, N3);
    block_t* block = block_assembly_block(assembly, b);
    if (block != NULL)
    {
      num_block_cells[b] = cubic_lattice_num_cells(lattices[b]);
      num_global_cells += num_block_cells[b];
    }
    else
      num_block_cells[b] = 0;
  }

  // Also, compute cell, face, edge, node offsets for each block.
  int cell_block_offsets[num_blocks], face_block_offsets[num_blocks],
      edge_block_offsets[num_blocks], node_block_offsets[num_blocks];
  for (int b = 0; b < num_blocks; ++b)
  {
    cell_block_offsets[b] = face_block_offsets[b] = 0;
    edge_block_offsets[b] = node_block_offsets[b] = 0;
    if (b == 0)
    {
      cubic_lattice_t* lattice = lattices[b];
      cell_block_offsets[b] = cubic_lattice_num_cells(lattice);
      face_block_offsets[b] = cubic_lattice_num_faces(lattice);
      edge_block_offsets[b] = cubic_lattice_num_edges(lattice);
      node_block_offsets[b] = cubic_lattice_num_nodes(lattice);
    }
    else
    {
      for (int bb = 0; bb < b; ++bb)
      {
        cell_block_offsets[b] += cell_block_offsets[bb];
        face_block_offsets[b] += face_block_offsets[bb];
        edge_block_offsets[b] += edge_block_offsets[bb];
        node_block_offsets[b] += node_block_offsets[bb];
      }
    }
  }

  // Create an index space for this set of cells, dividing the work 
  // evenly across all processes. This process will operate only on cells 
  // [is->low, is->high).
  index_space_t* is = index_space_from_naive_partitions(comm, num_global_cells);
  int num_cells = is->high - is->low;

  // Which blocks will the current process operate on? And which cells within
  // those blocks?
  int low_block = 0, low_block_cell = 0, 
      high_block = 0, high_block_cell = 0, offset = 0;
  for (int i = 0; i < num_blocks; ++i)
  {
    if (offset < is->low)
      low_block++;
    else
      low_block_cell = offset - is->low;
    offset += num_block_cells[i];
    if (offset <= is->high)
      high_block++;
    else 
      high_block_cell = is->high - offset;
  }

  // Initialize a record that allows us to connect mesh elements by index.
  // We'll use this to hook together everything once we've built the mesh.
  // Note that this array only contains information for local cells, and not 
  // ghost cells.
  hexahedral_cell_t* hexes = malloc(sizeof(hexahedral_cell_t)*num_cells);

  // Keep track of elements we've already counted.
  int_unordered_set_t *node_counted = int_unordered_set_new(),
                      *edge_counted = int_unordered_set_new(),
                      *face_counted = int_unordered_set_new();

  // Figure out the numbers of local faces, edges, nodes, and ghost cells. 
  // These are the things that are attached to locally-stored cells. 
  // For the moment, we do this by brute force.
  int num_faces = 0, num_edges = 0, num_nodes = 0, num_ghost_cells = 0;
  int cell_offset = 0;
  for (int b = low_block; b < high_block; ++b)
  {
    int N1 = resolutions[b].N1, N2 = resolutions[b].N2, N3 = resolutions[b].N3;
    cubic_lattice_t* lattice = lattices[b];

    // Go through each cell within this block that exists on this process and 
    // construct its stuff.
    int low_cell = (b == low_block) ? low_block_cell : 0;
    int high_cell = (b == (high_block - 1)) ? high_block_cell : num_block_cells[b];
    for (int c = low_cell; c < high_cell; ++c, ++cell_offset)
    {
      // Get the (i, j, k) position of this cell in the lattice.
      int i = c % N1, j = c / N1, k = c / (N1*N2);
      ASSERT(cubic_lattice_cell(lattice, i, j, k) == c);

      // Faces.
      int faces[6] = {cubic_lattice_x_face(lattice, i, j, k),
                      cubic_lattice_x_face(lattice, i+1, j, k),
                      cubic_lattice_y_face(lattice, i, j, k),
                      cubic_lattice_y_face(lattice, i, j+1, k),
                      cubic_lattice_z_face(lattice, i, j, k),
                      cubic_lattice_z_face(lattice, i, j, k+1)};
      for (int f = 0; f < 6; ++f)
      {
        int face_index = face_block_offsets[b] + faces[f];
        hexes[cell_offset].faces[f] = face_index;
        if (!int_unordered_set_contains(face_counted, face_index))
        {
          ++num_faces;
          int_unordered_set_insert(face_counted, face_index);
        }
      }

      // Edges.
      int edges[12] = {cubic_lattice_x_edge(lattice, i, j, k),
                       cubic_lattice_y_edge(lattice, i+1, j, k),
                       cubic_lattice_x_edge(lattice, i, j+1, k),
                       cubic_lattice_y_edge(lattice, i, j, k),
                       cubic_lattice_z_edge(lattice, i, j, k),
                       cubic_lattice_z_edge(lattice, i+1, j, k),
                       cubic_lattice_z_edge(lattice, i+1, j+1, k),
                       cubic_lattice_z_edge(lattice, i, j+1, k),
                       cubic_lattice_x_edge(lattice, i, j, k+1),
                       cubic_lattice_y_edge(lattice, i+1, j, k+1),
                       cubic_lattice_x_edge(lattice, i, j+1, k+1),
                       cubic_lattice_y_edge(lattice, i, j, k+1)};
      for (int e = 0; e < 12; ++e)
      {
        int edge_index = edge_block_offsets[b] + edges[e];
        hexes[cell_offset].edges[e] = edge_index;
        if (!int_unordered_set_contains(edge_counted, edge_index))
        {
          ++num_edges;
          int_unordered_set_insert(edge_counted, edge_index);
        }
      }

      // Nodes.
      int nodes[8] = {cubic_lattice_node(lattice, i, j, k),
                      cubic_lattice_node(lattice, i+1, j, k),
                      cubic_lattice_node(lattice, i+1, j+1, k),
                      cubic_lattice_node(lattice, i, j+1, k),
                      cubic_lattice_node(lattice, i, j, k+1),
                      cubic_lattice_node(lattice, i+1, j, k+1),
                      cubic_lattice_node(lattice, i+1, j+1, k+1),
                      cubic_lattice_node(lattice, i, j+1, k+1)};
      for (int n = 0; n < 8; ++n)
      {
        int node_index = node_block_offsets[b] + nodes[n];
        hexes[cell_offset].nodes[n] = node_index;
        if (!int_unordered_set_contains(node_counted, node_index))
        {
          ++num_nodes;
          int_unordered_set_insert(node_counted, node_index);
        }
      }

      // Ghost cells.
      if ((i == 0) || (i == (N1-1))) num_ghost_cells += num_ghosts;
      if ((j == 0) || (j == (N2-1))) num_ghost_cells += num_ghosts;
      if ((k == 0) || (k == (N3-1))) num_ghost_cells += num_ghosts;
    }
  }
  ASSERT(cell_offset == num_cells);

  // Construct a local mesh.
  mesh_t* mesh = mesh_new(num_cells, num_ghost_cells, num_faces, num_edges, num_nodes);

  // Create each of the blocks within the mesh for this process.
  for (int b = low_block; b < high_block; ++b)
  {
    block_t* block = block_assembly_block(assembly, b);
    if (block != NULL)
    {
      int low_cell = (b == low_block) ? low_block_cell : 0;
      int high_cell = (b == (high_block - 1)) ? high_block_cell : num_block_cells[b];
      create_block_grid(block, low_cell, high_cell, &hexes[cell_block_offsets[b]], mesh);
    }
  }

  // Tuck the index space into the mesh.
  mesh_set_property(mesh, "index_space", is, NULL);

  // Clean up.
  int_unordered_set_free(face_counted);
  int_unordered_set_free(edge_counted);
  int_unordered_set_free(node_counted);
  free(hexes);

  return mesh;
}

#ifdef __cplusplus
}
#endif

