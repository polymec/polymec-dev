#include "geometry/create_block_mesh.h"
#include "core/index_space.h"
#include "core/unordered_set.h"
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

mesh_t* create_block_mesh(MPI_Comm comm, 
                          block_assembly_t* assembly, 
                          block_grid_t* resolutions, 
                          int num_ghosts)
{
  // Compute the number of cells in the assembly.
  int num_blocks = block_assembly_num_blocks(assembly), num_global_cells = 0;
  int num_block_cells[num_blocks];
  for (int i = 0; i < num_blocks; ++i)
  {
    block_t* block = block_assembly_block(assembly, i);
    if (block != NULL)
    {
      num_block_cells[i] = resolutions[i].N1*resolutions[i].N2*resolutions[i].N3;
      num_global_cells += num_block_cells[i];
    }
    else
      num_block_cells[i] = 0;
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
    cubic_lattice_t* lattice = cubic_lattice_new(N1, N2, N3);

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
      memcpy(hexes[cell_offset].faces, faces, 6*sizeof(int));

      for (int f = 0; f < 6; ++f)
      {
        if (!int_unordered_set_contains(face_counted, faces[f]))
        {
          ++num_faces;
          int_unordered_set_insert(face_counted, faces[f]);
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
      memcpy(hexes[cell_offset].edges, edges, 12*sizeof(int));

      for (int e = 0; e < 12; ++e)
      {
        if (!int_unordered_set_contains(edge_counted, edges[e]))
        {
          ++num_edges;
          int_unordered_set_insert(edge_counted, edges[e]);
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
      memcpy(hexes[cell_offset].nodes, nodes, 8*sizeof(int));

      for (int n = 0; n < 8; ++n)
      {
        if (!int_unordered_set_contains(node_counted, nodes[n]))
        {
          ++num_nodes;
          int_unordered_set_insert(node_counted, nodes[n]);
        }
      }

      // Ghost cells.
      if ((i == 0) || (i == (N1-1))) num_ghost_cells += num_ghosts;
      if ((j == 0) || (j == (N2-1))) num_ghost_cells += num_ghosts;
      if ((k == 0) || (k == (N3-1))) num_ghost_cells += num_ghosts;
    }
  }

  // Construct a local mesh.
  mesh_t* mesh = mesh_new(num_cells, num_ghost_cells, num_faces, num_edges, num_nodes);

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

