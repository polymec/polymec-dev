#include "geometry/create_block_mesh.h"
#include "core/index_space.h"
#include "geometry/cubic_lattice.h"

#ifdef __cplusplus
extern "C" {
#endif

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

  // Figure out the numbers of ghost cells, local faces, edges, nodes.
  int num_ghost_cells, num_faces, num_edges, num_nodes;
  mesh_t* mesh = mesh_new(num_cells, num_ghost_cells, num_faces, num_edges, num_nodes);

  // Tuck the index space into the mesh.
  mesh_set_property(mesh, "index_space", is, NULL);

  return NULL;
}

#ifdef __cplusplus
}
#endif

