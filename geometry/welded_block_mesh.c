// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/unordered_map.h"
#include "geometry/welded_block_mesh.h"
#include "geometry/cubic_lattice.h"

mesh_t* welded_block_mesh(mesh_t** blocks, int num_blocks)
{
#ifndef NDEBUG
  for (int i = 0; i < num_blocks; ++i)
  {
    // Currently, we can only construct a block-structured mesh on a single 
    // process.
    ASSERT(blocks[i]->comm == MPI_COMM_SELF); 
    ASSERT(blocks[i]->num_ghost_cells == 0); 

    // Make sure that we can access boundary tags.
    ASSERT(mesh_property(blocks[i], "rectilinear_boundary_tags") != NULL);
  }
#endif

  // Make an inventory.

  // The welded mesh will have a number of cells equal to the sum total of 
  // all our blocks. Faces and nodes, meanwhile, will be less than this sum.
  int num_cells = 0, num_faces = 0, num_nodes = 0;
//  int_int_unordered_map_t* face_indices = int_int_unordered_map_new();
//  int_int_unordered_map_t* node_indices = int_int_unordered_map_new();
  int block_neighbors[num_blocks][6];
  int block_neighbor_faces[num_blocks][6];
  for (int i = 0; i < num_blocks; ++i)
  {
    num_cells += blocks[i]->num_cells;
    num_faces += blocks[i]->num_faces;
    num_nodes += blocks[i]->num_nodes;
    string_array_t* btags1 = mesh_property(blocks[i], "rectilinear_boundary_tags");

    cubic_lattice_t* lattice1 = mesh_property(blocks[i], "cubic_lattice");
    uint64_t nxf1 = cubic_lattice_num_x_faces(lattice1);
    uint64_t nyf1 = cubic_lattice_num_y_faces(lattice1);
    uint64_t nzf1 = cubic_lattice_num_z_faces(lattice1);
    uint64_t nfaces1[6] = {nyf1 * nzf1, nyf1 * nzf1, 
                           nxf1 * nzf1, nxf1 * nzf1, 
                           nxf1 * nzf1, nyf1 * nyf1};
    uint64_t nnodes1[6] = {(nyf1+1) * (nzf1+1), (nyf1+1) * (nzf1+1),
                           (nxf1+1) * (nzf1+1), (nxf1+1) * (nzf1+1),
                           (nxf1+1) * (nyf1+1), (nxf1+1) * (nyf1+1)};
    // Match up our boundary tags with the other blocks.
    for (int j = i+1; j < num_blocks; ++j)
    {
      cubic_lattice_t* lattice2 = mesh_property(blocks[j], "cubic_lattice");
      uint64_t nxf2 = cubic_lattice_num_x_faces(lattice2);
      uint64_t nyf2 = cubic_lattice_num_y_faces(lattice2);
      uint64_t nzf2 = cubic_lattice_num_z_faces(lattice2);
      uint64_t nfaces2[6] = {nyf2 * nzf2, nyf2 * nzf2, 
                             nxf2 * nzf2, nxf2 * nzf2, 
                             nxf2 * nzf2, nyf2 * nyf2};
      string_array_t* btags2 = mesh_property(blocks[j], "rectilinear_boundary_tags");
      for (int k = 0; k < 6; ++k)
      {
        for (int l = 0; l < 6; ++l)
        {
          if (strcmp(btags1->data[k], btags2->data[l]) == 0)
          {
            if (nfaces1[k] != nfaces2[l])
            {
              polymec_error("welded_block_mesh: Cannot weld boundary face %d of block %d\n"
                            "  to boundary face %d of block %d: %d faces != %d faces.", 
                            k, i, l, j, nfaces1[k], nfaces2[l]);
            }
            num_faces -= (int)nfaces1[k];
            num_nodes -= (int)nnodes1[k];
            block_neighbors[i][k] = j;
            block_neighbor_faces[i][k] = l;
          }
        }
      }
    }
  }
  mesh_t* block_mesh = mesh_new(MPI_COMM_SELF, num_cells, 0, 
                                num_faces, num_nodes);

  return block_mesh;
}

