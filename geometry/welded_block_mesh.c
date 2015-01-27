// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/unordered_map.h"
#include "core/kd_tree.h"
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

    // Make sure that we can access boundary tags. If not, blocks[i] is 
    // *not* a rectilinear mesh.
    ASSERT(mesh_property(blocks[i], "rectilinear_boundary_tags") != NULL);
  }
#endif

  // Make an inventory.

  // The welded mesh will have a number of cells equal to the sum total of 
  // all our blocks. Faces and nodes, meanwhile, will be less than this sum.
  int num_cells = 0, num_faces = 0, num_nodes = 0;
  int block_neighbors[num_blocks][6];
  int block_neighbor_faces[num_blocks][6];
  for (int i = 0; i < num_blocks; ++i)
  {
    for (int j = 0; j < num_blocks; ++j)
    {
      block_neighbors[i][j] = -1;
      block_neighbor_faces[i][j] = -1;
    }
    num_cells += blocks[i]->num_cells;
    num_faces += blocks[i]->num_faces;
    num_nodes += blocks[i]->num_nodes;
    string_array_t* btagsi = mesh_property(blocks[i], "rectilinear_boundary_tags");

    cubic_lattice_t* latticei = mesh_property(blocks[i], "cubic_lattice");
    ASSERT(latticei != NULL);
    index_t nxfi = cubic_lattice_num_x_faces(latticei);
    index_t nyfi = cubic_lattice_num_y_faces(latticei);
    index_t nzfi = cubic_lattice_num_z_faces(latticei);
    index_t nfacesi[6] = {nyfi * nzfi, nyfi * nzfi, 
                           nxfi * nzfi, nxfi * nzfi, 
                           nxfi * nzfi, nyfi * nyfi};
    index_t nnodesi[6] = {(nyfi+1) * (nzfi+1), (nyfi+1) * (nzfi+1),
                           (nxfi+1) * (nzfi+1), (nxfi+1) * (nzfi+1),
                           (nxfi+1) * (nyfi+1), (nxfi+1) * (nyfi+1)};
    // Match up our boundary tags with the other blocks.
    for (int j = i+1; j < num_blocks; ++j)
    {
      cubic_lattice_t* latticej = mesh_property(blocks[j], "cubic_lattice");
      ASSERT(latticej != NULL);
      index_t nxfj = cubic_lattice_num_x_faces(latticej);
      index_t nyfj = cubic_lattice_num_y_faces(latticej);
      index_t nzfj = cubic_lattice_num_z_faces(latticej);
      index_t nfacesj[6] = {nyfj * nzfj, nyfj * nzfj, 
                            nxfj * nzfj, nxfj * nzfj, 
                            nxfj * nzfj, nyfj * nyfj};
      string_array_t* btagsj = mesh_property(blocks[j], "rectilinear_boundary_tags");
      for (int k = 0; k < 6; ++k)
      {
        for (int l = 0; l < 6; ++l)
        {
          if (strcmp(btagsi->data[k], btagsj->data[l]) == 0)
          {
            if (nfacesi[k] != nfacesj[l])
            {
              polymec_error("welded_block_mesh: Cannot weld boundary face %d of block %d\n"
                            "  to boundary face %d of block %d: %d faces != %d faces.", 
                            k, i, l, j, nfacesi[k], nfacesj[l]);
            }
            num_faces -= (int)nfacesi[k];
            num_nodes -= (int)nnodesi[k];
            block_neighbors[i][k] = j;
            block_neighbor_faces[i][k] = l;
          }
        }
      }
    }
  }

  // Now set up the welded block mesh.
  mesh_t* block_mesh = mesh_new(MPI_COMM_SELF, num_cells, 0, 
                                num_faces, num_nodes);

  // Keep track of faces in the blocks that we've created.
  kd_tree_t* block_trees[num_blocks];
  for (int b = 0; b < num_blocks; ++b)
    block_trees[b] = NULL;

  // Go over each block and set up cells, nodes, etc. Keep track of faces 
  // and nodes that have been created by other blocks.
  index_t cell_offset = 0;
  for (int b = 0; b < num_blocks; ++b)
  {
    block_trees[b] = kd_tree_new(blocks[b]->face_centers, blocks[b]->num_faces);
    cubic_lattice_t* lattice = mesh_property(blocks[b], "cubic_lattice");

    for (int cell = 0; cell < blocks[b]->num_cells; ++c)
    {
      mesh_t* bmesh = blocks[b];
      index_t bcell = cell_offset + cell;
      index_t i, j, k;
      cubic_lattice_get_cell_triple(lattice, cell, &i, &j, &k);
          
      // Handle boundary faces.
      if ((i == 0) || (i == (nx-1)) || 
          (j == 0) || (j == (ny-1)) || 
          (k == 0) || (k == (nz­1)))
      {
        for (int n = 0; n < 6; ++n)
        {
          int nblock = block_neighbors[b][n];
          if (block_trees[nblock] != NULL)
          { 
            // Match the nearest pair of faces at this seam.
            int nearest = -1, pos = 0, face;
            real_t D = FLT_MAX;
            while (mesh_cell_next_face(bmesh, cell, &pos, &face))
            {
              point_t* xf = &bmesh->face_centers[face]; 
              int p = kd_tree_nearest(block_trees[nblock], xf);
              point_t* xp = &blocks[nblock]->face_centers[p];
              real_t Dp = point_distance(xf, xp);
              if (Dp < D)
              {
                nearest = p;
                D = Dp;
              }
            }
        }
          else
          {
          }
        }
      }
    }
    cell_offset += (index_t)blocks[b]->num_cells;
  }

  return block_mesh;
}

