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

mesh_t* welded_block_mesh(mesh_t** blocks, int num_blocks, real_t weld_tolerance)
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

  // Go through all the blocks and put them into kd-trees.
  kd_tree_t *face_trees[num_blocks], *node_trees[num_blocks];
  for (int i = 0; i < num_blocks; ++i)
  {
    face_trees[i] = kd_tree_new(blocks[i]->face_centers, blocks[i]->num_faces);
    node_trees[i] = kd_tree_new(blocks[i]->nodes, blocks[i]->num_nodes);
  }

  // Now weld together groups of faces/nodes in the same location.
  int_tuple_int_unordered_map_t* face_welds = int_tuple_int_unordered_map_new();
  int_tuple_int_unordered_map_t* node_welds = int_tuple_int_unordered_map_new();
  int total_num_cells = 0, total_num_faces = 0, total_num_nodes = 0, 
      num_face_welds = 0, num_node_welds = 0;
  for (int i = 0; i < num_blocks; ++i)
  {
    mesh_t* blocki = blocks[i];
    total_num_cells += blocki->num_cells;
    total_num_faces += blocki->num_faces;
    total_num_nodes += blocki->num_nodes;
    string_array_t* btagsi = mesh_property(blocks[i], "rectilinear_boundary_tags");

    for (int j = i+1; j < num_blocks; ++j)
    {
      mesh_t* blockj = blocks[j];
      string_array_t* btagsj = mesh_property(blocks[i], "rectilinear_boundary_tags");

      // Find the face tags common to blocks i and j.
      for (int ii = 0; ii < btagsi->size; ++ii)
      {
        char* btagi = btagsi->data[ii];
        for (int jj = 0; jj < btagsj->size; ++jj)
        {
          char* btagj = btagsi->data[ii];
          if (strcmp(btagi, btagj) == 0)
          {
            // For each face in btagi, find the closest face in btagj.
            int num_btagi_faces;
            int* btagi_faces = mesh_tag(blocki->face_tags, btagi, &num_btagi_faces);
            for (int f = 0; f < num_btagi_faces; ++f)
            {
              // Weld'em faces.
              int face = btagi_faces[f];
              point_t* xf = &(blocki->face_centers[face]);
              int nearest = kd_tree_nearest(face_trees[j], xf);
              ASSERT(nearest != -1);
              point_t* xfn = &(blockj->face_centers[nearest]);
              if (point_distance(xf, xfn) > weld_tolerance)
              {
                polymec_error("No match found for face %d of block %d within a distance of %g.\n"
                              "(Boundary tag: %s)", face, i, btagi);
              }

              int* tuplei = int_tuple_new(2);
              tuplei[0] = i;
              tuplei[1] = face;
              int* tuplej = int_tuple_new(2);
              tuplej[0] = j;
              tuplej[1] = nearest;
              int_tuple_int_unordered_map_insert_with_k_dtor(face_welds, tuplei, num_face_welds, int_tuple_free);
              int_tuple_int_unordered_map_insert_with_k_dtor(face_welds, tuplej, num_face_welds, int_tuple_free);
              ++num_face_welds;

              // Weld'em nodes.
              int pos = 0, node;
              while (mesh_face_next_node(blocki, face, &pos, &node))
              {
                point_t* xn = &(blocki->nodes[node]);
                int nearest = kd_tree_nearest(node_trees[j], xn);
                ASSERT(nearest != -1);
                point_t* xnn = &(blockj->nodes[nearest]);
                if (point_distance(xn, xnn) > weld_tolerance)
                {
                  polymec_error("No match found for node %d of block %d within a distance of %g.\n"
                                "(Boundary tag: %s)", node, i, btagi);
                }

                int* tuplei = int_tuple_new(2);
                tuplei[0] = i;
                tuplei[1] = node;
                int* tuplej = int_tuple_new(2);
                tuplej[0] = j;
                tuplej[1] = nearest;
                int_tuple_int_unordered_map_insert_with_k_dtor(node_welds, tuplei, num_node_welds, int_tuple_free);
                int_tuple_int_unordered_map_insert_with_k_dtor(node_welds, tuplej, num_node_welds, int_tuple_free);
                ++num_node_welds;
              }
            }
          }
        }
      }
    }
  }

  // Chop down'em trees.
  for (int i = 0; i < num_blocks; ++i)
  {
    kd_tree_free(face_trees[i]);
    kd_tree_free(node_trees[i]);
  }

  // Now set up the welded block mesh, which consists of hexahedra, which have 
  // 6 faces per cell and 4 nodes per face.
  mesh_t* block_mesh = mesh_new_with_cell_type(MPI_COMM_SELF, total_num_cells, 0, 
                                               total_num_faces - num_face_welds, 
                                               total_num_nodes - num_node_welds,
                                               6, 4);

  // We'll use these arrays to help keep books when indexing welded faces/nodes.
  int welded_face_indices[num_face_welds], welded_node_indices[num_node_welds];
  for (int i = 0; i < num_face_welds; ++i)
    welded_face_indices[i] = -1;
  for (int i = 0; i < num_node_welds; ++i)
    welded_node_indices[i] = -1;

  // Hook everything up.
  int next_face_index = 0;
  int* tuple = int_tuple_new(2);
  for (int i = 0; i < num_blocks; ++i)
  {
    mesh_t* block = blocks[i];
    for (int cell = 0; cell < block->num_cells; ++cell)
    {
      int pos = 0, face;
      while (mesh_cell_next_oriented_face(block, cell, &pos, &face))
      {
        int which_face = pos - 1;
        int actual_face = (face < 0) ? ~face : face;

        // Has this face already been constructed? We orient the face so
        // that its normal points outward w.r.t. its original cell.
        tuple[0] = i;
        tuple[1] = actual_face;
        int* weld_p = int_tuple_int_unordered_map_get(face_welds, tuple);
        if (weld_p == NULL)
        {
          int face_index = next_face_index;
          block_mesh->cell_faces[which_face] = face_index;
          ++next_face_index;
        }
        else
        {
          int weld_index = *weld_p;
          if (welded_face_indices[weld_index] == -1)
          {
            welded_face_indices[weld_index] = next_face_index;
            ++next_face_index;
          }
          int face_index = welded_face_indices[weld_index];
          int oriented_face = (face < 0) ? ~face_index : face_index;
          block_mesh->cell_faces[which_face] = oriented_face;
        }

        // Now set up the nodes.
        int npos = 0, node;
        while (mesh_face_next_node(block, actual_face, &npos, &node))
        {
          int which_node = pos - 1;

        }
      }
    }
  }

  // Create edges, compute geometry.
  mesh_construct_edges(block_mesh);
  mesh_compute_geometry(block_mesh);

  // Clean up.
  int_tuple_free(tuple);
  int_tuple_int_unordered_map_free(face_welds);
  int_tuple_int_unordered_map_free(node_welds);

  return block_mesh;
}

