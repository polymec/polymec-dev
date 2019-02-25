// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "geometry/partition_polymesh.h"

#if POLYMEC_HAVE_MPI
#include "core/array_utils.h"
#include "core/unordered_set.h"
#include "core/kd_tree.h"
#include "core/hilbert.h"
#include "core/timer.h"
#include "core/partitioning.h"

// This creates tags for a submesh whose elements belong to the given set.
static void create_submesh_tags(tagger_t* tagger, 
                                int_int_unordered_map_t* element_map)
{
  int pos = 0, *indices;
  size_t size;
  char* tag_name;
  int_array_t* sub_indices = int_array_new();
  while (tagger_next_tag(tagger, &pos, &tag_name, &indices, &size))
  {
    int* tag = tagger_create_tag(tagger, tag_name, sub_indices->size);
    if (tag != NULL) // Skip tags that already exist!
    {
      for (int i = 0; i < size; ++i)
      {
        int* elem_p = int_int_unordered_map_get(element_map, indices[i]);
        if (elem_p != NULL)
          int_array_append(sub_indices, *elem_p);
      }
      memcpy(tag, sub_indices->data, sizeof(int) * sub_indices->size);
    }
  }
  int_array_free(sub_indices);
}

// This helper constructs and returns a mesh from the cells with the given
// indices in the given mesh. The submesh is valid with the following 
// exception:
// The indices of ghost cells referenced in submesh->face_cells are 
// replaced with destination process ranks.
static polymesh_t* create_submesh(MPI_Comm comm, polymesh_t* mesh, 
                                  int64_t* partition, index_t* vtx_dist, 
                                  int* indices, size_t num_indices)
{
  START_FUNCTION_TIMER();

  // Figure out the destination process for this submesh.
  int dest_proc;
  MPI_Comm_rank(comm, &dest_proc);
  if (num_indices > 0)
    dest_proc = (int)partition[indices[0]];

  // Make a set of cells for querying membership in this submesh.
  int_unordered_set_t* cell_set = int_unordered_set_new();
  for (int i = 0; i < num_indices; ++i)
  {
    ASSERT(partition[indices[i]] == dest_proc);
    int_unordered_set_insert(cell_set, indices[i]);
  }

  // Count unique mesh elements (faces, nodes).
  int num_cells = (int)num_indices, num_ghost_cells = 0;
  int_unordered_set_t* face_indices = int_unordered_set_new();
  int_unordered_set_t* node_indices = int_unordered_set_new();
  for (int i = 0; i < num_indices; ++i)
  {
    int cell = indices[i], fpos = 0, face;
    while (polymesh_cell_next_face(mesh, cell, &fpos, &face))
    {
      int opp_cell = polymesh_face_opp_cell(mesh, face, cell);
      if ((opp_cell != -1) && !int_unordered_set_contains(cell_set, opp_cell))
      {
        int_unordered_set_insert(cell_set, opp_cell);
        ++num_ghost_cells;
      }

      int_unordered_set_insert(face_indices, face);
      int npos = 0, node;
      while (polymesh_face_next_node(mesh, face, &npos, &node))
        int_unordered_set_insert(node_indices, node);
    }
  }
  int_unordered_set_free(cell_set);

  // Create the submesh container.
  int num_faces = face_indices->size, num_nodes = node_indices->size;
  log_debug("create_submesh: Creating polymesh (%d cells, %d faces, %d nodes)", 
            num_cells, num_faces, num_nodes);
  polymesh_t* submesh = polymesh_new(comm, num_cells, num_ghost_cells, num_faces, num_nodes);

  // Create mappings for faces and nodes.
  // face_map maps a submesh face to its corresponding mesh face.
  int* face_map = polymec_malloc(sizeof(int) * num_faces);
  {
    int fpos = 0, face, f = 0;
    while (int_unordered_set_next(face_indices, &fpos, &face))
      face_map[f++] = face;
  }
  int_unordered_set_free(face_indices);

  // node_map maps a submesh node to its corresponding mesh node.
  int* node_map = polymec_malloc(sizeof(int) * num_nodes);
  {
    int npos = 0, node, n = 0;
    while (int_unordered_set_next(node_indices, &npos, &node))
      node_map[n++] = node;
  }
  int_unordered_set_free(node_indices);

  // Allocate cell faces and face nodes (and generate inverse maps).
  // inverse_cell_map maps a mesh cell to its corresponding submesh cell.
  int_int_unordered_map_t* inverse_cell_map = int_int_unordered_map_new();
  submesh->cell_face_offsets[0] = 0;
  for (int c = 0; c < submesh->num_cells; ++c)
  {
    int_int_unordered_map_insert(inverse_cell_map, indices[c], c);
    int num_cell_faces = mesh->cell_face_offsets[indices[c]+1] - mesh->cell_face_offsets[indices[c]];
    submesh->cell_face_offsets[c+1] = submesh->cell_face_offsets[c] + num_cell_faces;
  }

  // inverse_face_map maps a mesh face to its corresponding submesh face.
  int_int_unordered_map_t* inverse_face_map = int_int_unordered_map_new();
  submesh->face_node_offsets[0] = 0;
  for (int f = 0; f < submesh->num_faces; ++f)
  {
    int_int_unordered_map_insert(inverse_face_map, face_map[f], f);
    int num_face_nodes = mesh->face_node_offsets[face_map[f]+1] - mesh->face_node_offsets[face_map[f]];
    submesh->face_node_offsets[f+1] = submesh->face_node_offsets[f] + num_face_nodes;
  }

  // inverse_node_map maps a mesh node to its corresponding submesh node.
  int_int_unordered_map_t* inverse_node_map = int_int_unordered_map_new();
  for (int n = 0; n < submesh->num_nodes; ++n)
    int_int_unordered_map_insert(inverse_node_map, node_map[n], n);
  polymesh_reserve_connectivity_storage(submesh);

  // Copy cell faces.
  for (int c = 0; c < submesh->num_cells; ++c)
  {
    int num_cell_faces = submesh->cell_face_offsets[c+1] - submesh->cell_face_offsets[c];
    for (int f = 0; f < num_cell_faces; ++f)
    {
      int face_index = mesh->cell_faces[mesh->cell_face_offsets[indices[c]]+f];
      bool flipped = false;
      if (face_index < 0)
      {
        flipped = true;
        face_index = ~face_index;
      }
      int subface_index = *int_int_unordered_map_get(inverse_face_map, face_index);
      if (flipped)
        subface_index = ~subface_index;
      submesh->cell_faces[submesh->cell_face_offsets[c]+f] = subface_index;
    }
  }

  // Create a mapping of parallel boundary faces to global "ghost" cells.
  // We will use this to create an annotated "parallel_boundary_faces" tag 
  // in the submesh.
  int_int_unordered_map_t* parallel_bface_map = int_int_unordered_map_new();

  // Copy face cells and face nodes.
  for (int f = 0; f < submesh->num_faces; ++f)
  {
    // Identify the cells attached to the face and construct them within the submesh.
    int orig_mesh_face = face_map[f];
    int* cell_p = int_int_unordered_map_get(inverse_cell_map, mesh->face_cells[2*orig_mesh_face]);
    int* opp_cell_p = int_int_unordered_map_get(inverse_cell_map, mesh->face_cells[2*orig_mesh_face+1]);
    ASSERT((cell_p != NULL) || (opp_cell_p != NULL)); // There must be at least one internal cell attached to this face.
    int this_cell = -1, that_cell = -1;
    if ((cell_p != NULL) && (opp_cell_p != NULL))
    {
      this_cell = *cell_p;
      that_cell = *opp_cell_p;
    }
    else
    {
      // One of the cells attached to this face is a ghost cell, at least as far as this
      // submesh is concerned. Record the index of the "ghost cell" in the 
      // original mesh.
      int ghost_cell;
      if (cell_p != NULL)
      {
        this_cell = *cell_p;
        ghost_cell = mesh->face_cells[2*orig_mesh_face+1];
      }
      else
      {
        this_cell = *opp_cell_p;
        ghost_cell = mesh->face_cells[2*orig_mesh_face];
      }

      // We encode the destination process rank in the ghost cell and stash 
      // the original ghost index in our parallel boundary face map. The 
      // destination rank of the ghost cell is stored in its corresponding 
      // entry within our partition vector.
      if (ghost_cell != -1)
      {
        int ghost_cell_proc = (int)partition[ghost_cell];
        that_cell = -ghost_cell_proc - 2; // encoding of dest proc
        int_int_unordered_map_insert(parallel_bface_map, f, ghost_cell);
      }
    }
    submesh->face_cells[2*f] = this_cell;
    submesh->face_cells[2*f+1] = that_cell;

    // Copy over the nodes for this face.
    int num_face_nodes = submesh->face_node_offsets[f+1] - submesh->face_node_offsets[f];
    for (int n = 0; n < num_face_nodes; ++n)
    {
      int subnode = *int_int_unordered_map_get(inverse_node_map, mesh->face_nodes[mesh->face_node_offsets[orig_mesh_face]+n]);
      submesh->face_nodes[submesh->face_node_offsets[f]+n] = subnode;
    }
  }

  // Create a tag for boundary faces.
  {
    int* pbf_tag = polymesh_create_tag(submesh->face_tags, "parallel_boundary_faces", 
                                       parallel_bface_map->size);
    int_array_t* pbf_ghost_cells = int_array_new();
    int pos = 0, face, gcell, i = 0;
    while (int_int_unordered_map_next(parallel_bface_map, &pos, &face, &gcell))
    {
      pbf_tag[i] = face;
      int_array_append(pbf_ghost_cells, gcell);
      ++i;
    }

    // Stash the ghost cell indices in another tag.
    int* gci_tag = polymesh_create_tag(submesh->cell_tags, "ghost_cell_indices", 
                                       pbf_ghost_cells->size);
    memcpy(gci_tag, pbf_ghost_cells->data, sizeof(int) * pbf_ghost_cells->size);
    int_array_free(pbf_ghost_cells);

    // Also set up a copy of the array of global cell indices.
    gci_tag = polymesh_create_tag(submesh->cell_tags, "global_cell_indices", 
                                  num_indices);
    memcpy(gci_tag, indices, sizeof(int) * num_indices);
  }
  int_int_unordered_map_free(parallel_bface_map);

  // Copy node positions.
  for (int n = 0; n < submesh->num_nodes; ++n)
  {
    int orig_mesh_node = node_map[n];
    submesh->nodes[n] = mesh->nodes[orig_mesh_node];
  }

  // Construct edges.
  polymesh_construct_edges(submesh);

  // Do geometry.
  polymesh_compute_geometry(submesh);

  // Dump any tags in the existing global mesh in, too.
  // FIXME: Edge tags are not supported!
  create_submesh_tags(submesh->cell_tags, inverse_cell_map);
  create_submesh_tags(submesh->face_tags, inverse_face_map);
  create_submesh_tags(submesh->node_tags, inverse_node_map);

  // Clean up.
  int_int_unordered_map_free(inverse_cell_map);
  int_int_unordered_map_free(inverse_face_map);
  int_int_unordered_map_free(inverse_node_map);
  polymec_free(face_map);
  polymec_free(node_map);

  STOP_FUNCTION_TIMER();
  return submesh;
}

static void sort_global_cell_pairs(int* indices, int num_pairs)
{
  START_FUNCTION_TIMER();
  int data[3*num_pairs]; // (i, j, swapped) for each pair
  for (int i = 0; i < num_pairs; ++i)
  {
    if (indices[2*i] > indices[2*i+1])
    {
      data[3*i]   = indices[2*i+1];
      data[3*i+1] = indices[2*i];
      data[3*i+2] = 1;
    }
    else
    {
      data[3*i]   = indices[2*i];
      data[3*i+1] = indices[2*i+1];
      data[3*i+2] = 0;
    }
  }
  // We can sort using the ordinary integer pair comparator, since the 
  // swapped flag doesn't affect the ordering.
  qsort(data, (size_t)num_pairs, 3*sizeof(int), int_pair_bsearch_comp);
  for (int i = 0; i < num_pairs; ++i)
  {
    if (data[3*i+2] == 1) // swapped
    {
      indices[2*i]   = data[3*i+1];
      indices[2*i+1] = data[3*i];
    }
    else
    {
      indices[2*i]   = data[3*i];
      indices[2*i+1] = data[3*i+1];
    }
  }
  STOP_FUNCTION_TIMER();
}

// This helper creates an array containing an index map that removes "holes" 
// left by mapped duplicates as mapped in dup_map. The range [0, num_indices) 
// is mapped. This is used to remap faces and nodes in the fuse_submeshes() 
// helper below. A newly-allocated array of length num_indices is returned.
static int* create_index_map_with_dups_removed(int num_indices, 
                                               int_int_unordered_map_t* dup_map)
{
  START_FUNCTION_TIMER();

  // Count the number of holes before a given index.
  int num_holes[num_indices];
  num_holes[0] = 0;
  ASSERT(!int_int_unordered_map_contains(dup_map, 0));
  for (int i = 1; i < num_indices; ++i)
  {
    num_holes[i] = num_holes[i-1];
    if (int_int_unordered_map_contains(dup_map, i))
      ++num_holes[i];
  }

  int* map = polymec_malloc(sizeof(int) * num_indices);
  for (int i = 0; i < num_indices; ++i)
  {
    int* k_ptr = int_int_unordered_map_get(dup_map, i);
    if (k_ptr == NULL)
      map[i] = i - num_holes[i];
    else
    {
      // Follow any chain till we hit the end.
      int index = 0;
      while (k_ptr != NULL)
      {
        index = *k_ptr;
        k_ptr = int_int_unordered_map_get(dup_map, index);
      }
      map[i] = index - num_holes[index];
    }
    ASSERT(map[i] >= 0);
    ASSERT(map[i] < (num_indices - dup_map->size));
  }

  STOP_FUNCTION_TIMER();
  return map;
}

// This helper takes an array of submeshes and stitches them all together into 
// one single mesh (contiguous or not) on the current domain. The submeshes are 
// consumed in the process.
static polymesh_t* fuse_submeshes(polymesh_t** submeshes, 
                                  size_t num_submeshes)
{
  START_FUNCTION_TIMER();
  ASSERT(num_submeshes > 0);
  MPI_Comm comm = submeshes[0]->comm;
  int rank, nprocs;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &nprocs);

  // First we traverse each of the submeshes and count up all the internal 
  // and ghost cells, faces, nodes.
  int num_cells = 0, num_ghost_cells = 0, num_faces = 0, num_nodes = 0, 
      submesh_cell_offsets[num_submeshes+1], submesh_face_offsets[num_submeshes+1], 
      submesh_node_offsets[num_submeshes+1];
  submesh_cell_offsets[0] = submesh_face_offsets[0] = submesh_node_offsets[0] = 0;
  for (int m = 0; m < num_submeshes; ++m)
  {
    polymesh_t* submesh = submeshes[m];
    num_cells += submesh->num_cells;
    num_ghost_cells += submesh->num_ghost_cells;
    submesh_cell_offsets[m+1] = num_cells;
    num_faces += submesh->num_faces;
    submesh_face_offsets[m+1] = num_faces;
    num_nodes += submesh->num_nodes;
    submesh_node_offsets[m+1] = num_nodes;
  }

  // Next, we traverse these submeshes and construct sets of all the 
  // faces/nodes that make up the "seams" of the fused mesh--those whose 
  // would-be-ghost cells belong to other submeshes on this process. 
  int_unordered_set_t* seam_faces = int_unordered_set_new();
  int_unordered_set_t* seam_nodes = int_unordered_set_new();
  for (int m = 0; m < num_submeshes; ++m)
  {
    polymesh_t* submesh = submeshes[m];
    for (int f = 0; f < submesh->num_faces; ++f)
    {
      if (submesh->face_cells[2*f+1] == (-rank - 2)) // belongs to this domain
      {
        int_unordered_set_insert(seam_faces, submesh_face_offsets[m] + f);

        // These nodes need to be merged with their neighbors.
        int pos = 0, node;
        while (polymesh_face_next_node(submesh, f, &pos, &node))
          int_unordered_set_insert(seam_nodes, submesh_node_offsets[m] + node);
      }
      else if (submesh->face_cells[2*f+1] <= -1) // belongs to another domain
      {
        // These nodes may need to be merged with their neighbors.
        int pos = 0, node;
        while (polymesh_face_next_node(submesh, f, &pos, &node))
          int_unordered_set_insert(seam_nodes, submesh_node_offsets[m] + node);
      }
    }
  }
  ASSERT((seam_faces->size % 2) == 0); // Seam faces always come in pairs!
  num_ghost_cells -= seam_faces->size;

  // Construct mappings to remove duplicate faces by mapping redundant 
  // ones to "originals."
  int_int_unordered_map_t* dup_face_map = int_int_unordered_map_new();
  if (seam_faces->size > 0)
  {
    // Translate the seam faces to an array.
    int seam_face_array[seam_faces->size];
    {
      int pos = 0, i = 0, j;
      while (int_unordered_set_next(seam_faces, &pos, &j))
        seam_face_array[i++] = j;
    }

    // Make a kd-tree of seam face centers.
    point_t xf[seam_faces->size];
    for (int i = 0; i < seam_faces->size; ++i)
    {
      // Find the submesh that contains this face.
      int m = 0, j = seam_face_array[i];
      while (j >= submesh_face_offsets[m+1]) ++m;

      // Put the face center into this array.
      int f = j - submesh_face_offsets[m];
      ASSERT((f >= 0) && (f < submeshes[m]->num_faces));
      xf[i] = submeshes[m]->face_centers[f];
    }
    kd_tree_t* face_tree = kd_tree_new(xf, seam_faces->size);

    // Now examine each seam face and find the 2 nearest faces to it.
    // One will be the face itself, and the other will be a duplicate. 
    // Map the one with the higher index to the lower index. Note that 
    // the indices of the points in the kd-tree should correspond to those
    // in the seam face array.
    for (int i = 0; i < seam_faces->size; ++i)
    {
      int nearest[2], j = seam_face_array[i];

      // Get the face center we're considering (in the same way we did above).
      int m = 0;
      while (j >= submesh_face_offsets[m+1]) ++m;
      int f = j - submesh_face_offsets[m];
      point_t* x = &submeshes[m]->face_centers[f];

      // Find the 2 nearest faces, one of which is j itself, and the other 
      // which ostensibly lies on top of it.
      kd_tree_nearest_n(face_tree, x, 2, nearest);
      int index0 = seam_face_array[nearest[0]];
      int index1 = seam_face_array[nearest[1]];
      ASSERT(index0 != index1);
      ASSERT((index0 == j) || (index1 == j));
#ifndef NDEBUG
      static real_t epsilon = 1e-12;
      ASSERT(point_distance(&xf[nearest[0]], &xf[nearest[1]]) < epsilon);
#endif

      // Merge the faces by mapping the one with the higher index to the 
      // lower index. 
      int_int_unordered_map_insert(dup_face_map, MAX(index0, index1), MIN(index0, index1));

      // ALSO: Alter the submesh so that its interior ghost cells actually 
      // refer to the correct (flattened) interior cells. This may seem strange,
      // but it helps us to get the face->cell connectivity right later on.
      int j1 = (j == index0) ? index1 : index0;
      int m1 = 0;
      while (j1 >= submesh_face_offsets[m1+1]) ++m1;
      int f1 = j1 - submesh_face_offsets[m1];
      ASSERT((f1 >= 0) && (f1 < submeshes[m1]->num_faces));
      int flattened_cell1 = submesh_cell_offsets[m1] + submeshes[m1]->face_cells[2*f1];
      submeshes[m]->face_cells[2*f+1] = flattened_cell1;
    }

    // Clean up.
    kd_tree_free(face_tree);
  }
  int* face_map = create_index_map_with_dups_removed(num_faces, dup_face_map); 

  // Do the same for duplicate nodes.
  int_int_unordered_map_t* dup_node_map = int_int_unordered_map_new();
  if (seam_nodes->size > 0)
  {
    // Translate the seam nodes to an array.
    int seam_node_array[seam_nodes->size];
    {
      int pos = 0, i = 0, j;
      while (int_unordered_set_next(seam_nodes, &pos, &j))
        seam_node_array[i++] = j;
    }

    // Make a kd-tree of node positions.
    point_t xn[seam_nodes->size];
    for (int i = 0; i < seam_nodes->size; ++i)
    {
      // Find the submesh that contains this face.
      int m = 0, j = seam_node_array[i];
      while (j >= submesh_node_offsets[m+1]) ++m;

      // Put the node position into this array.
      int n = j - submesh_node_offsets[m];
      xn[i] = submeshes[m]->nodes[n];
    }
    kd_tree_t* node_tree = kd_tree_new(xn, seam_nodes->size);

    // Now examine each seam node and find ALL NODES that appear to be 
    // the same. NOTE: for now, we find all the nodes within a distance 
    // of epsilon, which is KNOWN NOT TO BE A ROBUST METHOD OF FUSING NODES.
    // It might be good enough for now, though...
    for (int i = 0; i < seam_nodes->size; ++i)
    {
      int j = seam_node_array[i];

      // Get the node position we're considering (in the same way we did above).
      int m = 0;
      while (j >= submesh_node_offsets[m+1]) ++m;
      int n = j - submesh_node_offsets[m];
      point_t* x = &submeshes[m]->nodes[n];

      // Find all the nodes within epsilon of this one.
      // Note that not all "seam nodes" will be merged with others, since we 
      // identify seam node candidates from a broader class of nodes including
      // nodes on the problem boundary that don't necessary belong to a seam 
      // face.
      static const real_t epsilon = 1e-12;
      int_array_t* same_nodes = kd_tree_within_radius(node_tree, x, epsilon);

      if (same_nodes->size > 1)
      {
        // Merge the nodes by mapping all those with higher indices to the 
        // lowest index. Recall that we have to map the "seam index" of the point
        // (which the kd tree uses) to the actual node index.
        int min_index = INT_MAX;
        for (int k = 0; k < same_nodes->size; ++k)
        {
          int seam_index = same_nodes->data[k];
          int index = seam_node_array[seam_index];
          if (index < min_index)
            min_index = index;
        }
        for (int k = 0; k < same_nodes->size; ++k)
        {
          int seam_index = same_nodes->data[k];
          int index = seam_node_array[seam_index];
          if (index != min_index)
            int_int_unordered_map_insert(dup_node_map, index, min_index);
        }
      }
      int_array_free(same_nodes);
    }

    // Clean up.
    kd_tree_free(node_tree);
  }
  int* node_map = create_index_map_with_dups_removed(num_nodes, dup_node_map); 

  // Reduce the number of faces and nodes in the fused mesh by the ones that 
  // have been merged to others.
  num_faces -= dup_face_map->size;
  num_nodes -= dup_node_map->size;

  // We're through with the seam faces/nodes and the duplicate node maps.
  // We still need the duplicate face map below.
  int_unordered_set_free(seam_faces);
  int_unordered_set_free(seam_nodes);
  int_int_unordered_map_free(dup_node_map);

  // Now we create the fused mesh and fill it with the contents of the submeshes.
  log_debug("fuse_submeshes: Creating polymesh (%d cells, %d faces, %d nodes)",
            num_cells, num_faces, num_nodes);
  polymesh_t* fused_mesh = polymesh_new(comm, num_cells, num_ghost_cells,
                                        num_faces, num_nodes);

  // Allocate storage for cell faces and face nodes.
  fused_mesh->cell_face_offsets[0] = 0;
  fused_mesh->face_node_offsets[0] = 0;
  int cell = 0;
  for (int m = 0; m < num_submeshes; ++m)
  {
    polymesh_t* submesh = submeshes[m];
    for (int c = 0; c < submesh->num_cells; ++c, ++cell)
    {
      int num_cell_faces = submesh->cell_face_offsets[c+1] - submesh->cell_face_offsets[c];
      fused_mesh->cell_face_offsets[cell+1] = fused_mesh->cell_face_offsets[cell] + num_cell_faces;
    }
    for (int f = 0; f < submesh->num_faces; ++f)
    {
      int flattened_face = submesh_face_offsets[m] + f;
      int face = face_map[flattened_face];
      int num_face_nodes = submesh->face_node_offsets[f+1] - submesh->face_node_offsets[f];
      fused_mesh->face_node_offsets[face+1] = fused_mesh->face_node_offsets[face] + num_face_nodes;
    }
  }
  polymesh_reserve_connectivity_storage(fused_mesh);

  // Copy cell faces.
  cell = 0;
  for (int m = 0; m < num_submeshes; ++m)
  {
    polymesh_t* submesh = submeshes[m];
    for (int c = 0; c < submesh->num_cells; ++c, ++cell)
    {
      int num_cell_faces = submesh->cell_face_offsets[c+1] - submesh->cell_face_offsets[c];
      for (int f = 0; f < num_cell_faces; ++f)
      {
        int subface_index = submesh->cell_faces[submesh->cell_face_offsets[c] + f];
        bool flipped = false;
        if (subface_index < 0)
        {
          flipped = true;
          subface_index = ~subface_index;
        }
        int flattened_face = submesh_face_offsets[m] + subface_index;
        int face_index = face_map[flattened_face];
        ASSERT(face_index < fused_mesh->num_faces);
        if (flipped)
          face_index = ~face_index;
        fused_mesh->cell_faces[fused_mesh->cell_face_offsets[cell]+f] = face_index;
      }
    }
  }
  int_int_unordered_map_free(dup_face_map);

  // Hook up all the interior face cells.
  for (int icell = 0; icell < fused_mesh->num_cells; ++icell)
  {
    int fpos = 0, face;
    while (polymesh_cell_next_face(fused_mesh, icell, &fpos, &face))
    {
      if (fused_mesh->face_cells[2*face] == -1)
        fused_mesh->face_cells[2*face] = icell;
      else
      {
        ASSERT(fused_mesh->face_cells[2*face+1] == -1);
        fused_mesh->face_cells[2*face+1] = icell;
      }
    }
  }

  // Copy over the remaining face cells and do the face nodes.
  for (int m = 0; m < num_submeshes; ++m)
  {
    polymesh_t* submesh = submeshes[m];
    for (int f = 0; f < submesh->num_faces; ++f)
    {
      int flattened_face = submesh_face_offsets[m] + f;
      int face = face_map[flattened_face];

      // This face already has one cell attached.
      ASSERT(fused_mesh->face_cells[2*face] != -1);
      if (fused_mesh->face_cells[2*face+1] == -1)
      {
        // If the second cell is a ghost cell (< -1), it has its owning 
        // process encoded. If so, we simply copy it into place and 
        // deal with it later.
        if (submesh->face_cells[2*f+1] < -1)
          fused_mesh->face_cells[2*face+1] = submesh->face_cells[2*f+1];
      }

      // Set up the nodes of the face.
      int num_face_nodes = submesh->face_node_offsets[f+1] - submesh->face_node_offsets[f];
      for (int n = 0; n < num_face_nodes; ++n)
      {
        int subnode_index = submesh->face_nodes[submesh->face_node_offsets[f] + n];
        int flattened_node = submesh_node_offsets[m] + subnode_index;
        int node_index = node_map[flattened_node];
        ASSERT(node_index < fused_mesh->num_nodes);
        fused_mesh->face_nodes[fused_mesh->face_node_offsets[face]+n] = node_index;
      }
    }
  }
  polymec_free(face_map);

  // Copy node positions.
  for (int m = 0; m < num_submeshes; ++m)
  {
    polymesh_t* submesh = submeshes[m];
    for (int n = 0; n < submesh->num_nodes; ++n)
    {
      int flattened_node = submesh_node_offsets[m] + n;
      int node = node_map[flattened_node];
      fused_mesh->nodes[node] = submesh->nodes[n];
    }
  }
  polymec_free(node_map);

  // Construct edges.
  polymesh_construct_edges(fused_mesh);

  // Now compute geometry.
  polymesh_compute_geometry(fused_mesh);

  // Consume the submeshes.
  for (int m = 0; m < num_submeshes; ++m)
    polymesh_free(submeshes[m]);

  // Now fill the exchanger for the fused mesh with data.
  exchanger_proc_map_t* send_map = exchanger_proc_map_new(); 
  exchanger_proc_map_t* recv_map = exchanger_proc_map_new(); 
  int ghost_cell = fused_mesh->num_cells;
  for (int f = 0; f < fused_mesh->num_faces; ++f)
  {
    if (fused_mesh->face_cells[2*f+1] < -1)
    {
      // Found a ghost cell. Get its owning process.
      int proc = -fused_mesh->face_cells[2*f+1] - 2;
      ASSERT(proc != rank);
      ASSERT(proc >= 0);
      ASSERT(proc < nprocs);

      // Set up the sends and receives associated with this process.
      exchanger_proc_map_add_index(send_map, proc, fused_mesh->face_cells[2*f]);
      exchanger_proc_map_add_index(recv_map, proc, ghost_cell);

      // Make sure we attach the correct ghost cell to the face.
      fused_mesh->face_cells[2*f+1] = ghost_cell;
      ++ghost_cell;
    }
  }
  ASSERT(ghost_cell == (fused_mesh->num_cells + fused_mesh->num_ghost_cells));
  exchanger_t* fused_ex = polymesh_exchanger(fused_mesh, POLYMESH_CELL);
  exchanger_set_sends(fused_ex, send_map);
  exchanger_set_receives(fused_ex, recv_map);

  // Now that we've corrected all of our face->cell connections, we can verify the 
  // topological correctness of the fused mesh.
  ASSERT(polymesh_verify_topology(fused_mesh, polymec_error));

  // Return the final fused mesh.
  STOP_FUNCTION_TIMER();
  return fused_mesh;
}

#endif

bool partition_polymesh(polymesh_t** mesh, 
                        MPI_Comm comm, 
                        int* weights, 
                        real_t imbalance_tol,
                        polymesh_field_t** fields,
                        size_t num_fields)
{
#if POLYMEC_HAVE_MPI
  START_FUNCTION_TIMER();
  ASSERT((weights == NULL) || (imbalance_tol > 0.0));
  ASSERT((weights == NULL) || (imbalance_tol <= 1.0));
  ASSERT((*mesh == NULL) || ((*mesh)->comm == MPI_COMM_SELF));

  int nprocs, rank;
  MPI_Comm_size(comm, &nprocs);
  MPI_Comm_rank(comm, &rank);
  ASSERT((rank != 0) || (*mesh != NULL));

  // On a single process, partitioning has no meaning, but we do replace the communicator
  // if needed. NOTE: the migrator will still have its original communicator, but this 
  // shouldn't matter in any practical sense.
  if (nprocs == 1)
  {
    if (comm != (*mesh)->comm)
      (*mesh)->comm = comm;
    STOP_FUNCTION_TIMER();
    return true;
  }

  log_debug("partition_mesh: Partitioning mesh into %d subdomains.", nprocs);

  // If meshes on rank != 0 are not NULL, we delete them.
  polymesh_t* m = *mesh;
  if ((rank != 0) && (m != NULL))
  {
    polymesh_free(m);
    *mesh = m = NULL; 
  }

  // Generate a global adjacency graph for the mesh.
  adj_graph_t* global_graph = (m != NULL) ? graph_from_polymesh_cells(m) : NULL;

#ifndef NDEBUG
  // Make sure there are enough cells to go around for the processes we're given.
  if (rank == 0)
    ASSERT((*mesh)->num_cells > nprocs);
#endif

  // Map the graph to the different domains, producing a local partition vector.
  int64_t* global_partition = partition_graph(global_graph, comm, weights, imbalance_tol, true);
  if (global_graph != NULL)
    adj_graph_free(global_graph);
  if (global_partition == NULL)
  {
    STOP_FUNCTION_TIMER();
    return false;
  }

  // Distribute the mesh.
  log_debug("partition_mesh: Distributing mesh and %d fields to %d processes.", (int)num_fields, nprocs);
  distribute_polymesh(mesh, comm, global_partition, fields, num_fields);

  // Clean up.
  if (global_partition != NULL)
    polymec_free(global_partition);

  // Return the migrator.
  STOP_FUNCTION_TIMER();
  return true;
#else
  // Replace the communicator if needed.
  if (comm != (*mesh)->comm)
    (*mesh)->comm = comm;
  return true;
#endif
}

int64_t* partition_vector_from_polymesh(polymesh_t* global_mesh, 
                                        MPI_Comm comm, 
                                        int* weights, 
                                        real_t imbalance_tol,
                                        bool broadcast)
{
#if POLYMEC_HAVE_MPI
  START_FUNCTION_TIMER();
  ASSERT((weights == NULL) || (imbalance_tol > 0.0));
  ASSERT((weights == NULL) || (imbalance_tol <= 1.0));

  int nprocs, rank;
  MPI_Comm_size(comm, &nprocs);
  MPI_Comm_rank(comm, &rank);
  ASSERT((rank != 0) || (global_mesh != NULL));

  // On a single process, partitioning has no meaning.
  if (nprocs == 1)
  {
    // Dumb, but correct.
    int64_t* global_partition = polymec_calloc(global_mesh->num_cells, sizeof(int64_t));
    STOP_FUNCTION_TIMER();
    return global_partition;
  }

  // Generate a global adjacency graph for the mesh.
  adj_graph_t* global_graph = (rank == 0) ? graph_from_polymesh_cells(global_mesh) : NULL;

#ifndef NDEBUG
  // Make sure there are enough cells to go around for the processes we're given.
  if (rank == 0)
  {
    ASSERT(global_mesh->num_cells > nprocs);
  }
#endif

  // Map the graph to the different domains, producing a local partition vector.
  int64_t* global_partition = partition_graph(global_graph, comm, weights, imbalance_tol, broadcast);

  // Get rid of the graph.
  if (global_graph != NULL)
    adj_graph_free(global_graph);

  STOP_FUNCTION_TIMER();
  return global_partition;

#else
  // This is dumb, but we were asked for it.
  int64_t* global_partition = polymec_calloc(global_mesh->num_cells, sizeof(int64_t));
  return global_partition;
#endif
}

void distribute_polymesh(polymesh_t** mesh, 
                         MPI_Comm comm,
                         int64_t* global_partition,
                         polymesh_field_t** fields,
                         size_t num_fields)
{
#if POLYMEC_HAVE_MPI
  ASSERT((*mesh == NULL) || ((*mesh)->comm == MPI_COMM_SELF));

  int nprocs, rank;
  MPI_Comm_size(comm, &nprocs);
  MPI_Comm_rank(comm, &rank);
  ASSERT((rank != 0) || (global_partition != NULL));
  ASSERT((rank != 0) || (*mesh != NULL));

  // On a single process, distributing has no meaning.
  if (nprocs == 1)
    return;

  START_FUNCTION_TIMER();

  // Make sure we're all here.
  MPI_Barrier(comm);

  // We use the int_array serializer, so we need to make sure it's 
  // registered on all processes. See serializer.h for details.
  {
    serializer_t* s = int_array_serializer();
    (void)s;
  }

  polymesh_t* global_mesh = *mesh;
  polymesh_t* local_mesh = NULL;
  uint64_t vtx_dist[nprocs+1];
  int num_cells[nprocs];
  if (rank == 0)
  {
    // Take stock of how many cells we'll have per process.
    memset(num_cells, 0, sizeof(int) * nprocs);
    for (int i = 0; i < global_mesh->num_cells; ++i)
      num_cells[global_partition[i]]++;

    // Construct the distribution of vertices for the partitioning.
    vtx_dist[0] = 0;
    for (int p = 0; p < nprocs; ++p)
      vtx_dist[p+1] = vtx_dist[p] + num_cells[p];

    // Carve out the portion of the mesh that will stick around on process 0.
    {
      int indices[num_cells[0]], k = 0;
      for (int i = 0; i < global_mesh->num_cells; ++i)
      {
        if (global_partition[i] == rank)
          indices[k++] = i;
      }
      local_mesh = create_submesh(comm, global_mesh, global_partition, NULL, indices, num_cells[0]);
    }

    // Now do the other processes.
    serializer_t* ser = polymesh_serializer();
    byte_array_t* bytes = byte_array_new();
    for (int p = 1; p < nprocs; ++p)
    {
      // Share the vertex distribution.
      MPI_Send(vtx_dist, nprocs+1, MPI_UINT64_T, p, p, comm);

      // Create the pth submesh.
      int indices[num_cells[p]], k = 0;
      for (int i = 0; i < global_mesh->num_cells; ++i)
      {
        if (global_partition[i] == p)
          indices[k++] = i;
      }
      polymesh_t* p_mesh = create_submesh(comm, global_mesh, global_partition, NULL, indices, num_cells[p]);

      // Serialize it and send its size (and it) to process p.
      size_t offset = 0;
      serializer_write(ser, p_mesh, bytes, &offset);
      MPI_Send(&bytes->size, 1, MPI_INT, p, p, comm);
      MPI_Send(bytes->data, (int)bytes->size, MPI_BYTE, p, p, comm);

      // Clean up.
      byte_array_clear(bytes);
      polymesh_free(p_mesh);
    }
    byte_array_free(bytes);
  }
  else
  {
    // Receive the vertex distribution of the incoming mesh.
    MPI_Status status;
    MPI_Recv(vtx_dist, nprocs+1, MPI_UINT64_T, 0, rank, comm, &status);

    // Receive the size of the incoming mesh.
    int mesh_size;
    MPI_Recv(&mesh_size, 1, MPI_INT, 0, rank, comm, &status);

    // Now receive the mesh.
    byte_array_t* bytes = byte_array_new();
    byte_array_resize(bytes, mesh_size);

    MPI_Recv(bytes->data, mesh_size, MPI_BYTE, 0, rank, comm, &status);
    serializer_t* ser = polymesh_serializer();
    size_t offset = 0;
    local_mesh = serializer_read(ser, bytes, &offset);
    
    byte_array_free(bytes);
  }

  *mesh = local_mesh;

  // Clean up.
  if (global_mesh != NULL)
    polymesh_free(global_mesh);

  // Extract the boundary faces and the original (global) ghost cell indices
  // associated with them.
  ASSERT(polymesh_has_tag(local_mesh->face_tags, "parallel_boundary_faces"));
  ASSERT(polymesh_has_tag(local_mesh->cell_tags, "ghost_cell_indices"));
  size_t num_pbfaces, num_pbgcells;
  int* pbfaces = polymesh_tag(local_mesh->face_tags, "parallel_boundary_faces", &num_pbfaces);
  int* pbgcells = polymesh_tag(local_mesh->cell_tags, "ghost_cell_indices", &num_pbgcells);

  // Retrieve our global cell indices.
  ASSERT(polymesh_has_tag(local_mesh->cell_tags, "global_cell_indices"));
  size_t num_global_cell_indices;
  int* global_cell_indices = polymesh_tag(local_mesh->cell_tags, "global_cell_indices", &num_global_cell_indices);

  // Here's an inverse map for global cell indices to the local ones we have.
  int_int_unordered_map_t* inverse_cell_map = int_int_unordered_map_new();
  for (int i = 0; i < local_mesh->num_cells; ++i)
    int_int_unordered_map_insert(inverse_cell_map, global_cell_indices[i], i);

  // Now we create pairwise global cell indices for the exchanger.
  int_ptr_unordered_map_t* ghost_cell_indices = int_ptr_unordered_map_new();
  int next_ghost_index = local_mesh->num_cells;
  for (int i = 0; i < num_pbfaces; ++i)
  {
    int f = pbfaces[i];
    ASSERT(local_mesh->face_cells[2*f+1] < -1);
    // Get the destination process for this ghost cell.
    int proc = -local_mesh->face_cells[2*f+1] - 2;
    ASSERT(proc >= 0);
    ASSERT(proc < nprocs);

    // Generate a mapping from a global index to its ghost cell.
    int global_ghost_index = pbgcells[i];
    int local_ghost_index;
    int* ghost_index_p = int_int_unordered_map_get(inverse_cell_map, global_ghost_index);
    if (ghost_index_p == NULL)
    {
      local_ghost_index = next_ghost_index++;
      int_int_unordered_map_insert(inverse_cell_map, global_ghost_index, local_ghost_index);
    }
    else
      local_ghost_index = *ghost_index_p;
    local_mesh->face_cells[2*f+1] = local_ghost_index;

    if (!int_ptr_unordered_map_contains(ghost_cell_indices, proc))
      int_ptr_unordered_map_insert_with_v_dtor(ghost_cell_indices, proc, int_array_new(), DTOR(int_array_free));
    int_array_t* indices = *int_ptr_unordered_map_get(ghost_cell_indices, proc);
    int local_cell = local_mesh->face_cells[2*f];
    int global_cell = global_cell_indices[local_cell];
    int_array_append(indices, global_cell);
    int_array_append(indices, global_ghost_index);
  }

  // Make sure everything lines up.
  ASSERT(next_ghost_index - local_mesh->num_cells == local_mesh->num_ghost_cells);

  exchanger_t* ex = polymesh_exchanger(local_mesh, POLYMESH_CELL);
  int pos = 0, proc;
  int_array_t* indices;
  while (int_ptr_unordered_map_next(ghost_cell_indices, &pos, &proc, (void**)&indices))
  {
    if (proc != rank)
    {
      ASSERT((indices->size > 0) && ((indices->size % 2) == 0));
      int num_pairs = (int)indices->size/2;

      // Sort the indices array lexicographically by pairs so that all of the 
      // exchanger send/receive transactions have the same order across 
      // processes. This requires a specialized sort, since we have to 
      // arrange the integers within the pairs in ascending order, sort 
      // them, and then switch them back.
      sort_global_cell_pairs(indices->data, num_pairs);
      
      int send_indices[num_pairs], recv_indices[num_pairs];
      for (int i = 0; i < num_pairs; ++i)
      {
        send_indices[i] = *int_int_unordered_map_get(inverse_cell_map, indices->data[2*i]);
        ASSERT(send_indices[i] < local_mesh->num_cells);
        recv_indices[i] = *int_int_unordered_map_get(inverse_cell_map, indices->data[2*i+1]);
        ASSERT(recv_indices[i] >= local_mesh->num_cells);
      }
      exchanger_set_send(ex, proc, send_indices, num_pairs, true);
      exchanger_set_receive(ex, proc, recv_indices, num_pairs, true);
    }
  }

  // Clean up again.
  int_ptr_unordered_map_free(ghost_cell_indices);
  int_int_unordered_map_free(inverse_cell_map);

  // Remove the tags we used for ghost cells and global cell indices.
  polymesh_delete_tag(local_mesh->face_tags, "parallel_boundary_faces");
  polymesh_delete_tag(local_mesh->cell_tags, "ghost_cell_indices");
  polymesh_delete_tag(local_mesh->cell_tags, "global_cell_indices");

  // Now handle field data.
  polymesh_field_t* local_fields[num_fields];
  if (rank == 0)
  {
    // Do the local portion.
    for (size_t i = 0; i < num_fields; ++i)
    {
      polymesh_field_t* global_field = fields[i];
      ASSERT(global_field->centering == POLYMESH_CELL); // only cell-centered fields supported at the moment!
      polymesh_field_t* local_field = polymesh_field_new(local_mesh, POLYMESH_CELL, global_field->num_components);
      int num_comps = local_field->num_components;
      for (size_t j = 0; j < local_field->num_local_values; ++j)
        for (int c = 0; c < num_comps; ++c)
          local_field->data[num_comps*j+c] = global_field->data[num_comps*indices->data[j]+c];
      local_fields[i] = local_field;
    }
    // Distribute the remote portions.
    for (int p = 1; p < nprocs; ++p)
    {
      for (size_t i = 0; i < num_fields; ++i)
      {
        polymesh_field_t* global_field = fields[i];
        ASSERT(global_field->centering == POLYMESH_CELL); // only cell-centered fields supported at the moment!
        int num_comps = global_field->num_components;
        MPI_Send(&num_comps, 1, MPI_SIZE_T, p, p, comm);
        real_t local_vals[num_comps*num_cells[p]];
        for (size_t j = 0; j < num_cells[p]; ++j)
          for (int c = 0; c < num_comps; ++c)
            local_vals[num_comps*j+c] = global_field->data[num_comps*indices->data[j]+c];
        MPI_Send(local_vals, (int)(num_comps * num_cells[p]), MPI_REAL_T, p, p, comm);
      }
    }
  }
  else
  {
    MPI_Status status;
    for (size_t i = 0; i < num_fields; ++i)
    {
      int num_comps;
      MPI_Recv(&num_comps, 1, MPI_SIZE_T, 0, rank, comm, &status);
      polymesh_field_t* field = polymesh_field_new(local_mesh, POLYMESH_CELL, num_comps);
      size_t size = num_comps * field->num_local_values;
      MPI_Recv(field->data, (int)size, MPI_REAL_T, 0, rank, comm, &status);
      local_fields[i] = field;
    }
  }

  // Replace the global fields with the local fields.
  for (size_t i = 0; i < num_fields; ++i)
  {
    if (fields[i] != NULL)
      polymesh_field_free(fields[i]);
    fields[i] = local_fields[i];
  }

  STOP_FUNCTION_TIMER();
#endif
}

//------------------------------------------------------------------------
//                            Mesh repartitioning
//------------------------------------------------------------------------

#if POLYMEC_HAVE_MPI
static void redistribute_polymesh_with_graph(polymesh_t** mesh, 
                                             int64_t* local_partition, 
                                             adj_graph_t* local_graph, 
                                             polymesh_field_t** fields, 
                                             size_t num_fields)
{
#ifndef NDEBUG
  // Only cell-centered fields can be redistributed at the moment.
  for (size_t i = 0; i < num_fields; ++i)
  {
    ASSERT(fields[i]->centering == POLYMESH_CELL);
  }
#endif

  polymesh_t* m = *mesh;
  if (log_level() == LOG_DEBUG)
  {
    char P_string[m->num_cells * 16];
    sprintf(P_string, "P = [ ");
    for (int i = 0; i < m->num_cells; ++i)
    {
      char Pi[17];
      snprintf(Pi, 16, "%" PRIi64 " ", local_partition[i]);
      strcat(P_string, Pi);
    }
    strcat(P_string, "]");
    log_debug("redistribute_polymesh: %s", P_string);
  }

  START_FUNCTION_TIMER();
  index_t* vtx_dist = adj_graph_vertex_dist(local_graph);

  // Get redistribution data from the local partition vector.
  redistribution_t* redist = redistribution_from_partition((*mesh)->comm, 
                                                           local_partition, 
                                                           (*mesh)->num_cells);

  // Make a parallel-aware version of our partition vector.
  int64_t partition[m->num_cells + m->num_ghost_cells];
  memcpy(partition, local_partition, sizeof(int64_t) * m->num_cells);
  exchanger_t* mesh_ex = polymesh_exchanger(m, POLYMESH_CELL);
  exchanger_exchange(mesh_ex, partition, 1, 0, MPI_INT64_T);

  // Post receives for buffer sizes.
  size_t num_receives = redist->receive_procs->size;
  size_t num_sends = redist->send_procs->size;
  int receive_buffer_sizes[num_receives];
  MPI_Request requests[num_receives + num_sends];
  for (size_t p = 0; p < num_receives; ++p)
  {
    int proc = redist->receive_procs->data[p];
    MPI_Irecv(&receive_buffer_sizes[p], 1, MPI_INT, proc, 0, m->comm, &requests[p]);
  }

  // Build meshes to send to other processes.
  int_unordered_set_t* sent_cells = int_unordered_set_new();
  serializer_t* ser = polymesh_serializer();
  byte_array_t* send_buffers[num_sends];
  for (size_t p = 0; p < num_sends; ++p)
  {
    int proc = redist->send_procs->data[p];
    int_array_t* indices = redist->send_indices[p];
    byte_array_t* bytes = byte_array_new();

    // Add the indices of the cells we are sending.
    for (size_t i = 0; i < indices->size; ++i)
    {
      ASSERT(!int_unordered_set_contains(sent_cells, indices->data[i]));
      int_unordered_set_insert(sent_cells, indices->data[i]);
    }

    // Create the mesh to send. Recall that the submesh encodes the destination
    // process rank in the ghost cells referenced in submesh->face_cells.
    polymesh_t* submesh = create_submesh(m->comm, m, partition, vtx_dist, 
                                         indices->data, indices->size);
    log_debug("redistribute_polymesh: Redistributing %d cells to process %d.", 
              submesh->num_cells, proc);

    // Serialize the submesh.
    size_t offset = 0;
    serializer_write(ser, submesh, bytes, &offset);

    // Pack the relevant contents of the fields into the buffer.
    for (size_t i = 0; i < num_fields; ++i)
    {
      size_t num_comps = fields[i]->num_components;
      byte_array_write_ints(send_buffers[p], 1, &submesh->num_cells, &offset);
      real_t field_data[num_comps*submesh->num_cells];
      for (int j = 0; j < submesh->num_cells; ++j)
        for (size_t c = 0; c < num_comps; ++c)
          field_data[num_comps*j+c] = fields[i]->data[num_comps*indices->data[j]+c];
      byte_array_write_real_ts(send_buffers[p], num_comps*submesh->num_cells, field_data, &offset);
    }

    // Send the buffer size.
    MPI_Isend(&bytes->size, 1, MPI_INT, proc, 0, m->comm, &requests[p + num_receives]);

    // Clean up.
    polymesh_free(submesh);
    send_buffers[p] = bytes;
  }

  // Wait for the buffer sizes to be transmitted.
  MPI_Status statuses[num_receives + num_sends];
  MPI_Waitall((int)(num_receives + num_sends), requests, statuses);

  // Post receives for the actual messages.
  byte_array_t* receive_buffers[num_receives];
  for (size_t i = 0; i < num_receives; ++i)
  {
    receive_buffers[i] = byte_array_new();
    byte_array_resize(receive_buffers[i], receive_buffer_sizes[i]);
    MPI_Irecv(receive_buffers[i]->data, receive_buffer_sizes[i], MPI_BYTE, 
              redist->receive_procs->data[i], 0, m->comm, &requests[i]);
  }

  // Send the actual meshes and wait for receipt.
  for (size_t i = 0; i < num_sends; ++i)
  {
    MPI_Isend(send_buffers[i]->data, (int)send_buffers[i]->size, MPI_BYTE, 
              redist->send_procs->data[i], 0, m->comm, &requests[num_receives + i]);
  }
  MPI_Waitall((int)(num_receives + num_sends), requests, statuses);

  // Unpack the meshes.
  size_t receive_offsets[num_receives];
  polymesh_t* submeshes[1+num_receives];
  for (size_t i = 0; i < num_receives; ++i)
  {
    receive_offsets[i] = 0;
    submeshes[i+1] = serializer_read(ser, receive_buffers[i], &receive_offsets[i]);
  }

  // Clean up the send buffers and the serializer. We still need the 
  // receive buffer.
  for (size_t i = 0; i < num_sends; ++i)
    byte_array_free(send_buffers[i]);

  // Construct a local submesh and store it in submeshes[0]. This submesh
  // consists of all cells not sent to other processes.
  size_t num_cells = adj_graph_num_vertices(local_graph);
  size_t num_local_cells = num_cells - sent_cells->size;
  int local_cells[num_local_cells], j = 0;
  for (size_t i = 0; i < num_cells; ++i)
  {
    if (!int_unordered_set_contains(sent_cells, (int)i))
      local_cells[j++] = (int)i;
  }
  submeshes[0] = create_submesh(m->comm, m, partition, vtx_dist, 
                                local_cells, num_local_cells);

  // Fuse all the submeshes into a single mesh.
  *mesh = fuse_submeshes(submeshes, 1+num_receives);

  // Unpack the migrated field data.
  for (size_t i = 0; i < num_fields; ++i)
  {
    int num_comps = fields[i]->num_components;
    polymesh_field_t* new_field = polymesh_field_new(*mesh, POLYMESH_CELL, num_comps);

    // Local portion of the field.
    for (int k = 0; k < num_local_cells; ++k)
      for (size_t c = 0; c < num_comps; ++c)
        new_field->data[num_comps*k+c] = fields[i]->data[num_comps*local_cells[k]+c];
    size_t offset = num_comps*num_local_cells;

    // Remote portions.
    for (size_t r = 0; r < num_receives; ++r)
    {
      int num_values;
      byte_array_read_ints(receive_buffers[r], 1, &num_values, &receive_offsets[r]);
      byte_array_read_real_ts(receive_buffers[r], num_comps*num_values, &(new_field->data[offset]), &receive_offsets[r]);
      offset += num_comps*num_values;
    }

    // Out with the old, in with the new!
    polymesh_field_free(fields[i]);
    fields[i] = new_field;
  }

  // Clean up the rest.
  for (size_t i = 0; i < num_receives; ++i)
    byte_array_free(receive_buffers[i]);
  int_unordered_set_free(sent_cells);
  redistribution_free(redist);
  polymesh_free(m);

  STOP_FUNCTION_TIMER();
}
#endif

bool repartition_polymesh(polymesh_t** mesh, 
                          int* weights, 
                          real_t imbalance_tol,
                          polymesh_field_t** fields,
                          size_t num_fields)
{
#if POLYMEC_HAVE_MPI
  START_FUNCTION_TIMER();
  ASSERT(imbalance_tol > 0.0);
  ASSERT(imbalance_tol <= 1.0);
  polymesh_t* m = *mesh;

  int nprocs;
  MPI_Comm_size(m->comm, &nprocs);

  // On a single process, repartitioning has no meaning.
  if (nprocs == 1)
  {
    STOP_FUNCTION_TIMER();
    return true;
  }

  // Generate a local adjacency graph for the mesh.
  adj_graph_t* local_graph = graph_from_polymesh_cells(m);

#ifndef NDEBUG
  // Make sure there are enough cells to go around for the processes we're given.
  index_t* vtx_dist = adj_graph_vertex_dist(local_graph);
  index_t total_num_cells = vtx_dist[nprocs];
  ASSERT(total_num_cells >= nprocs);
#endif

  log_debug("repartition_mesh: Repartitioning polymesh (%d cells) and %d fields on %d domains...", 
            m->num_cells, (int)num_fields, nprocs);

  // Get the exchanger for the mesh.
  exchanger_t* mesh_ex = polymesh_exchanger(m, POLYMESH_CELL);

  // Map the graph to the different domains, producing a local partition vector.
  int64_t* local_partition = repartition_graph(local_graph, mesh_ex, m->num_ghost_cells, weights, imbalance_tol);
  if (local_partition == NULL)
    return false;

  // Redistribute the polymesh and its fields.
  redistribute_polymesh_with_graph(mesh, local_partition, local_graph, fields, num_fields);

  // Clean up.
  adj_graph_free(local_graph);
  polymec_free(local_partition);

  STOP_FUNCTION_TIMER();
  return true;
#else
  return true;
#endif
}

void redistribute_polymesh(polymesh_t** mesh, 
                           int64_t* local_partition, 
                           polymesh_field_t** fields, 
                           size_t num_fields)
{
#if POLYMEC_HAVE_MPI
  adj_graph_t* local_graph = graph_from_polymesh_cells(*mesh);
  redistribute_polymesh_with_graph(mesh, local_partition, local_graph, fields, num_fields);
  adj_graph_free(local_graph);
#endif
}

