// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/array.h"
#include "core/unordered_set.h"
#include "core/unordered_map.h"
#include "core/slist.h"
#include "geometry/create_rectilinear_polymesh.h"
#include "geometry/cubic_lattice.h"

polymesh_t* create_rectilinear_polymesh(MPI_Comm comm,
                                        real_t* xs, int nxs,
                                        real_t* ys, int nys,
                                        real_t* zs, int nzs)
{
  ASSERT(nxs > 1);
  ASSERT(nys > 1);
  ASSERT(nzs > 1);

#ifndef NDEBUG
  for (int i = 1; i < nxs; ++i)
    ASSERT(xs[i]-xs[i-1] > 0.0);
  for (int j = 1; j < nys; ++j)
    ASSERT(ys[j]-ys[j-1] > 0.0);
  for (int k = 1; k < nzs; ++k)
    ASSERT(zs[k]-zs[k-1] > 0.0);
#endif

  // Numbers of cells in each direction.
  int nx = nxs-1, ny = nys-1, nz = nzs-1;

  // Create a cubic lattice object for indexing.
  cubic_lattice_t* lattice = cubic_lattice_new(nx, ny, nz);

  // We start out with the naive partitioning in index space, and count
  // mesh entities.
  int nproc, rank;
  MPI_Comm_size(comm, &nproc);
  MPI_Comm_rank(comm, &rank);
  index_t total_num_cells = nx * ny * nz;
  index_t cells_per_proc = total_num_cells / nproc;
  index_int_unordered_map_t* local_faces = index_int_unordered_map_new();
  index_int_unordered_map_t* local_nodes = index_int_unordered_map_new();
  int num_cells = (int)((rank == nproc-1) ? total_num_cells - cells_per_proc*rank
                                          : cells_per_proc);
  int num_ghost_cells = 0, num_faces = 0, num_nodes = 0;
  for (int c = 0; c < num_cells; ++c)
  {
    // Figure out (i, j, k) indices for this cell.
    index_t global_cell_index = rank*cells_per_proc + c, i, j, k;
    cubic_lattice_get_cell_triple(lattice, global_cell_index, &i, &j, &k);

    // Count up faces and nodes and generate global-to-local mappings.
    index_t faces[6] = {cubic_lattice_xface(lattice, i, j, k),
                        cubic_lattice_xface(lattice, i+1, j, k),
                        cubic_lattice_yface(lattice, i, j, k),
                        cubic_lattice_yface(lattice, i, j+1, k),
                        cubic_lattice_zface(lattice, i, j, k),
                        cubic_lattice_zface(lattice, i, j, k+1)};
    for (int f = 0; f < 6; ++f)
    {
      if (!index_int_unordered_map_contains(local_faces, faces[f]))
        index_int_unordered_map_insert(local_faces, faces[f], num_faces++);
    }

    index_t nodes[8] = {cubic_lattice_node(lattice, i, j, k),
                        cubic_lattice_node(lattice, i, j, k+1),
                        cubic_lattice_node(lattice, i, j+1, k),
                        cubic_lattice_node(lattice, i, j+1, k+1),
                        cubic_lattice_node(lattice, i+1, j, k),
                        cubic_lattice_node(lattice, i+1, j, k+1),
                        cubic_lattice_node(lattice, i+1, j+1, k),
                        cubic_lattice_node(lattice, i+1, j+1, k+1)};
    for (int n = 0; n < 8; ++n)
    {
      if (!index_int_unordered_map_contains(local_nodes, nodes[n]))
        index_int_unordered_map_insert(local_nodes, nodes[n], num_nodes++);
    }

    // Ghost cells.
    index_t neighboring_cells[6];
    neighboring_cells[0] = (i == 0) ? -1 : cubic_lattice_cell(lattice, i-1, j, k);
    neighboring_cells[1] = (i == nx-1) ? -1 : cubic_lattice_cell(lattice, i+1, j, k);
    neighboring_cells[2] = (j == 0) ? -1 : cubic_lattice_cell(lattice, i, j-1, k);
    neighboring_cells[3] = (j == ny-1) ? -1 : cubic_lattice_cell(lattice, i, j+1, k);
    neighboring_cells[4] = (k == 0) ? -1 : cubic_lattice_cell(lattice, i, j, k-1);
    neighboring_cells[5] = (k == nz-1) ? -1 : cubic_lattice_cell(lattice, i, j, k+1);
    for (int ii = 0; ii < 6; ++ii)
    {
      if ((neighboring_cells[ii] != -1) &&
          ((neighboring_cells[ii] < cells_per_proc*rank) ||
          (neighboring_cells[ii] >= cells_per_proc*(rank+1))))
      {
        if (! ((rank == (nproc-1)) &&
               (neighboring_cells[ii] >= cells_per_proc*rank) &&
               (neighboring_cells[ii] < (cells_per_proc*rank + num_cells))))
          ++num_ghost_cells;
      }
    }
  }

  // Create the mesh.
  int num_faces_per_cell = 6, num_nodes_per_face = 4;
  polymesh_t* mesh = polymesh_new_with_cell_type(comm, num_cells, num_ghost_cells,
                                                 num_faces, num_nodes,
                                                 num_faces_per_cell,
                                                 num_nodes_per_face);

  int_unordered_set_t* processed_nodes = int_unordered_set_new();
  exchanger_proc_map_t* send_map = exchanger_proc_map_new();
  exchanger_proc_map_t* recv_map = exchanger_proc_map_new();

  int ghost_cell_index = num_cells;
  for (int cell = 0; cell < num_cells; ++cell)
  {
    // Figure out (i, j, k) indices for this cell.
    index_t global_cell_index = rank*cells_per_proc + cell, i, j, k;
    cubic_lattice_get_cell_triple(lattice, global_cell_index, &i, &j, &k);

    // Hook up the cell and faces.

    // Note that we define the nodes of the faces to be ordered such that
    // the right-hand rule produces a normal vector that always points
    // in the positive x, y, or z direction. Thus, the faces whose normal
    // vectors point in the opposite directions get their ones complement.
    if (comm == MPI_COMM_SELF)
    {
      // Use the dirt-simple indexing scheme!
      mesh->cell_faces[6*cell+0] = ~((int)cubic_lattice_xface(lattice, i, j, k));
      mesh->cell_faces[6*cell+1] =  (int)cubic_lattice_xface(lattice, i+1, j, k);
      mesh->cell_faces[6*cell+2] = ~((int)cubic_lattice_yface(lattice, i, j, k));
      mesh->cell_faces[6*cell+3] =  (int)cubic_lattice_yface(lattice, i, j+1, k);
      mesh->cell_faces[6*cell+4] = ~((int)cubic_lattice_zface(lattice, i, j, k));
      mesh->cell_faces[6*cell+5] =  (int)cubic_lattice_zface(lattice, i, j, k+1);
    }
    else
    {
      mesh->cell_faces[6*cell+0] = ~(*index_int_unordered_map_get(local_faces, cubic_lattice_xface(lattice, i, j, k)));
      mesh->cell_faces[6*cell+1] = *index_int_unordered_map_get(local_faces, cubic_lattice_xface(lattice, i+1, j, k));
      mesh->cell_faces[6*cell+2] = ~(*index_int_unordered_map_get(local_faces, cubic_lattice_yface(lattice, i, j, k)));
      mesh->cell_faces[6*cell+3] = *index_int_unordered_map_get(local_faces, cubic_lattice_yface(lattice, i, j+1, k));
      mesh->cell_faces[6*cell+4] = ~(*index_int_unordered_map_get(local_faces, cubic_lattice_zface(lattice, i, j, k)));
      mesh->cell_faces[6*cell+5] = *index_int_unordered_map_get(local_faces, cubic_lattice_zface(lattice, i, j, k+1));
    }

    // Hook up each face to its nodes.
    int nodes[6][4];

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
    if (comm == MPI_COMM_SELF)
    {
      // Use the dirt-simple indexing scheme!
      node_indices[0] = (int)cubic_lattice_node(lattice, i, j, k);
      node_indices[1] = (int)cubic_lattice_node(lattice, i+1, j, k);
      node_indices[2] = (int)cubic_lattice_node(lattice, i+1, j+1, k);
      node_indices[3] = (int)cubic_lattice_node(lattice, i, j+1, k);
      node_indices[4] = (int)cubic_lattice_node(lattice, i, j, k+1);
      node_indices[5] = (int)cubic_lattice_node(lattice, i+1, j, k+1);
      node_indices[6] = (int)cubic_lattice_node(lattice, i+1, j+1, k+1);
      node_indices[7] = (int)cubic_lattice_node(lattice, i, j+1, k+1);
    }
    else
    {
      node_indices[0] = *index_int_unordered_map_get(local_nodes, cubic_lattice_node(lattice, i, j, k));
      node_indices[1] = *index_int_unordered_map_get(local_nodes, cubic_lattice_node(lattice, i+1, j, k));
      node_indices[2] = *index_int_unordered_map_get(local_nodes, cubic_lattice_node(lattice, i+1, j+1, k));
      node_indices[3] = *index_int_unordered_map_get(local_nodes, cubic_lattice_node(lattice, i, j+1, k));
      node_indices[4] = *index_int_unordered_map_get(local_nodes, cubic_lattice_node(lattice, i, j, k+1));
      node_indices[5] = *index_int_unordered_map_get(local_nodes, cubic_lattice_node(lattice, i+1, j, k+1));
      node_indices[6] = *index_int_unordered_map_get(local_nodes, cubic_lattice_node(lattice, i+1, j+1, k+1));
      node_indices[7] = *index_int_unordered_map_get(local_nodes, cubic_lattice_node(lattice, i, j+1, k+1));
    }

    // Face 0 (-x) -- backward traversal
    nodes[0][0] = node_indices[7];
    nodes[0][1] = node_indices[4];
    nodes[0][2] = node_indices[0];
    nodes[0][3] = node_indices[3];

    // Face 1 (+x) -- forward traversal
    nodes[1][0] = node_indices[1];
    nodes[1][1] = node_indices[2];
    nodes[1][2] = node_indices[6];
    nodes[1][3] = node_indices[5];

    // Face 2 (-y) -- backward traversal
    nodes[2][0] = node_indices[4];
    nodes[2][1] = node_indices[5];
    nodes[2][2] = node_indices[1];
    nodes[2][3] = node_indices[0];

    // Face 3 (+y) -- forward traversal
    nodes[3][0] = node_indices[2];
    nodes[3][1] = node_indices[3];
    nodes[3][2] = node_indices[7];
    nodes[3][3] = node_indices[6];

    // Face 4 (-z) -- backward traversal
    nodes[4][0] = node_indices[0];
    nodes[4][1] = node_indices[1];
    nodes[4][2] = node_indices[2];
    nodes[4][3] = node_indices[3];

    // Face 5 (+z) -- forward traversal
    nodes[5][0] = node_indices[4];
    nodes[5][1] = node_indices[5];
    nodes[5][2] = node_indices[6];
    nodes[5][3] = node_indices[7];

    // Hook everything up.
    for (int f = 0; f < 6; ++f)
    {
      int face = mesh->cell_faces[6*cell+f];
      if (face < 0) face = ~face;
      if (mesh->face_cells[2*face] == -1)
        mesh->face_cells[2*face] = cell;
      else if (mesh->face_cells[2*face+1] == -1)
        mesh->face_cells[2*face+1] = cell;

      for (int n = 0; n < 4; ++n)
        mesh->face_nodes[4*face+n] = nodes[f][n];
    }

    // Assign the node positions for the cell.
    static const int i_offsets[] = {0, 1, 1, 0, 0, 1, 1, 0};
    static const int j_offsets[] = {0, 0, 1, 1, 0, 0, 1, 1};
    static const int k_offsets[] = {0, 0, 0, 0, 1, 1, 1, 1};
    for (int n = 0; n < 8; ++n)
    {
      if (!int_unordered_set_contains(processed_nodes, node_indices[n]))
      {
        point_t* node = &mesh->nodes[node_indices[n]];
        node->x = xs[i + i_offsets[n]];
        node->y = ys[j + j_offsets[n]];
        node->z = zs[k + k_offsets[n]];
        int_unordered_set_insert(processed_nodes, node_indices[n]);
      }
    }
  }

  // Hook up ghost cells and send/receive stuff. Note that we reverse the
  // order of the nested loops to make things consistent across processes.
  for (int ii = 0; ii < 6; ++ii)
  {
    for (int cell = 0; cell < num_cells; ++cell)
    {
      // Figure out (i, j, k) indices for this cell.
      index_t global_cell_index = rank*cells_per_proc + cell, i, j, k;
      cubic_lattice_get_cell_triple(lattice, global_cell_index, &i, &j, &k);
      index_t neighboring_cells[6];
      neighboring_cells[0] = (i == 0) ? -1 : cubic_lattice_cell(lattice, i-1, j, k);
      neighboring_cells[1] = (i == nx-1) ? -1 : cubic_lattice_cell(lattice, i+1, j, k);
      neighboring_cells[2] = (j == 0) ? -1 : cubic_lattice_cell(lattice, i, j-1, k);
      neighboring_cells[3] = (j == ny-1) ? -1 : cubic_lattice_cell(lattice, i, j+1, k);
      neighboring_cells[4] = (k == 0) ? -1 : cubic_lattice_cell(lattice, i, j, k-1);
      neighboring_cells[5] = (k == nz-1) ? -1 : cubic_lattice_cell(lattice, i, j, k+1);

      if ((neighboring_cells[ii] != -1) &&
          ((neighboring_cells[ii] < cells_per_proc*rank) ||
           (neighboring_cells[ii] >= cells_per_proc*(rank+1))))
      {
        int face = mesh->cell_faces[6*cell+ii];
        if (face < 0) face = ~face;
        ASSERT(mesh->face_cells[2*face] != -1);

        // First of all, let's make sure that the neighboring cell
        // actually belongs on a different process.
        if ((rank == (nproc-1)) &&
            (neighboring_cells[ii] >= cells_per_proc*rank) &&
            (neighboring_cells[ii] < (cells_per_proc*rank + num_cells)))
        {
//          mesh->face_cells[2*face+1] = -1;
          continue;
        }

        // Okay, it's a proper ghost cell.
        mesh->face_cells[2*face+1] = ghost_cell_index;

        // Generate the send mapping.
        int ghost_proc = (int)(MIN(neighboring_cells[ii] / cells_per_proc, nproc-1));
        exchanger_proc_map_add_index(send_map, ghost_proc, cell);

        // Generate the receive mapping.
        exchanger_proc_map_add_index(recv_map, ghost_proc, ghost_cell_index);
        ++ghost_cell_index;
      }
    }
  }
  ASSERT(ghost_cell_index == num_cells + num_ghost_cells);

  // Now we set up the exchanger.
  exchanger_t* ex = polymesh_exchanger(mesh, POLYMESH_CELL);
  exchanger_set_sends(ex, send_map);
  exchanger_set_receives(ex, recv_map);

  // Clean up.
  release_ref(lattice);
  index_int_unordered_map_free(local_faces);
  index_int_unordered_map_free(local_nodes);
  int_unordered_set_free(processed_nodes);

  // Construct edge information.
  polymesh_construct_edges(mesh);

  // Compute mesh geometry.
  polymesh_compute_geometry(mesh);

  return mesh;
}

polymesh_t* create_rectilinear_polymesh_on_rank(MPI_Comm comm,
                                                int rank,
                                                real_t* xs, int nxs,
                                                real_t* ys, int nys,
                                                real_t* zs, int nzs)
{
  ASSERT(comm != MPI_COMM_SELF);
  ASSERT(rank >= 0);

  int my_rank, nprocs;
  MPI_Comm_rank(comm, &my_rank);
  MPI_Comm_rank(comm, &nprocs);

  if (rank > nprocs)
    polymec_error("create_rectilinear_polymesh_on_rank: invalid rank: %d", rank);

  polymesh_t* mesh = NULL;
  if (my_rank == rank)
    mesh = create_rectilinear_polymesh(comm, xs, nxs, ys, nys, zs, nzs);
  return mesh;
}

void tag_rectilinear_polymesh_faces(polymesh_t* mesh,
                                    const char* x1_tag,
                                    const char* x2_tag,
                                    const char* y1_tag,
                                    const char* y2_tag,
                                    const char* z1_tag,
                                    const char* z2_tag)
{
  // Find the min/max extents of the mesh (in terms of face centers.
  bbox_t bbox = {.x1 = REAL_MAX, .x2 = -REAL_MAX,
                 .y1 = REAL_MAX, .y2 = -REAL_MAX,
                 .z1 = REAL_MAX, .z2 = -REAL_MAX};
  for (int f = 0; f < mesh->num_faces; ++f)
    bbox_grow(&bbox, &(mesh->face_centers[f]));

  // Figure out the boundary faces of the mesh by checking their face centers
  // against the boundary coordinates. We don't have to be too clever, since
  // the face centers should be basically right on top of the bounding box
  // coordinates (in a floating point sense).
  int_slist_t* x1_faces = int_slist_new();
  int_slist_t* x2_faces = int_slist_new();
  int_slist_t* y1_faces = int_slist_new();
  int_slist_t* y2_faces = int_slist_new();
  int_slist_t* z1_faces = int_slist_new();
  int_slist_t* z2_faces = int_slist_new();
  real_t face_tol = 1e-8;
  for (int f = 0; f < mesh->num_faces; ++f)
  {
    if (mesh->face_cells[2*f+1] == -1)
    {
      if (reals_nearly_equal(mesh->face_centers[f].x, bbox.x1, face_tol))
        int_slist_append(x1_faces, f);
      else if (reals_nearly_equal(mesh->face_centers[f].x, bbox.x2, face_tol))
        int_slist_append(x2_faces, f);
      if (reals_nearly_equal(mesh->face_centers[f].y, bbox.y1, face_tol))
        int_slist_append(y1_faces, f);
      else if (reals_nearly_equal(mesh->face_centers[f].y, bbox.y2, face_tol))
        int_slist_append(y2_faces, f);
      if (reals_nearly_equal(mesh->face_centers[f].z, bbox.z1, face_tol))
        int_slist_append(z1_faces, f);
      else if (reals_nearly_equal(mesh->face_centers[f].z, bbox.z2, face_tol))
        int_slist_append(z2_faces, f);
    }
  }

  // Now create the boundary tags and populate them.
  int_slist_node_t* iter = NULL;
  int f = 0, face;
  int* x1tag = polymesh_create_tag(mesh->face_tags, x1_tag, x1_faces->size);
  while (int_slist_next(x1_faces, &iter, &face))
    x1tag[f++] = face;

  iter = NULL;
  f = 0;
  int* x2tag = polymesh_create_tag(mesh->face_tags, x2_tag, x2_faces->size);
  while (int_slist_next(x2_faces, &iter, &face))
    x2tag[f++] = face;

  iter = NULL;
  f = 0;
  int* y1tag = polymesh_create_tag(mesh->face_tags, y1_tag, y1_faces->size);
  while (int_slist_next(y1_faces, &iter, &face))
    y1tag[f++] = face;

  iter = NULL;
  f = 0;
  int* y2tag = polymesh_create_tag(mesh->face_tags, y2_tag, y2_faces->size);
  while (int_slist_next(y2_faces, &iter, &face))
    y2tag[f++] = face;

  iter = NULL;
  f = 0;
  int* z1tag = polymesh_create_tag(mesh->face_tags, z1_tag, z1_faces->size);
  while (int_slist_next(z1_faces, &iter, &face))
    z1tag[f++] = face;

  iter = NULL;
  f = 0;
  int* z2tag = polymesh_create_tag(mesh->face_tags, z2_tag, z2_faces->size);
  while (int_slist_next(z2_faces, &iter, &face))
    z2tag[f++] = face;

  // Clean up.
  int_slist_free(x1_faces);
  int_slist_free(x2_faces);
  int_slist_free(y1_faces);
  int_slist_free(y2_faces);
  int_slist_free(z1_faces);
  int_slist_free(z2_faces);
}

