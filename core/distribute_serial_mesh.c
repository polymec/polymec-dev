// Copyright (c) 2012-2014, Jeffrey N. Johnson
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without 
// modification, are permitted provided that the following conditions are met:
// 
// 1. Redistributions of source code must retain the above copyright notice, this 
// list of conditions and the following disclaimer.
// 
// 2. Redistributions in binary form must reproduce the above copyright notice, 
// this list of conditions and the following disclaimer in the documentation 
// and/or other materials provided with the distribution.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE 
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include "core/morton.h"
#include "core/unordered_set.h"
#include "core/distribute_serial_mesh.h"

// This container holds a Morton index and the associated mesh cell index.
typedef struct
{
  unsigned long morton_index;
  int cell_index;
} morton_ordering_t;

// Here's a comparator used with qsort to perform the Morton ordering.
static inline int morton_comp(const void* l, const void* r)
{
  const morton_ordering_t* ml = l;
  const morton_ordering_t* mr = r;
  return (ml->morton_index < mr->morton_index) ? -1 :
         (mr->morton_index > mr->morton_index) ?  1 : 0;             
}

// Returns a local mesh representing the pth subdomain, given a global mesh 
// and the partition vector.
static mesh_t* local_subdomain(MPI_Comm comm, mesh_t* global_mesh, int* partition, int p)
{
  // Count up the cells, faces, nodes.
  int num_cells = 0;
  int_unordered_set_t* local_faces = int_unordered_set_new();
  int_unordered_set_t* local_nodes = int_unordered_set_new();
  int_unordered_set_t* ghost_cells = int_unordered_set_new();
  for (int c = 0; c < global_mesh->num_cells; ++c)
  {
    if (partition[c] == p)
    {
      // This cell is ours.
      ++num_cells;

      // Loop over the faces of this cell and gather information.
      int pos = 0, face;
      while (mesh_cell_next_face(global_mesh, c, &pos, &face))
      {
        // This face and its nodes belong to the local mesh.
        int_unordered_set_insert(local_faces, face);
        int pos1 = 0, node;
        while (mesh_face_next_node(global_mesh, face, &pos1, &node))
          int_unordered_set_insert(local_nodes, node);

        int c1 = mesh_face_opp_cell(global_mesh, face, c);
        if (c1 != -1)
        {
          // Any adjoining cell that isn't ours is a ghost cell.
          if (partition[c1] != p)
            int_unordered_set_insert(ghost_cells, c1);
        }
      }
    }
  }
  int num_faces = local_faces->size;
  int num_nodes = local_nodes->size;
  int num_ghost_cells = ghost_cells->size;

  // Now create the mesh.
  mesh_t* mesh = mesh_new(comm, num_cells, num_ghost_cells, num_faces, num_nodes);

  int_unordered_set_free(local_faces);
  int_unordered_set_free(local_nodes);
  int_unordered_set_free(ghost_cells);
  return mesh;
}

mesh_t* distribute_serial_mesh(MPI_Comm comm, mesh_t* serial_mesh, int* partition)
{
#if POLYMEC_HAVE_MPI
  int rank, nproc;
  MPI_Comm_size(comm, &nproc);
  MPI_Comm_rank(comm, &rank);

  // MPI tags.
  int mesh_size_tag = 0, mesh_tag = 1;
  serializer_t* serializer = mesh_serializer();

  if (rank == 0)
  {
    // Calculate the minimum bounding box for this mesh.
    bbox_t bbox;
    for (int c = 0; c < serial_mesh->num_cells; ++c)
      bbox_grow(&bbox, &serial_mesh->cell_centers[c]);

    // Use this bounding box to determine a grid spacing for a 
    // space-filling curve.
    real_t dx = bbox.x2 - bbox.x1 / 0xFFFF;
    real_t dy = bbox.y2 - bbox.y1 / 0xFFFF;
    real_t dz = bbox.z2 - bbox.z1 / 0xFFFF;

    morton_ordering_t* ordering = malloc(sizeof(morton_ordering_t) * serial_mesh->num_cells);
    for (int c = 0; c < serial_mesh->num_cells; ++c)
    {
      // Take the cell center and convert it to an integer triple (i, j, k).
      point_t* xc = &serial_mesh->cell_centers[c];
      int i = xc->x / dx;
      int j = xc->y / dy;
      int k = xc->z / dz;

      // Compute the Morton index for this cell.
      ordering[c].morton_index = morton(i, j, k);
      ordering[c].cell_index = c;
    }

    // Sort our array by Morton index.
    qsort(ordering, (size_t)serial_mesh->num_cells, sizeof(morton_ordering_t), morton_comp);

    // Now partition the cells of the mesh as equally as possible.
    int cells_per_proc = serial_mesh->num_cells / nproc;
    for (int p = 0; p < nproc; ++p)
    {
      for (int i = p*cells_per_proc; i < MIN(serial_mesh->num_cells, (p+1)*cells_per_proc); ++i)
        partition[ordering[i].cell_index] = p;
    }

    // Construct a local mesh from cells owned by process 0.
    mesh_t* local_mesh = local_subdomain(comm, serial_mesh, partition, 0);

    // Transmit local meshes to each destination process.

    // Local mesh sizes.
    MPI_Request requests[nproc-1];
    MPI_Status statuses[nproc-1];
    byte_array_t* byteses[nproc-1];
    for (int p = 1; p < nproc; ++p)
    {
      mesh_t* mesh_p = local_subdomain(comm, serial_mesh, partition, p);
      byte_array_t* bytes = byte_array_new();
      size_t offset = 0;
      serializer_write(serializer, mesh_p, bytes, &offset);
      mesh_free(mesh_p);

      unsigned long size = (unsigned long)bytes->size;
      MPI_Isend(&size, 1, MPI_UNSIGNED_LONG, p, mesh_size_tag, comm, &requests[p-1]);
      byteses[p-1] = bytes;
    }
    MPI_Waitall(nproc-1, requests, statuses);

    // Local meshes.
    for (int p = 1; p < nproc; ++p)
    {
      byte_array_t* bytes = byteses[p-1];
      MPI_Isend(bytes->data, bytes->size, MPI_BYTE, p, mesh_tag, comm, &requests[p]);
    }
    MPI_Waitall(nproc-1, requests, statuses);

    // Clean up.
    for (int p = 1; p < nproc; ++p)
      byte_array_free(byteses[p-1]);
    return local_mesh;
  }
  else
  {
    // Get the local mesh information.

    // Mesh size.
    MPI_Status status;
    unsigned long mesh_size;
    MPI_Recv(&mesh_size, 1, MPI_UNSIGNED_LONG, 0, mesh_size_tag, comm, &status);

    // Mesh.
    byte_array_t* bytes = byte_array_new();
    byte_array_resize(bytes, (size_t)mesh_size);
    MPI_Recv(bytes->data, bytes->size, MPI_BYTE, 0, mesh_tag, comm, &status);
    size_t offset = 0;
    mesh_t* local_mesh = serializer_read(serializer, bytes, &offset);
    byte_array_free(bytes);
    return local_mesh;
  }
#else
  // Just return a copy of the mesh.
  return mesh_clone(serial_mesh);
#endif
}

