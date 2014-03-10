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

// This structure describes a chunk of mesh.
typedef struct
{
} mesh_chunk_t;

// This constructs a mesh chunk from a whole mesh and a partition vector.
static mesh_chunk_t* mesh_chunk_new(mesh_t* whole_mesh, int* partition)
{
  return NULL;
}

// This destroys a mesh chunk.
static void mesh_chunk_free(mesh_chunk_t* chunk)
{
}

// Returns the size of the given mesh chunk in bytes.
static int mesh_chunk_size(mesh_chunk_t* chunk)
{
  return 0;
}

// Converts the mesh chunk to a newly-allocated character buffer.
static char* mesh_chunk_to_bytes(mesh_chunk_t* chunk)
{
  return NULL;
}

// Reconstructs a mesh chunk from a character buffer.
static mesh_chunk_t* mesh_chunk_from_bytes(char* bytes, int size)
{
  return NULL;
}

// This constructs a local mesh on the given communicator from a chunk.
static mesh_t* mesh_from_chunk(MPI_Comm comm, mesh_chunk_t* chunk)
{
  return NULL;
}

mesh_t* distribute_serial_mesh(MPI_Comm comm, mesh_t* serial_mesh, int* partition)
{
  int rank, nproc;
  MPI_Comm_size(comm, &nproc);
  MPI_Comm_rank(comm, &rank);

#if POLYMEC_HAVE_MPI
  // MPI tags.
  int chunk_size_tag = 0, chunk_tag = 0;
#endif

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

    // Local mesh.
    mesh_chunk_t* chunk = mesh_chunk_new(serial_mesh, partition);
    mesh_t* mesh = mesh_from_chunk(comm, chunk);
    mesh_chunk_free(chunk);

#if POLYMEC_HAVE_MPI
    // Transmit information to each destination process.

    // Chunk sizes.
    MPI_Request requests[nproc-1];
    MPI_Status statuses[nproc-1];
    int sizes[nproc-1];
    char* byteses[nproc-1];
    for (int p = 1; p < nproc; ++p)
    {
      mesh_chunk_t* chunk = mesh_chunk_new(serial_mesh, partition);
      sizes[p-1] = mesh_chunk_size(chunk);
      byteses[p-1] = mesh_chunk_to_bytes(chunk);
      mesh_chunk_free(chunk);
      MPI_Isend(&sizes[p-1], 1, MPI_INT, p, chunk_size_tag, comm, &requests[p-1]);
    }
    MPI_Waitall(comm, requests, statuses);

    // Chunks.
    for (int p = 1; p < nproc; ++p)
    {
      char* bytes = byteses[p-1];
      int size = sizes[p-1];
      MPI_Isend(bytes, size, MPI_CHAR, p, chunk_tag, comm, &requests[p]);
    }
    MPI_Waitall(comm, requests, statuses);

    // Clean up.
    for (int p = 1; p < nproc; ++p)
      free(byteses[p-1]);
#endif
    return mesh;
  }
  else
  {
#if POLYMEC_HAVE_MPI
    // Get the chunk information.

    // Chunk size.
    MPI_Status status;
    int chunk_size;
    MPI_Recv(&chunk_size, 1, MPI_INT, 0, chunk_size_tag, comm, &status);

    // Chunk.
    char* bytes = malloc(sizeof(char) * chunk_size);
    MPI_Recv(&bytes, chunk_size, MPI_CHAR, 0, chunk_tag, comm, &status);
    mesh_chunk_t* chunk = mesh_chunk_from_bytes(buffer, chunk_size);
    free(bytes);

    // Construct the local mesh.
    mesh_t* mesh = mesh_from_chunk(comm, chunk);
    mesh_chunk_free(chunk);
    return mesh;
#else
    return NULL; // Never reached.
#endif
  }
}

