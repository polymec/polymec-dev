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

#include <gc/gc.h>
#include "core/array.h"
#include "core/point.h"

point_t* point_new(real_t x, real_t y, real_t z)
{
  point_t* p = GC_MALLOC(sizeof(point_t));
  p->x = x, p->y = y, p->z = z;
  return p;
}

void point_fprintf(point_t* x, FILE* stream)
{
  fprintf(stream, "(%g, %g, %g)\n", x->x, x->y, x->z);
}

vector_t* vector_new(real_t vx, real_t vy, real_t vz)
{
  vector_t* v = GC_MALLOC(sizeof(vector_t));
  v->x = vx, v->y = vy, v->z = vz;
  return v;
}

bbox_t* bbox_new(real_t x1, real_t x2, real_t y1, real_t y2, real_t z1, real_t z2)
{
  ASSERT(x1 < x2);
  ASSERT(y1 < y2);
  ASSERT(z1 < z2);
  bbox_t* b = GC_MALLOC(sizeof(bbox_t));
  b->x1 = x1;
  b->x2 = x2;
  b->y1 = y1;
  b->y2 = y2;
  b->z1 = z1;
  b->z2 = z2;
  return b;
}

void compute_orthonormal_basis(vector_t* e1, vector_t* e2, vector_t* e3)
{
  ASSERT(fabs(vector_mag(e1) - 1.0) < 1e-14);

  // Pick an arbitrary vector, e2, that is perpendicular to e1. One of these
  // should work.
  if (e1->x != 0.0)
  {
    e2->y = 1.0, e2->z = 1.0; 
    e2->x = -(e1->y + e1->z) / e1->x;
  }
  else if (e1->y != 0.0)
  {
    e2->x = 1.0, e2->z = 1.0; 
    e2->y = -(e1->x + e1->z) / e1->y;
  }
  else if (e1->z != 0.0)
  {
    e2->x = 1.0, e2->y = 1.0; 
    e2->z = -(e1->x + e1->y) / e1->z;
  }
  vector_normalize(e2);

  // e3 = e1 x e2.
  vector_cross(e1, e2, e3);
  ASSERT(vector_mag(e3) > 1e-14);
}

void bbox_grow(bbox_t* box, point_t* p)
{
  if (box->x1 > p->x)
    box->x1 = p->x;
  if (box->x2 < p->x)
    box->x2 = p->x;
  if (box->y1 > p->y)
    box->y1 = p->y;
  if (box->y2 < p->y)
    box->y2 = p->y;
  if (box->z1 > p->z)
    box->z1 = p->z;
  if (box->z2 < p->z)
    box->z2 = p->z;
}

void bbox_fprintf(bbox_t* box, FILE* stream)
{
  fprintf(stream, "bounding box [%g, %g] x [%g, %g] x [%g, %g]\n",
          box->x1, box->x2, box->y1, box->y2, box->z1, box->z2);
}

int* bbox_intersecting_processes(bbox_t* bbox, MPI_Comm comm, int* num_procs)
{
#if POLYMEC_HAVE_MPI
  int rank, nprocs;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &nprocs);

  if (nprocs == 1)
  {
    *num_procs = 0;
    return NULL;
  }

  // Gather all the data for the bounding boxes of all processes on this 
  // communicator.
  real_t my_data[6] = {bbox->x1, bbox->x2, bbox->y1, bbox->y2, bbox->z1, bbox->z2};
  real_t* all_data = polymec_malloc(sizeof(real_t) * 6 * nprocs);
  MPI_Allgather(my_data, 6, MPI_REAL_T, all_data, 6, MPI_REAL_T, comm);

  // Now just go over the boxes and find the ones that intersect.
  int_array_t* proc_array = int_array_new();
  for (int p = 0; p < nprocs; ++p)
  {
    bbox_t bbox_p = {.x1 = all_data[6*p  ], .x2 = all_data[6*p+1],
                     .y1 = all_data[6*p+2], .y2 = all_data[6*p+3],
                     .z1 = all_data[6*p+4], .z2 = all_data[6*p+5]};
    if ((p != rank) && bbox_intersects_bbox(bbox, &bbox_p))
      int_array_append(proc_array, p);
  }

  // Clean up.
  polymec_free(all_data);
  int* procs = proc_array->data;
  *num_procs = proc_array->size;
  int_array_release_data_and_free(proc_array);
  return procs;

#else
  *num_procs = 0;
  return NULL;
#endif
}

extern double orient3d(double* pa, double* pb, double* pc, double* pd);
bool points_are_coplanar(point_t* p1, point_t* p2, point_t* p3, point_t* p4)
{
#if POLYMEC_HAVE_DOUBLE_PRECISION
  return (orient3d((double*)p1, (double*)p2, (double*)p3, (double*)p4) == 0.0);
#else
  double da[3], db[3], dc[3], dd[3];
  da[0] = (double)p1->x; da[1] = (double)p1->y; da[2] = (double)p1->z;
  db[0] = (double)p2->x; db[1] = (double)p2->y; db[2] = (double)p2->z;
  dc[0] = (double)p3->x; dc[1] = (double)p3->y; dc[2] = (double)p3->z;
  dd[0] = (double)p4->x; dd[1] = (double)p4->y; dd[2] = (double)p4->z;
  return (orient3d(da, db, dc, dd) == 0.0);
#endif
}

bool all_points_are_coplanar(point_t* points, int num_points)
{
  ASSERT(num_points >= 3);
  if (num_points == 3)
    return true;
  else
  {
    for (int i = 3; i < num_points; ++i)
    {
      if (!points_are_coplanar(&points[0], &points[1], &points[2], &points[i]))
        return false;
    }
    return true;
  }
}

