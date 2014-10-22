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

#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <string.h>
#include "cmockery.h"
#include "core/kd_tree.h"
#include "geometry/create_point_lattice.h"
#include "model/point_weight_function.h"
#include "model/neighbor_pairing.h"

// These functions implement a hat-shaped weight function with constant 
// spatial extent h.
static real_t hat_W(void* context, vector_t* y)
{
  real_t h = *((real_t*)context); // Spatial extent of W.
  real_t D = vector_mag(y);
  return (D < h) ? 1.0 - D/h : 0.0;
}

static vector_t hat_grad(void* context, vector_t* y)
{
  real_t h = *((real_t*)context); // Spatial extent of W.
  real_t D = vector_mag(y);
  static vector_t zero = {.x = 0.0, .y = 0.0, .z = 0.0};
  vector_t g = {.x = -y->x/(D*h), .y = -y->y/(D*h), .z = -y->z/(D*h)};
  return (D < h) ? g : zero;
}

static void hat_dtor(void* context)
{
  polymec_free(context);
}

static point_weight_function_t* hat_function_new(real_t h)
{
  real_t* Wh = polymec_malloc(sizeof(real_t));
  *Wh = h;
  point_weight_function_vtable vtable = {.value = hat_W, 
                                         .gradient = hat_grad, 
                                         .dtor = hat_dtor};
  return point_weight_function_new("Hat function", Wh, vtable);
}

static neighbor_pairing_t* create_simple_pairing(point_cloud_t* cloud,
                                                 real_t h)
{
  point_weight_function_t* W = hat_function_new(h);

  // Toss all the points into a kd-tree.
  kd_tree_t* tree = kd_tree_new(cloud->points, cloud->num_points);

  // Now do a neighbor search -- everything that falls within h of a point 
  // is a neighbor.
  int_array_t* pairs = int_array_new();
  real_array_t* weights = real_array_new();
  for (int i = 0; i < cloud->num_points; ++i)
  {
    point_t* xi = &cloud->points[i];
    int_slist_t* neighbors = kd_tree_within_radius(tree, &cloud->points[i], h);
    int_slist_node_t* n = NULL;
    int j;
    while (int_slist_next(neighbors, &n, &j))
    {
      int_array_append(pairs, i);
      int_array_append(pairs, j);
      point_t* xj = &cloud->points[j];
      vector_t y;
      point_displacement(xj, xi, &y);
      real_array_append(weights, point_weight_function_value(W, &y));
    }
    int_slist_free(neighbors);
  }
  kd_tree_free(tree);
  point_weight_function_free(W);

  exchanger_t* ex = exchanger_new(MPI_COMM_SELF);
  neighbor_pairing_t* pairing = neighbor_pairing_new("simple pairing", 
                                                     pairs->size/2,
                                                     pairs->data,
                                                     weights->data,
                                                     ex);
  int_array_release_data_and_free(pairs);
  real_array_release_data_and_free(weights);
  return pairing;
}

void test_point_lattice(void** state, 
                        MPI_Comm comm,
                        int nx, int ny, int nz, real_t h,
                        int num_interior_neighbors, 
                        int num_boundary_neighbors, 
                        int num_edge_neighbors,
                        int num_corner_neighbors)
{
  bbox_t bbox = {.x1 = 0.0, .x2 = 1.0, .y1 = 0.0, .y2 = 1.0, .z1 = 0.0, .z2 = 1.0};
  point_cloud_t* cloud = create_uniform_point_lattice(comm, nx, ny, nz, &bbox);
  neighbor_pairing_t* pairing = create_simple_pairing(cloud, h);
  // FIXME
  neighbor_pairing_free(pairing);
  point_cloud_free(cloud);
}

void test_serial_1x1x1_lattice(void** state)
{
  test_point_lattice(state, MPI_COMM_SELF, 1, 1, 1, 0.1, 0, 0, 0, 0);
}

void test_serial_10x1x1_lattice(void** state)
{
  test_point_lattice(state, MPI_COMM_SELF, 10, 1, 1, 0.15, 0, 0, 2, 1);
}

void test_serial_10x10x1_lattice(void** state)
{
  test_point_lattice(state, MPI_COMM_SELF, 10, 10, 1, 0.15, 4, 4, 3, 2);
}

void test_serial_10x10x10_lattice(void** state)
{
  test_point_lattice(state, MPI_COMM_SELF, 10, 10, 10, 0.15, 6, 5, 4, 3);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const UnitTest tests[] = 
  {
    unit_test(test_serial_1x1x1_lattice),
    unit_test(test_serial_10x1x1_lattice),
    unit_test(test_serial_10x10x1_lattice),
    unit_test(test_serial_10x10x10_lattice)
  };
  return run_tests(tests);
}
