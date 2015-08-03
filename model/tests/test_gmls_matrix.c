// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <string.h>
#include "cmockery.h"
#include "poisson_gmls_functional.h"
#include "model/gmls_matrix.h"

void test_gmls_matrix_ctor(void** state)
{
  point_cloud_t* points = NULL; // FIXME
  stencil_t* stencil = NULL; // FIXME
  point_weight_function_t* W = NULL; // FIXME
  gmls_matrix_t* matrix = stencil_based_gmls_matrix_new(W, points, stencil);

  // Clean up.
  gmls_matrix_free(matrix);
  point_cloud_free(points);
  stencil_free(stencil);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const UnitTest tests[] = 
  {
    unit_test(test_gmls_matrix_ctor)
  };
  return run_tests(tests);
}
