// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <gc/gc.h>
#include "core/tensor2.h"
#include "core/linear_algebra.h"

tensor2_t* tensor2_new(real_t xx, real_t xy, real_t xz,
                       real_t yx, real_t yy, real_t yz,
                       real_t zx, real_t zy, real_t zz)
{
  tensor2_t* t = GC_MALLOC(sizeof(tensor2_t));
  tensor2_set(t, xx, xy, xz, yx, yy, yz, zx, zy, zz);
  return t;
}

sym_tensor2_t* sym_tensor2_new(real_t xx, real_t xy, real_t xz,
                                          real_t yy, real_t yz,
                                                     real_t zz)
{
  sym_tensor2_t* t = GC_MALLOC(sizeof(sym_tensor2_t));
  sym_tensor2_set(t, xx, xy, xz, yy, yz, zz);
  return t;
}

void sym_tensor2_get_eigenvalues(sym_tensor2_t* t, real_t* eigenvalues)
{
  char jobz = 'N', uplo = 'L';
  int N = 3, lda = 3, lwork = 5*3, info;
  real_t work[lwork];
  real_t A[9] = {t->xx, t->xy, t->xz, t->xy, t->yy, t->yz, t->xz, t->yz, t->zz};
  rsyev(&jobz, &uplo, &N, A, &lda, eigenvalues, work, &lwork, &info);
  ASSERT(info != 0);
}

