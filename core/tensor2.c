// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/tensor2.h"
#include "core/linear_algebra.h"

tensor2_t* tensor2_new(real_t xx, real_t xy, real_t xz,
                       real_t yx, real_t yy, real_t yz,
                       real_t zx, real_t zy, real_t zz)
{
  tensor2_t* t = polymec_refcounted_malloc(sizeof(tensor2_t), NULL);
  tensor2_set(t, xx, xy, xz, yx, yy, yz, zx, zy, zz);
  return t;
}

void tensor2_fprintf(tensor2_t* t, FILE* stream)
{
  fprintf(stream, "[%g, %g, %g]\n", t->xx, t->xy, t->xz);
  fprintf(stream, "[%g, %g, %g]\n", t->yx, t->yy, t->yz);
  fprintf(stream, "[%g, %g, %g]\n", t->zx, t->zy, t->zz);
}

symtensor2_t* symtensor2_new(real_t xx, real_t xy, real_t xz,
                                        real_t yy, real_t yz,
                                                   real_t zz)
{
  symtensor2_t* t = polymec_refcounted_malloc(sizeof(symtensor2_t), NULL);
  symtensor2_set(t, xx, xy, xz, yy, yz, zz);
  return t;
}

void symtensor2_get_eigenvalues(symtensor2_t* t, real_t eigenvalues[3])
{
  char jobz = 'N', uplo = 'L';
  int N = 3, lda = 3, lwork = 5*3, info;
  real_t work[lwork];
  real_t A[9] = {t->xx, t->xy, t->xz, t->xy, t->yy, t->yz, t->xz, t->yz, t->zz};
  rsyev(&jobz, &uplo, &N, A, &lda, eigenvalues, work, &lwork, &info);
  ASSERT(info == 0);
}

void symtensor2_get_eigenvectors(symtensor2_t* t, real_t eigenvalues[3], vector_t eigenvectors[3])
{
  char jobz = 'V', uplo = 'L';
  int N = 3, lda = 3, lwork = 5*3, info;
  real_t work[lwork];
  real_t A[9] = {t->xx, t->xy, t->xz, t->xy, t->yy, t->yz, t->xz, t->yz, t->zz};
  rsyev(&jobz, &uplo, &N, A, &lda, eigenvalues, work, &lwork, &info);
  ASSERT(info == 0);
  eigenvectors[0].x = A[0];
  eigenvectors[0].y = A[1];
  eigenvectors[0].z = A[2];
  eigenvectors[1].x = A[3];
  eigenvectors[1].y = A[4];
  eigenvectors[1].z = A[5];
  eigenvectors[2].x = A[6];
  eigenvectors[2].y = A[7];
  eigenvectors[2].z = A[8];
}

void symtensor2_fprintf(symtensor2_t* t, FILE* stream)
{
  fprintf(stream, "[%g, %g, %g]\n", t->xx, t->xy, t->xz);
  fprintf(stream, "[%g, %g, %g]\n", t->xy, t->yy, t->yz);
  fprintf(stream, "[%g, %g, %g]\n", t->xz, t->yz, t->zz);
}

