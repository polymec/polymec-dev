// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/linear_algebra.h"

// NOTE: All LAPACK functions are implemented within the library, not here.
// However, the real-valued prototypes are implemented here as 
// dispatches to their LAPACK functions.

void rgemv(char* trans, int* m, int* n, real_t* alpha,
           void* a, int* lda, void* x, int* incx,
           real_t* beta, void* y, int* incy)
{
#if POLYMEC_HAVE_DOUBLE_PRECISION
  dgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy);
#else
  sgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy);
#endif
}

void rgemm(char* transa, char* transB, int* m, int* n, int* k, real_t* alpha,
           real_t* A, int* lda, real_t* B, int* ldb, real_t* beta, real_t* C, 
           int* ldc)
{
#if POLYMEC_HAVE_DOUBLE_PRECISION
  dgemm(transa, transB, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
#else
  sgemm(transa, transB, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
#endif
}

void rgesv(int* n, int* nrhs, real_t* A, int* lda, int* ipiv, 
           real_t* b, int* ldb, int* info)
{
#if POLYMEC_HAVE_DOUBLE_PRECISION
  dgesv(n, nrhs, A, lda, ipiv, b, ldb, info);
#else
  sgesv(n, nrhs, A, lda, ipiv, b, ldb, info);
#endif
}

void rgesvx(char* fact, char* trans, int* n, int* nrhs, real_t* A,
            int* lda, real_t* af, int* ldaf, int* ipiv, char* equed,
            real_t* r, real_t* c, real_t* b, int* ldb, real_t* x, int* ldx,
            real_t* rcond, real_t* ferr, real_t* berr, real_t* work, 
            int* iwork, int* info)
{
#if POLYMEC_HAVE_DOUBLE_PRECISION
  dgesvx(fact, trans, n, nrhs, A, lda, af, ldaf, ipiv, equed, r, c, b, ldb, 
         x, ldx, rcond, ferr, berr, work, iwork, info);
#else
  sgesvx(fact, trans, n, nrhs, A, lda, af, ldaf, ipiv, equed, r, c, b, ldb, 
         x, ldx, rcond, ferr, berr, work, iwork, info);
#endif
}

void rgetrf(int* n, int* nrhs, real_t* A, int* lda, int* ipiv, int* info)
{
#if POLYMEC_HAVE_DOUBLE_PRECISION
  dgetrf(n, nrhs, A, lda, ipiv, info);
#else
  sgetrf(n, nrhs, A, lda, ipiv, info);
#endif
}

void rgetrs(char* trans, int* n, int* nrhs, real_t* A, 
            int* lda, int* ipiv, real_t* b, int* ldb, int* info)
{
#if POLYMEC_HAVE_DOUBLE_PRECISION
  dgetrs(trans, n, nrhs, A, lda, ipiv, b, ldb, info);
#else
  sgetrs(trans, n, nrhs, A, lda, ipiv, b, ldb, info);
#endif
}

void rposv(char* uplo, int* n, real_t* nrhs, double* A, int* lda, 
           real_t* b, int* ldb, int* info)
{
#if POLYMEC_HAVE_DOUBLE_PRECISION
  dposv(uplo, n, nrhs, A, lda, b, ldb, info);
#else
  sposv(uplo, n, nrhs, A, lda, b, ldb, info);
#endif
}

void rposvx(char* fact, char* uplo, int* n, int* nrhs, real_t* A,
            int* lda, real_t* af, int* ldaf, char* equed, real_t* s, 
            real_t* b, int* ldb, real_t* x, int* ldx,
            real_t* rcond, real_t* ferr, real_t* berr, real_t* work, 
            int* iwork, int* info)
{
#if POLYMEC_HAVE_DOUBLE_PRECISION
  dposvx(fact, uplo, n, nrhs, A, lda, af, ldaf, equed, s, b, ldb, 
         x, ldx, rcond, ferr, berr, work, iwork, info);
#else
  sposvx(fact, uplo, n, nrhs, A, lda, af, ldaf, equed, s, b, ldb, 
         x, ldx, rcond, ferr, berr, work, iwork, info);
#endif
}

void rpotrf(char* uplo, int* n, real_t* A, int* lda, int* info)
{
#if POLYMEC_HAVE_DOUBLE_PRECISION
  dpotrf(uplo, n, A, lda, info);
#else
  spotrf(uplo, n, A, lda, info);
#endif
}

void rpotrs(char* uplo, int* n, int* nrhs, real_t* A, 
            int* lda, real_t* b, int* ldb, int* info)
{
#if POLYMEC_HAVE_DOUBLE_PRECISION
  dpotrs(uplo, n, nrhs, A, lda, b, ldb, info);
#else
  spotrs(uplo, n, nrhs, A, lda, b, ldb, info);
#endif
}

void rgeqrf(int* m, int* n, real_t* A, int* lda, real_t* tau, real_t* work,
            int* lwork, int* info)
{
#if POLYMEC_HAVE_DOUBLE_PRECISION
  dgeqrf(m, n, A, lda, tau, work, lwork, info);
#else
  sgeqrf(m, n, A, lda, tau, work, lwork, info);
#endif
}

void rorgqr(int* m, int* n, int* k, real_t* A, int* lda, real_t* tau,
            real_t* work, int* lwork, int* info)
{
#if POLYMEC_HAVE_DOUBLE_PRECISION
  dorgqr(m, n, k, A, lda, tau, work, lwork, info);
#else
  sorgqr(m, n, k, A, lda, tau, work, lwork, info);
#endif
}

void rormqr(char* side, char* trans, int* m, int* n, int* k, real_t* A,
            int* lda, real_t* tau, real_t* C, int* ldc, real_t* work, 
            int* lwork, int* info)
{
#if POLYMEC_HAVE_DOUBLE_PRECISION
  dormqr(side, trans, m, n, k, A, lda, tau, C, ldc, work, lwork, info);
#else
  sormqr(side, trans, m, n, k, A, lda, tau, C, ldc, work, lwork, info);
#endif
}

void rgesvd(char* jobU, char* jobVT, int* m, int* n, 
           real_t *A, int* lda, real_t* S, real_t* U, int* ldu, 
           real_t* VT, int* ldvt, real_t *work, int* lwork, int* info)
{
#if POLYMEC_HAVE_DOUBLE_PRECISION
  dgesvd(jobU, jobVT, m, n, A, lda, S, U, ldu, VT, ldvt, work, lwork, info);
#else
  sgesvd(jobU, jobVT, m, n, A, lda, S, U, ldu, VT, ldvt, work, lwork, info);
#endif
}

void rgels(char* trans, int* m, int* n, int* nrhs, real_t* A, int* lda, real_t* B, int* ldb, 
           real_t* work, int* lwork, int* info)
{
#if POLYMEC_HAVE_DOUBLE_PRECISION
  dgels(trans, m, n, nrhs, A, lda, B, ldb, work, lwork, info);
#else
  sgels(trans, m, n, nrhs, A, lda, B, ldb, work, lwork, info);
#endif
}

void rgelsy(int* m, int* n, int* nrhs, real_t* A, int* lda, real_t* B, int* ldb, 
            int* jpvt, real_t* rcond, int* rank, real_t* work, int* lwork, 
            int* info)
{
#if POLYMEC_HAVE_DOUBLE_PRECISION
  dgelsy(m, n, nrhs, A, lda, B, ldb, jpvt, rcond, rank, work, lwork, info);
#else
  sgelsy(m, n, nrhs, A, lda, B, ldb, jpvt, rcond, rank, work, lwork, info);
#endif
}

void rgelss(int* m, int* n, int* nrhs, real_t* A, int* lda, real_t* B, int* ldb, 
            real_t* S, real_t* rcond, int* rank, real_t* work, int* lwork, 
            int* info)
{
#if POLYMEC_HAVE_DOUBLE_PRECISION
  dgelss(m, n, nrhs, A, lda, B, ldb, S, rcond, rank, work, lwork, info);
#else
  sgelss(m, n, nrhs, A, lda, B, ldb, S, rcond, rank, work, lwork, info);
#endif
}

void rgelsd(int* m, int* n, int* nrhs, real_t* A, int* lda, real_t* B, int* ldb, 
            real_t* S, real_t* rcond, int* rank, real_t* work, int* lwork, 
            int* info)
{
#if POLYMEC_HAVE_DOUBLE_PRECISION
  dgelsd(m, n, nrhs, A, lda, B, ldb, S, rcond, rank, work, lwork, info);
#else
  sgelsd(m, n, nrhs, A, lda, B, ldb, S, rcond, rank, work, lwork, info);
#endif
}

void matrix_fprintf(real_t* matrix, int nr, int nc, FILE* stream)
{
  if (stream == NULL) return;
  fprintf(stream, "[ ");
  for (int i = 0; i < nr; ++i)
  {
    if (i > 0)
      fprintf(stream, "  ");
    for (int j = 0; j < nc; ++j)
      fprintf(stream, "%g ", matrix[nr*j+i]);
    if (i < (nr -1))
      fprintf(stream, ";\n");
  }
  fprintf(stream, "]");
}

void vector_fprintf(real_t* vec, int nr, FILE* stream)
{
  if (stream == NULL) return;
  fprintf(stream, "[");
  for (int i = 0; i < nr; ++i)
    fprintf(stream, "%g ", vec[i]);
  fprintf(stream, "]");
}

real_t matrix2_det(real_t* matrix)
{
  return matrix[0]*matrix[3] - matrix[1]*matrix[2];
}

real_t matrix3_det(real_t* matrix)
{
  return matrix[0]*(matrix[4]*matrix[8] - matrix[5]*matrix[7]) -
         matrix[3]*(matrix[1]*matrix[8] - matrix[2]*matrix[7]) + 
         matrix[6]*(matrix[1]*matrix[5] - matrix[2]*matrix[4]);
}

void solve_2x2(real_t* A, real_t* b, real_t* x)
{
  real_t b0 = b[0], b1 = b[1];
  real_t inv_det_A = 1.0 / matrix2_det(A);
  x[0] = inv_det_A * ( A[3]*b0 - A[2]*b1);
  x[1] = inv_det_A * (-A[1]*b0 + A[0]*b1);
}

void solve_3x3(real_t* A, real_t* b, real_t* x)
{
  real_t b0 = b[0], b1 = b[1], b2 = b[2];

  // x = Ainv * b.
  real_t inv_det_A = 1.0 / matrix3_det(A);

  x[0] = inv_det_A * 
         ((A[8]*A[4]-A[5]*A[7]) * b0 - 
          (A[8]*A[3]-A[5]*A[6]) * b1 + 
          (A[7]*A[3]-A[4]*A[6]) * b2);

  x[1] = inv_det_A * 
         (-(A[8]*A[1]-A[2]*A[7]) * b0 +
           (A[8]*A[0]-A[2]*A[6]) * b1 -
           (A[7]*A[0]-A[1]*A[6]) * b2);

  x[2] = inv_det_A * 
         ((A[5]*A[1]-A[2]*A[4]) * b0 -
          (A[5]*A[0]-A[2]*A[3]) * b1 +
          (A[4]*A[0]-A[1]*A[3]) * b2);
}

