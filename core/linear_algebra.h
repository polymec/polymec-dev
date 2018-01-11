// Copyright (c) 2012-2018, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_LINEAR_ALGEBRA_H
#define POLYMEC_LINEAR_ALGEBRA_H

#include "core/polymec.h"

// Some LAPACK prototypes.

#ifdef LINUX
#define dgemv dgemv_
#define dgemm dgemm_
#define dgesv dgesv_
#define dgesvx dgesvx_
#define dgetrf dgetrf_
#define dgetrs dgetrs_
#define dposv dposv_
#define dposvx dposvx_
#define dpotrf dpotrf_
#define dpotrs dpotrs_
#define dgeqrf dgeqrf_
#define dorgqr dorgqr_
#define dormqr dormqr_
#define dgesvd dgesvd_
#define dgels dgels_
#define dgelsy dgelsy_
#define dgelss dgelss_
#define dgelsd dgelsd_
#define dlange dlange_
#define dgeev dgeev_
#define dsyev dsyev_
#define sgemv sgemv_
#define sgemm sgemm_
#define sgesv sgesv_
#define sgesvx sgesvx_
#define sgetrf sgetrf_
#define sgetrs sgetrs_
#define sposv sposv_
#define sposvx sposvx_
#define spotrf spotrf_
#define spotrs spotrs_
#define sgeqrf sgeqrf_
#define sorgqr sorgqr_
#define sormqr sormqr_
#define sgesvd sgesvd_
#define sgels  sgels_
#define sgelsy sgelsy_
#define sgelss sgelss_
#define sgelsd sgelsd_
#define slange slange_
#define sgeev sgeev_
#define ssyev ssyev_
#endif

// Matrix-vector multiplication: y := alpha*A*x + beta*y.
void dgemv(char* trans, int* m, int* n, double* alpha,
           double* a, int* lda, double* x, int* incx,
           double* beta, double* y, int* incy);
void sgemv(char* trans, int* m, int* n, float* alpha,
           float* a, int* lda, float* x, int* incx,
           float* beta, float* y, int* incy);
void rgemv(char* trans, int* m, int* n, real_t* alpha,
           real_t* a, int* lda, real_t* x, int* incx,
           real_t* beta, real_t* y, int* incy);

// Matrix-matrix multiplication: C := alpha*op(A)*op(B) + beta*C.
void dgemm(char* transa, char* transB, int* m, int* n, int* k, double* alpha,
           double* A, int* lda, double* B, int* ldb, double* beta, double* C, 
           int* ldc);
void sgemm(char* transa, char* transB, int* m, int* n, int* k, float* alpha,
           float* A, int* lda, float* B, int* ldb, float* beta, float* C, 
           int* ldc);
void rgemm(char* transa, char* transB, int* m, int* n, int* k, real_t* alpha,
           real_t* A, int* lda, real_t* B, int* ldb, real_t* beta, real_t* C, 
           int* ldc);

// Solves a linear system using LU factorization (in one step).
void dgesv(int* n, int* nrhs, double* A, int* lda, int* ipiv, 
           double* b, int* ldb, int* info);
void sgesv(int* n, int* nrhs, float* A, int* lda, int* ipiv, 
           float* b, int* ldb, int* info);
void rgesv(int* n, int* nrhs, real_t* A, int* lda, int* ipiv, 
           real_t* b, int* ldb, int* info);

// Solves a linear system using LU factorization (existing or not), providing 
// error bounds and an estimate of the condition number of the matrix A.
void dgesvx(char* fact, char* trans, int* n, int* nrhs, double* A,
            int* lda, double* af, int* ldaf, int* ipiv, char* equed,
            double* r, double* c, double* b, int* ldb, double* x, int* ldx,
            double* rcond, double* ferr, double* berr, double* work, 
            int* iwork, int* info);
void sgesvx(char* fact, char* trans, int* n, int* nrhs, float* A,
            int* lda, float* af, int* ldaf, int* ipiv, char* equed,
            float* r, float* c, float* b, int* ldb, float* x, int* ldx,
            float* rcond, float* ferr, float* berr, float* work, 
            int* iwork, int* info);
void rgesvx(char* fact, char* trans, int* n, int* nrhs, real_t* A,
            int* lda, real_t* af, int* ldaf, int* ipiv, char* equed,
            real_t* r, real_t* c, real_t* b, int* ldb, real_t* x, int* ldx,
            real_t* rcond, real_t* ferr, real_t* berr, real_t* work, 
            int* iwork, int* info);

// Computes an LU factorization.
void dgetrf(int* n, int* nrhs, double* A, int* lda, int* ipiv, int* info);
void sgetrf(int* n, int* nrhs, float* A, int* lda, int* ipiv, int* info);
void rgetrf(int* n, int* nrhs, real_t* A, int* lda, int* ipiv, int* info);

// Solves a linear system using an existing LU factorization.
void dgetrs(char* trans, int* n, int* nrhs, double* A, 
            int* lda, int* ipiv, double* b, int* ldb, int* info);
void sgetrs(char* trans, int* n, int* nrhs, float* A, 
            int* lda, int* ipiv, float* b, int* ldb, int* info);
void rgetrs(char* trans, int* n, int* nrhs, real_t* A, 
            int* lda, int* ipiv, real_t* b, int* ldb, int* info);

// Solves a linear system using Cholesky factorization (in one step).
// Note that A should be symmetric and positive definite.
void dposv(char* uplo, int* n, int* nrhs, double* A, int* lda, 
           double* b, int* ldb, int* info);
void sposv(char* uplo, int* n, int* nrhs, double* A, int* lda, 
           float* b, int* ldb, int* info);
void rposv(char* uplo, int* n, int* nrhs, double* A, int* lda, 
           real_t* b, int* ldb, int* info);

// Solves a linear system using Cholesky factorization (existing or not), 
// providing error bounds and an estimate of the condition number of the 
// (symmetric positive definite) matrix A.
void dposvx(char* fact, char* uplo, int* n, int* nrhs, double* A,
            int* lda, double* af, int* ldaf, char* equed, double* s, 
            double* b, int* ldb, double* x, int* ldx,
            double* rcond, double* ferr, double* berr, double* work, 
            int* iwork, int* info);
void sposvx(char* fact, char* uplo, int* n, int* nrhs, float* A,
            int* lda, float* af, int* ldaf, char* equed, float* s, 
            float* b, int* ldb, float* x, int* ldx,
            float* rcond, float* ferr, float* berr, float* work, 
            int* iwork, int* info);
void rposvx(char* fact, char* uplo, int* n, int* nrhs, real_t* A,
            int* lda, real_t* af, int* ldaf, char* equed, real_t* s, 
            real_t* b, int* ldb, real_t* x, int* ldx,
            real_t* rcond, real_t* ferr, real_t* berr, real_t* work, 
            int* iwork, int* info);

// Computes a Cholesky factorization.
// Note that A should be symmetric and positive definite.
void dpotrf(char* uplo, int* n, double* A, int* lda, int* info);
void spotrf(char* uplo, int* n, float* A, int* lda, int* info);
void rpotrf(char* uplo, int* n, real_t* A, int* lda, int* info);

// Solves a linear system using an existing Cholesky factorization.
// Note that A should be symmetric and positive definite.
void dpotrs(char* uplo, int* n, int* nrhs, double* A, 
            int* lda, double* b, int* ldb, int* info);
void spotrs(char* uplo, int* n, int* nrhs, float* A, 
            int* lda, float* b, int* ldb, int* info);
void rpotrs(char* uplo, int* n, int* nrhs, real_t* A, 
            int* lda, real_t* b, int* ldb, int* info);

// QR factorization.
// The matrix Q is represented as a product of elementary reflectors
//    Q = H(1) H(2) . . . H(k), where k = min(m,n).
//  Each H(i) has the form
//     H(i) = I - tau * v * vt
//  where tau is a real scalar, and v is a real vector with
//  v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i),
//  and tau in tau(i).
void dgeqrf(int* m, int* n, double* A, int* lda, double* tau, double* work,
            int* lwork, int* info);
void sgeqrf(int* m, int* n, float* A, int* lda, float* tau, float* work,
            int* lwork, int* info);
void rgeqrf(int* m, int* n, real_t* A, int* lda, real_t* tau, real_t* work,
            int* lwork, int* info);

// Generates an m x n real matrix Q with orthogonal columns, which is defined
// as the first N columns of a product of K elementary reflectors of order M
//    Q = H(1) H(2) . . . H(k)
//  as returned by DGEQRF.
void dorgqr(int* m, int* n, int* k, double* A, int* lda, double* tau,
            double* work, int* lwork, int* info);
void sorgqr(int* m, int* n, int* k, float* A, int* lda, float* tau,
            float* work, int* lwork, int* info);
void rorgqr(int* m, int* n, int* k, real_t* A, int* lda, real_t* tau,
            real_t* work, int* lwork, int* info);

// Overwrites the general real m x n matrix C with
//                  SIDE = 'L'     SIDE = 'R'
// TRANS = 'N':     Q * C          C * Q
// TRANS = 'T':     Qt * C         C * Qt
// where Q is a real orthogonal matrix defined as the product of k
// elementary reflectors
//       Q = H(1) H(2) . . . H(k)
// as returned by dgeqrf. Q is of order m if side = 'L' and of order n
// if side = 'R'.
void dormqr(char* side, char* trans, int* m, int* n, int* k, double* A,
            int* lda, double* tau, double* C, int* ldc, double* work, 
            int* lwork, int* info);
void sormqr(char* side, char* trans, int* m, int* n, int* k, float* A,
            int* lda, float* tau, float* C, int* ldc, float* work, 
            int* lwork, int* info);
void rormqr(char* side, char* trans, int* m, int* n, int* k, real_t* A,
            int* lda, real_t* tau, real_t* C, int* ldc, real_t* work, 
            int* lwork, int* info);

// Singular value decomposition: computes the SVD of an MxN matrix A:
// A = U * Sigma * transpose(V)
// where Sigma is a diagonal MxN matrix containing the singular values of A,
// U is an MxM orthogonal matrix, and V is an NxN orthogonal matrix. The columns 
// of U and V are the left and right singular vectors of A.
// Arguments:
// jobU - specifies options for computing the matrix U:
//        'A': all M columns of U are returned in array U:
//        'S':  the first min(m,n) columns of U (the left singular
//              vectors) are returned in the array U;
//        'O':  the first min(m,n) columns of U (the left singular
//              vectors) are overwritten on the array A;
//        'N':  no columns of U (no left singular vectors) are
//              computed.
// jobVT - specifies options for computing the matrix VT (V transpose):
//        'A': all M columns of VT are returned in array VT:
//        'S':  the first min(m,n) columns of VT (the right singular
//              vectors) are returned in the array VT;
//        'O':  the first min(m,n) columns of VT (the right singular
//              vectors) are overwritten on the array VT;
//        'N':  no columns of VT (no right singular vectors) are
//              computed.
// NOTE: jobU and jobVT cannot both be 'O'.
void dgesvd(char* jobU, char* jobVT, int* m, int* n, 
            double *A, int* lda, double* S, double* U, int* ldu, 
            double* VT, int* ldvt, double *work, int* lwork, int* info);
void sgesvd(char* jobU, char* jobVT, int* m, int* n, 
            float *A, int* lda, float* S, float* U, int* ldu, 
            float* VT, int* ldvt, float *work, int* lwork, int* info);
void rgesvd(char* jobU, char* jobVT, int* m, int* n, 
            real_t *A, int* lda, real_t* S, real_t* U, int* ldu, 
            real_t* VT, int* ldvt, real_t *work, int* lwork, int* info);

// DGELSY computes the minimum-norm solution to a real linear least
// squares problem:
//       minimize || A * X - B ||
// using a QR factorization of A. A is an M-by-N
// matrix which must not be rank-deficient. See LAPACK documentation for details.
void dgels(char* trans, int* m, int* n, int* nrhs, double* A, int* lda, double* B, int* ldb, 
           double* work, int* lwork, int* info);
void sgels(char* trans, int* m, int* n, int* nrhs, float* A, int* lda, float* B, int* ldb, 
           float* work, int* lwork, int* info);
void rgels(char* trans, int* m, int* n, int* nrhs, real_t* A, int* lda, real_t* B, int* ldb, 
           real_t* work, int* lwork, int* info);

// DGELSY computes the minimum-norm solution to a real linear least
// squares problem:
//       minimize || A * X - B ||
// using a complete orthogonal factorization of A. A is an M-by-N
// matrix which may be rank-deficient. See LAPACK documentation for details.
void dgelsy(int* m, int* n, int* nrhs, double* A, int* lda, double* B, int* ldb, 
            int* jpvt, double* rcond, int* rank, double* work, int* lwork, 
            int* info);
void sgelsy(int* m, int* n, int* nrhs, float* A, int* lda, float* B, int* ldb, 
            int* jpvt, float* rcond, int* rank, float* work, int* lwork, 
            int* info);
void rgelsy(int* m, int* n, int* nrhs, real_t* A, int* lda, real_t* B, int* ldb, 
            int* jpvt, real_t* rcond, int* rank, real_t* work, int* lwork, 
            int* info);

// DGELSS computes the minimum-norm solution to a real linear least
// squares problem:
//       minimize || A * X - B ||
// using the singular value decomposition of A. A is an M-by-N
// matrix which may be rank-deficient. See LAPACK documentation for details.
void dgelss(int* m, int* n, int* nrhs, double* A, int* lda, double* B, int* ldb, 
            double* S, double* rcond, int* rank, double* work, int* lwork, 
            int* info);
void sgelss(int* m, int* n, int* nrhs, float* A, int* lda, float* B, int* ldb, 
            float* S, float* rcond, int* rank, float* work, int* lwork, 
            int* info);
void rgelss(int* m, int* n, int* nrhs, real_t* A, int* lda, real_t* B, int* ldb, 
            real_t* S, real_t* rcond, int* rank, real_t* work, int* lwork, 
            int* info);

// DGELSD computes the minimum-norm solution to a real linear least
// squares problem:
//       minimize || A * X - B ||
// using the singular value decomposition of A with a divide-and-conquer 
// method. A is an M-by-N matrix which may be rank-deficient. See LAPACK 
// documentation for details.
void dgelsd(int* m, int* n, int* nrhs, double* A, int* lda, double* B, int* ldb, 
            double* S, double* rcond, int* rank, double* work, int* lwork, 
            int* info);
void sgelsd(int* m, int* n, int* nrhs, float* A, int* lda, float* B, int* ldb, 
            float* S, float* rcond, int* rank, float* work, int* lwork, 
            int* info);
void rgelsd(int* m, int* n, int* nrhs, real_t* A, int* lda, real_t* B, int* ldb, 
            real_t* S, real_t* rcond, int* rank, real_t* work, int* lwork, 
            int* info);

// Matrix norm: 'M' -> max column sum, 'I' -> Infinity, '1' -> 1-norm,
// 'F' -> Frobenius norm.
double dlange(char* norm, int* m, int* n, double* A, int* lda, double* work);
float slange(char* norm, int* m, int* n, float* A, int* lda, float* work);
real_t rlange(char* norm, int* m, int* n, real_t* A, int* lda, real_t* work);

// Eigenvalues and eigenvectors of a matrix A.
void dgeev(char* jobvl, char* jobvr, int* n, double* A, int* lda, double* wr, double* wi, 
           double* vl, int* ldvl, double* vr, int* ldvr, 
           double* work, int* lwork, int* info);
void sgeev(char* jobvl, char* jobvr, int* n, float* A, int* lda, float* wr, float* wi, 
           float* vl, int* ldvl, float* vr, int* ldvr, 
           float* work, int* lwork, int* info);
void rgeev(char* jobvl, char* jobvr, int* n, real_t* A, int* lda, real_t* wr, real_t* wi, 
           real_t* vl, int* ldvl, real_t* vr, int* ldvr, 
           real_t* work, int* lwork, int* info);

// Eigenvalues and eigenvectors of a real symmetric matrix A.
void dsyev(char* jobz, char* uplo, int* n, double* A, int* lda, double* W,
           double* work, int* lwork, int* info);
void ssyev(char* jobz, char* uplo, int* n, float* A, int* lda, float* W,
           float* work, int* lwork, int* info);
void rsyev(char* jobz, char* uplo, int* n, real_t* A, int* lda, real_t* W,
           real_t* work, int* lwork, int* info);

// Print a (column-major-ordered) matrix to the given file stream.
void matrix_fprintf(real_t* matrix, int nr, int nc, FILE* stream);

// Print a vector to the given file stream.
void vector_fprintf(real_t* vec, int nr, FILE* stream);

// Computes the determinant of the given 2x2 matrix.
real_t matrix2_det(real_t* matrix);

// Computes the determinant of the given 3x3 matrix.
real_t matrix3_det(real_t* matrix);

// Solves a 2x2 linear system Ax = b. 
// Note that b and x CAN point to the same vector.
void solve_2x2(real_t* A, real_t* b, real_t* x);

// Solves a 3x3 linear system Ax = b.
// Note that b and x CAN point to the same vector.
void solve_3x3(real_t* A, real_t* b, real_t* x);

#endif
