// Copyright (c) 2012-2013, Jeffrey N. Johnson
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

#ifndef POLYMEC_LINEAR_ALGEBRA_H
#define POLYMEC_LINEAR_ALGEBRA_H

#include "core/polymec.h"

// Some LAPACK prototypes.

#ifdef LINUX
#define dgemv dgemv_
#define dgemm dgemm_
#define dgesv dgesv_
#define dgetrf dgetrf_
#define dgetrs dgetrs_
#define dgeqrf dgeqrf_
#define dorgqr dorgqr_
#define dormqr dormqr_
#define dgesvd dgesvd_
#define dgelsy dgelsy_
#endif

// Matrix-vector multiplication: y := alpha*A*x + beta*y.
void dgemv(char* trans, int* m, int* n, double* alpha,
           void* a, int* lda, void* x, int* incx,
           double* beta, void* y, int* incy);

// Matrix-matrix multiplication: C := alpha*op(A)*op(B) + beta*C.
void dgemm(char* transa, char* transB, int* m, int* n, int* k, double* alpha,
           double* A, int* lda, double* B, int* ldb, double* beta, double* C, 
           int* ldc);

// Solves a linear system using LU factorization (in one step).
void dgesv(int* n, int* nrhs, double* A, int* lda, int* ipiv, 
           double* b, int* ldb, int* info);

// Computes an LU factorization.
void dgetrf(int* n, int* nrhs, double* A, int* lda, int* ipiv, int* info);

// Solves a linear system using an existing LU factorization.
void dgetrs(char* trans, int* n, int* nrhs, double* A, 
            int* lda, int* ipiv, double* b, int* ldb, int* info);

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

// Generates an m x n real matrix Q with orthogonal columns, which is defined
// as the first N columns of a product of K elementary reflectors of order M
//    Q = H(1) H(2) . . . H(k)
//  as returned by DGEQRF.
void dorgqr(int* m, int* n, int* k, double* A, int* lda, double* tau,
            double* work, int* lwork, int* info);

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
int dgesvd(char* jobU, char* jobVT, int* m, int* n, 
           double *A, int* lda, double* S, double* U, int* ldu, 
           double* VT, int* ldvt, double *work, int* lwork, int* info);

// DGELSY computes the minimum-norm solution to a real linear least
// squares problem:
//       minimize || A * X - B ||
// using a complete orthogonal factorization of A. A is an M-by-N
// matrix which may be rank-deficient. See LAPACK documentation for details.
void dgelsy(int* m, int* n, int* nrhs, double* A, int* lda, double* B, int* ldb, 
            int* jpvt, double* rcond, int* rank, double* work, int* lwork, 
            int* info);

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
