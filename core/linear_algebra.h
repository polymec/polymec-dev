#ifndef POLYMEC_LINEAR_ALGEBRA_H
#define POLYMEC_LINEAR_ALGEBRA_H

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

// Some LAPACK prototypes.

// LU factorization.
void dgetrf(int *N, int *NRHS, double *A, int *LDA, int *IPIV, int *INFO);

// Solution of a linear system using an LU factorization.
void dgetrs(char *TRANS, int *N, int *NRHS, double *A, 
            int *LDA, int *IPIV, double *B, int *LDB, int *INFO);

// Matrix-vector multiplication: y := alpha*A*x + beta*y.
void dgemv(char *trans, int *m, int *n, double *alpha,
           void *a, int *lda, void *x, int *incx,
           double *beta, void *y, int *incy);

// Matrix-matrix multiplication: C := alpha*op(A)*op(B) + beta*C.
void dgemm(char *transa, char* transB, int *m, int *n, int *k, double *alpha,
           double *A, int *lda, double *B, int *ldb, double *beta, double *C, 
           int *ldc);

// Print a (column-major-ordered) matrix to the given file stream.
void matrix_fprintf(double* matrix, int nr, int nc, FILE* stream);

// Print a vector to the given file stream.
void vector_fprintf(double* vec, int nr, FILE* stream);

// Computes the determinant of the given 3x3 matrix.
double matrix3_det(double* matrix);

#ifdef __cplusplus
}
#endif

#endif
