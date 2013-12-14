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

#include "core/linear_algebra.h"

// NOTE: All LAPACK functions are implemented within the library, not here.

void matrix_fprintf(double* matrix, int nr, int nc, FILE* stream)
{
  fprintf(stream, "[");
  for (int i = 0; i < nr; ++i)
  {
    for (int j = 0; j < nc; ++j)
      fprintf(stream, "%g ", matrix[nr*j+i]);
    fprintf(stream, "; ");
  }
  fprintf(stream, "]");
}

void vector_fprintf(double* vec, int nr, FILE* stream)
{
  fprintf(stream, "[");
  for (int i = 0; i < nr; ++i)
    fprintf(stream, "%g ", vec[i]);
  fprintf(stream, "]");
}

double matrix2_det(double* matrix)
{
  return matrix[0]*matrix[3] - matrix[1]*matrix[2];
}

double matrix3_det(double* matrix)
{
  return matrix[0]*(matrix[4]*matrix[8] - matrix[5]*matrix[7]) -
         matrix[3]*(matrix[1]*matrix[8] - matrix[2]*matrix[7]) + 
         matrix[6]*(matrix[1]*matrix[5] - matrix[2]*matrix[4]);
}

void solve_2x2(double* A, double* b, double* x)
{
  double b0 = b[0], b1 = b[1];
  double inv_det_A = 1.0 / matrix2_det(A);
  x[0] = inv_det_A * ( A[3]*b0 - A[2]*b1);
  x[1] = inv_det_A * (-A[1]*b0 + A[0]*b1);
}

void solve_3x3(double* A, double* b, double* x)
{
  double b0 = b[0], b1 = b[1], b2 = b[2];

  // x = Ainv * b.
  double inv_det_A = 1.0 / matrix3_det(A);

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

