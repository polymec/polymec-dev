#include "core/linear_algebra.h"

#ifdef __cplusplus
extern "C" {
#endif

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

double matrix3_det(double* matrix)
{
  return matrix[0]*(matrix[4]*matrix[8] - matrix[5]*matrix[7]) -
         matrix[3]*(matrix[1]*matrix[8] - matrix[2]*matrix[7]) + 
         matrix[6]*(matrix[1]*matrix[5] - matrix[2]*matrix[4]);
}

#ifdef __cplusplus
}
#endif

