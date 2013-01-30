#include "core/hypre_helpers.h"

#ifdef __cplusplus
extern "C" {
#endif

void HYPRE_IJMatrixSetRowSizesFromTable(HYPRE_IJMatrix matrix, index_space_t* index_space, double_table_t* table)
{
  int N = index_space->high - index_space->low;
  int nnz[N];
  int rpos = 0, i;
  double_table_row_t* row_data;
  while (double_table_next_row(table, &rpos, &i, &row_data))
    nnz[i] = row_data->size;
  HYPRE_IJMatrixSetRowSizes(matrix, nnz);
}

static void HYPRE_IJMatrixModifyValuesFromTable(HYPRE_IJMatrix matrix, index_space_t* index_space, double_table_t* table, int (*modify_values)(HYPRE_IJMatrix, int, int*, const int*, const int*, const double*))
{
  HYPRE_IJMatrixInitialize(matrix);
  int rpos = 0, i;
  double_table_row_t* row_data;
  while (double_table_next_row(table, &rpos, &i, &row_data))
  {
    int cpos = 0, j, num_cols = row_data->size;
    int indices[num_cols], offset = 0;
    double values[num_cols], Aij;
    while (double_table_row_next(row_data, &cpos, &j, &Aij))
    {
      indices[offset] = j;
      values[offset] = Aij;
      offset++;
    }
    modify_values(matrix, 1, &num_cols, &i, indices, values);
  }
  HYPRE_IJMatrixAssemble(matrix);
}

void HYPRE_IJMatrixSetValuesFromTable(HYPRE_IJMatrix matrix, index_space_t* index_space, double_table_t* table)
{
  HYPRE_IJMatrixModifyValuesFromTable(matrix, index_space, table, HYPRE_IJMatrixSetValues);
}

void HYPRE_IJMatrixAddToValuesFromTable(HYPRE_IJMatrix matrix, index_space_t* index_space, double_table_t* table)
{
  HYPRE_IJMatrixModifyValuesFromTable(matrix, index_space, table, HYPRE_IJMatrixAddToValues);
}

static void HYPRE_IJVectorModifyValuesFromArray(HYPRE_IJVector vector, index_space_t* index_space, double* array, int (*modify_values)(HYPRE_IJVector, int, const int*, const double*))
{
  HYPRE_IJVectorInitialize(vector);
  int N = index_space->high - index_space->low;
  int indices[N];
  for (int i = 0; i < N; ++i)
    indices[i] = index_space->low + i;
  modify_values(vector, N, indices, array);
  HYPRE_IJVectorAssemble(vector);
}

void HYPRE_IJVectorSetValuesFromArray(HYPRE_IJVector vector, index_space_t* index_space, double* array)
{
  HYPRE_IJVectorModifyValuesFromArray(vector, index_space, array, HYPRE_IJVectorSetValues);
}

void HYPRE_IJVectorAddToValuesFromArray(HYPRE_IJVector vector, index_space_t* index_space, double* array)
{
  HYPRE_IJVectorModifyValuesFromArray(vector, index_space, array, HYPRE_IJVectorAddToValues);
}

#ifdef __cplusplus
}
#endif

