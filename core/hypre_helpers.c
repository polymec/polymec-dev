#include "core/hypre_helpers.h"

#ifdef __cplusplus
extern "C" {
#endif

HYPRE_IJMatrix HYPRE_IJMatrixNew(index_space_t* index_space)
{
  HYPRE_IJMatrix A;
  MPI_Comm comm = index_space->comm;
  int low = index_space->low;
  int high = index_space->high;
  int err = HYPRE_IJMatrixCreate(comm, low, high, low, high, &A);
  ASSERT(err == 0);
  err = HYPRE_IJMatrixSetObjectType(A, HYPRE_PARCSR);
  ASSERT(err == 0);
  return A;
}

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

HYPRE_IJVector HYPRE_IJVectorNew(index_space_t* index_space)
{
  HYPRE_IJVector x;
  MPI_Comm comm = index_space->comm;
  int low = index_space->low;
  int high = index_space->high;
  int err = HYPRE_IJVectorCreate(comm, low, high, &x);
  ASSERT(err == 0);
  err = HYPRE_IJVectorSetObjectType(x, HYPRE_PARCSR);
  ASSERT(err == 0);

  // Initialize the vector to zero.
  int N = index_space->high - index_space->low;
  int indices[N];
  double values[N];
  for (int i = 0; i < N; ++i)
  {
    indices[i] = index_space->low + i;
    values[i] = 0.0;
  }
  err = HYPRE_IJVectorInitialize(x);
  ASSERT(err == 0);
  HYPRE_IJVectorSetValues(x, N, indices, values);
  err = HYPRE_IJVectorAssemble(x);
  ASSERT(err == 0);

  return x;
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

