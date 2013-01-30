#ifndef POLYMEC_HYPRE_HELPERS_H
#define POLYMEC_HYPRE_HELPERS_H

#include "HYPRE_IJ_mv.h"
#include "HYPRE_krylov.h"
#include "HYPRE_parcsr_mv.h"
#include "HYPRE_parcsr_ls.h"
#include "core/index_space.h"
#include "core/table.h"

// These helper functions obey HYPRE's style conventions and are meant to 
// add support for Polymec's ansilary data structures.

#ifdef __cplusplus
extern "C" {
#endif

// Preallocates the nonzero entries in the given IJ matrix using the 
// sparsity pattern in the given table.
void HYPRE_IJMatrixSetRowSizesFromTable(HYPRE_IJMatrix matrix, index_space_t* index_space, double_table_t* table);

// Inserts the values from the given table into the given IJ matrix.
void HYPRE_IJMatrixSetValuesFromTable(HYPRE_IJMatrix matrix, index_space_t* index_space, double_table_t* table);

// Adds the values from the given table into the given IJ matrix.
void HYPRE_IJMatrixAddToValuesFromTable(HYPRE_IJMatrix matrix, index_space_t* index_space, double_table_t* table);

#ifdef __cplusplus
}
#endif

#endif


