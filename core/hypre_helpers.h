// Copyright 2012-2013 Jeffrey Johnson.
// 
// This file is part of Polymec, and is licensed under the Apache License, 
// Version 2.0 (the "License"); you may not use this file except in 
// compliance with the License. You may may find the text of the license in 
// the LICENSE file at the top-level source directory, or obtain a copy of 
// it at
// 
// http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef POLYMEC_HYPRE_HELPERS_H
#define POLYMEC_HYPRE_HELPERS_H

#include "HYPRE.h"
#include "HYPRE_IJ_mv.h"
#include "HYPRE_krylov.h"
#include "HYPRE_parcsr_mv.h"
#include "HYPRE_parcsr_ls.h"
#include "HYPRE_utilities.h"
#include "core/index_space.h"
#include "core/table.h"

// These helper functions obey HYPRE's style conventions and are meant to 
// add support for Polymec's ansilary data structures.

// Creates a HYPRE IJ matrix using the given index space.
HYPRE_IJMatrix HYPRE_IJMatrixNew(index_space_t* index_space);

// Preallocates the nonzero entries in the given IJ matrix using the 
// sparsity pattern in the given table.
void HYPRE_IJMatrixSetRowSizesFromTable(HYPRE_IJMatrix matrix, index_space_t* index_space, double_table_t* table);

// Inserts the values from the given table into the given IJ matrix.
void HYPRE_IJMatrixSetValuesFromTable(HYPRE_IJMatrix matrix, index_space_t* index_space, double_table_t* table);

// Adds the values from the given table into the given IJ matrix.
void HYPRE_IJMatrixAddToValuesFromTable(HYPRE_IJMatrix matrix, index_space_t* index_space, double_table_t* table);

// Creates a HYPRE IJ vector using the given index space.
HYPRE_IJVector HYPRE_IJVectorNew(index_space_t* index_space);

// Inserts the values from the given table into the given IJ vector.
void HYPRE_IJVectorSetValuesFromArray(HYPRE_IJVector vector, index_space_t* index_space, double* array);

// Adds the values from the given table into the given IJ vector.
void HYPRE_IJVectorAddToValuesFromArray(HYPRE_IJVector vector, index_space_t* index_space, double* array);

// Gets values from the given vector, placing them into the given array.
void HYPRE_IJVectorGetValuesToArray(HYPRE_IJVector vector, index_space_t* index_space, double* array);

#endif


