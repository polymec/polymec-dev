#ifndef POLYMEC_TABLE_H
#define POLYMEC_TABLE_H

#include "core/unordered_map.h"

// A table is a sparse table (like a spreadsheet) that has rows consisting 
// of sparse columns, each of which contains an element. One defines a 
// table using DEFINE_TABLE(map_name, element_type).

// Interface for a type x_table_t (with element type E) defined with
// DEFINE_TABLE(x_table, E):
// 
// x_table_t* x_table_new() - Creates a new empty table.
// x_table_t* x_table_new_with_capacity(int N) - Creates a new table with 
//                                               initial capacity N (rows).
// x_table_t* x_table_copy(x_table_t* t) - Creates a (shallow) copy of t.
// void x_table_free(x_table_t* table) - Destroys the table.
// void x_table_clear(x_table_t* table) - Empties the table.
// bool x_table_get_row(x_table_t* table, int row, int* num_cols, int** cols, x_table_value_t** values) - Retrieves the given row in the table, returning true if successful, false if not.
// x_table_value_t* x_table_get(x_table_t* table, int row, int col) - Returns the element at the given row/column, or NULL.
// bool x_table_contains_row(x_table_t* table, int row) - Returns true if the table contains the given row, false if not.
// bool x_table_contains(x_table_t* table, int row, int col) - Returns true if the table contains an element at row/column, false if not.
// void x_table_insert_row(x_table_t* table, int row, int num_cols, int* cols, x_table_value_t* values) - Inserts a row into the table.
// void x_table_insert_row_with_dtor(x_table_t* table, int row, int num_cols, int* cols, x_table_value_t* values, x_table_value_dtor dtor) - Inserts a row into the table with a destructor.
// void x_table_insert(x_table_t* table, int row, int col, x_table_value_t value) - Inserts a value at the given row/column.
// void x_table_insert_with_dtor(x_table_t* table, int row, int col, x_table_value_t value, x_table_value_dtor dtor) - Inserts a value at the given row/column with a destructor.
// void x_table_copy_row(x_table_t* table, int src_row, int dest_row) - Copies the contents of src_row into dest_row.
// void x_table_delete_row(x_table_t* table, int row) - Deletes the given row from the table.
// void x_table_delete(x_table_t* table, int row, col) - Deletes the given element from the table.
// bool x_table_next_row(x_table_t* map, int* pos, int* row, int* num_cols, int** cols, x_table_value_t** values) - Allows the traversal of the table rows.

#define DEFINE_TABLE(table_name, value_type) \
typedef value_type table_name##_value_t; \
typedef void (*table_name##_value_dtor)(table_name##_value_t); \
typedef struct \
{ \
  int num_cols; \
  int* cols; \
  table_name##_value_t* values; \
  table_name##_value_dtor dtors; \
} table_name##_row_t; \
\
static inline table_name##_row_t* table_name##_row_new(int num_cols, int* cols, table_name##_value_t* values, table_name##_value_dtor* dtors) \
{ \
  ASSERT(num_cols > 0); \
  ASSERT(cols != NULL); \
  ASSERT(values != NULL); \
  table_name##_row_t* row = malloc(sizeof(table_name##_row_t)); \
  row->num_cols = num_cols; \
  row->cols = malloc(sizeof(int)*num_cols); \
  for (int c = 0; c < num_cols; ++c) \
  { \
    row->cols[c] = cols[c]; \
    row->values[c] = values[c]; \
    row->dtors[c] = (dtors != NULL) ? dtors[c] : NULL; \
  } \
  return row; \
} \
\
static inline void table_name##_row_free(table_name##_row_t* row) \
{ \
  for (int i = 0; i < row->num_cols; ++i) \
  { \
    if (row->dtors[i] != NULL) \
      row->dtors[i](row->values[i]); \
  } \
  free(row->cols); \
  free(row->values); \
  free(row->dtors); \
  free(row); \
} \
\
DEFINE_UNORDERED_MAP(table_name##_map_t, int, table_name##_row_t*, int_hash, int_equals) \
\
typedef struct \
{ \
  table_name##_map_t map; \
  int num_rows; \
} table_name##_t; \
\
static inline table_name##_t* table_name##_new_with_capacity(int N) \
{ \
  table_name##_t* table = malloc(sizeof(table_name##_t)); \
  table->map = table_name##_map_new_with_capacity(N); \
  table->num_rows = 0; \
  return table; \
} \
\
static inline table_name##_t* table_name##_new() \
{ \
  return table_name##_new_with_capacity(32); \
} \
\
static inline void table_name##_clear(table_name##_t* table) \
{ \
  table_name##_map_clear(table->map); \
  table->num_rows = 0; \
} \
\
static inline void table_name##_free(table_name##_t* table) \
{ \
  table_name##_clear(table); \
  free(table); \
} \
\
static inline bool table_name##_get_row(table_name##_t* table, int row, int* num_cols, int** cols, table_name##_value_t** values) \
{ \
  table_name##_row_t** row = table_name##_map_get(table->map, row); \
  if (row == NULL) \
    return false; \
  *num_cols = *row->num_cols; \
  *cols = *row->cols; \
  *values = *row->values; \
} \
\
static inline table_name##_value_t* table_name##_get(table_name##_t* table, int row, int col) \
{ \
  table_name##_row_t** row = table_name##_map_get(table->map, row); \
  if (row == NULL) \
    return NULL; \
  table_name##_row_t* r = *row; \
  for (int c = 0; c < r->num_cols; ++c) \
  { \
    if (r->cols[c] == col) \
      return &r->values[c]; \
  } \
  return NULL; \
} \
\
static inline bool table_name##_contains_row(table_name##_t* table, int row) \
{ \
  return table_name##_map_contains(table->map, row); \
} \
\
static inline bool table_name##_contains(table_name##_t* table, int row, int col) \
{ \
  return (table_name##_get(table, row, col) != NULL); \
} \
\
static inline void table_name##_insert_row_with_dtor(table_name##_t* table, int row, int num_cols, int* cols, table_name##_value_t* values, table_name##_value_dtor dtor) \
{ \
  table_name##_value_dtor* dtors; \
  if (dtor != NULL) \
  { \
    dtors = malloc(sizeof(table_name##_value_dtor)); \
    for (int c = 0; c < num_cols; ++c) \
      dtors[c] = dtor; \
  } \
  table_name##_row_t* r = table_name##_row_new(num_cols, cols, values, dtors); \
  table_name##_map_insert_with_dtor(table->map, row, r, table_name##_row_free); \
} \
\
static inline void table_name##_insert_row(table_name##_t* table, int row, int num_cols, int* cols, table_name##_value_t* values) \
{ \
  table_name##_insert_row_with_dtor(table, row, num_cols, cols, values, NULL); \
} \
\
static inline void table_name##_insert_with_dtor(table_name##_t* table, int row, int col, table_name##_value_t value, table_name##_value_dtor dtor) \
{ \
} \
static inline void table_name##_insert(table_name##_t* map, int row, int col, table_name##_value_t value) \
{ \
  table_name##_insert_with_dtor(map, row, col, value, NULL); \
} \
\
static inline void table_name##_copy_row(table_name##_t* table, int src_row, int dest_row) \
{ \
  table_name##_row_t** row = table_name##_map_get(table->map, src_row); \
  if (row == NULL) return; \
  table_name##_row_t* r = *row; \
  table_name##_insert_row(table, dest_row, r->num_cols, r->cols, r->values); \
} \
\
static inline void table_name##_delete_row(table_name##_t* table, int row) \
{ \
  table_name##_map_delete(table->map, row); \
} \
\
static inline void table_name##_delete(table_name##_t* table, int row, int col) \
{ \
  table_name##_row_t** row = table_name##_map_get(table->map, src_row); \
  if (row == NULL) return; \
  table_name##_row_delete_col(*row, col); \
} \
\
static inline bool table_name##_next_row(table_name##_t* table, int* pos, int* row, int* num_cols, int** cols, table_name##_value_t** values) \
{ \
  table_name##_row_t* r; \
  bool result = table_name##_map_next(table->map, pos, row, &r); \
  if (result) \
  { \
    *num_cols = r->num_cols; \
    *cols = r->cols; \
    *values = r->values; \
  } \
  return result; \
} \
\
static inline table_name##_t* table_name##_copy(table_name##_t* table) \
{ \
  table_name##_t* t = malloc(sizeof(table_name##_t)); \
  t->map = map_name##_copy(table->map); \
  t->num_rows = table->num_rows; \
  return t; \
} \

// Define some tables.
DEFINE_UNORDERED_MAP(int_table, int)
DEFINE_UNORDERED_MAP(double_table, double)

#endif
