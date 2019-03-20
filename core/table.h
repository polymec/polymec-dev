// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_TABLE_H
#define POLYMEC_TABLE_H

#include "core/unordered_map.h"

/// \addtogroup core core
///@{

/// \def DEFINE_TABLE(map_name, element_type).
/// Defines a table with a given element type. A table is a sparse table
/// (like a spreadsheet) that has rows consisting of sparse columns, each of
/// which contains a primitive element. Elements in a table must be primitive
/// data types. The following interface is defined for a table with name
/// `x_table`.
/// * `x_table_t* x_table_new()` - Creates a new empty table.
/// * `x_table_t* x_table_new_with_capacity(int N)` - Creates a new table with initial capacity N (rows).
/// * `x_table_t* x_table_clone(x_table_t* t)` - Creates a deep copy of t.
/// * `void x_table_free(x_table_t* table)` - Destroys the table.
/// * `void x_table_clear(x_table_t* table)` - Empties the table.
/// * `x_table_row_t* x_table_get_row(x_table_t* table, int row)` - Returns the row in the given table (itself an unordered map).
/// * `x_table_value_t* x_table_get(x_table_t* table, int row, int col)` - Returns the element at the given row/column, or NULL.
/// * `bool x_table_contains_row(x_table_t* table, int row)` - Returns true if the table contains the given row, false if not.
/// * `bool x_table_contains(x_table_t* table, int row, int col)` - Returns true if the table contains an element at row/column, false if not.
/// * `void x_table_insert_row(x_table_t* table, int row, int num_cols, int* cols, x_table_value_t* values)` - Inserts a row into the table.
/// * `void x_table_insert_row_with_dtor(x_table_t* table, int row, int num_cols, int* cols, x_table_value_t* values, x_table_value_dtor dtor)` - Inserts a row into the table with a destructor.
/// * `void x_table_insert(x_table_t* table, int row, int col, x_table_value_t value)` - Inserts a value at the given row/column.
/// * `void x_table_insert_with_dtor(x_table_t* table, int row, int col, x_table_value_t value, x_table_value_dtor dtor)` - Inserts a value at the given row/column with a destructor.
/// * `void x_table_delete_row(x_table_t* table, int row)` - Deletes the given row from the table.
/// * `void x_table_delete(x_table_t* table, int row, col)` - Deletes the given element from the table.
/// * `bool x_table_next_row(x_table_t* table, int* pos, int* row, x_table_row_t** row_data)` - Allows the traversal of the table rows.
/// * `bool x_table_next_cell(x_table_t* table, x_table_cell_pos_t* pos, int* row, int* col, x_table_value_t* value)` - Allows the traversal of the table values.
/// * `x_table_cell_pos_t x_table_start(x_table_t* table)` - Returns a new value position for use with x_table_next.
///
/// Member data for a table `table`:
/// * `table->num_rows` - The number of rows in the table.
///
/// \param table_name The name of the table.
/// \param value_type The data type stored in the table.
#define DEFINE_TABLE(table_name, value_type) \
typedef value_type table_name##_value_t; \
typedef void (*table_name##_value_dtor)(table_name##_value_t); \
DEFINE_UNORDERED_MAP(table_name##_row, int, value_type, int_hash, int_equals) \
DEFINE_UNORDERED_MAP(table_name##_map, int, table_name##_row_t*, int_hash, int_equals) \
typedef struct \
{ \
  int row_pos; \
  int current_row; \
  int col_pos; \
} table_name##_cell_pos_t; \
\
typedef struct \
{ \
  table_name##_map_t* map; \
  int num_rows; \
} table_name##_t; \
\
static inline table_name##_t* table_name##_new_with_capacity(int N) \
{ \
  table_name##_t* table = (table_name##_t*)polymec_malloc(sizeof(table_name##_t)); \
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
  table_name##_map_free(table->map); \
  polymec_free(table); \
} \
\
static inline table_name##_row_t** table_name##_get_row(table_name##_t* table, int row) \
{ \
  return table_name##_map_get(table->map, row); \
} \
\
static inline table_name##_value_t* table_name##_get(table_name##_t* table, int row, int col) \
{ \
  table_name##_row_t** row_data = table_name##_map_get(table->map, row); \
  if (row_data == NULL) \
    return NULL; \
  table_name##_row_t* r = *row_data; \
  return table_name##_row_get(r, col); \
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
  table_name##_row_t* r = table_name##_row_new(); \
  for (int c = 0; c < num_cols; ++c) \
    table_name##_row_insert_with_v_dtor(r, cols[c], values[c], dtor); \
  table_name##_map_insert_with_v_dtor(table->map, row, r, table_name##_row_free); \
  table->num_rows = table->map->size; \
} \
\
static inline void table_name##_insert_row(table_name##_t* table, int row, int num_cols, int* cols, table_name##_value_t* values) \
{ \
  table_name##_insert_row_with_dtor(table, row, num_cols, cols, values, NULL); \
} \
\
static inline void table_name##_insert_with_dtor(table_name##_t* table, int row, int col, table_name##_value_t value, table_name##_value_dtor dtor) \
{ \
  table_name##_row_t** row_data = table_name##_map_get(table->map, row); \
  table_name##_row_t* r; \
  if (row_data != NULL) \
  { \
    r = *row_data; \
  } \
  else \
  { \
    r = table_name##_row_new(); \
    table_name##_map_insert_with_v_dtor(table->map, row, r, table_name##_row_free); \
    table->num_rows = table->map->size; \
  } \
  table_name##_row_insert_with_v_dtor(r, col, value, dtor); \
} \
\
static inline void table_name##_insert(table_name##_t* map, int row, int col, table_name##_value_t value) \
{ \
  table_name##_insert_with_dtor(map, row, col, value, NULL); \
} \
\
static inline void table_name##_delete_row(table_name##_t* table, int row) \
{ \
  table_name##_map_delete(table->map, row); \
  table->num_rows = table->map->size; \
} \
\
static inline void table_name##_delete(table_name##_t* table, int row, int col) \
{ \
  table_name##_row_t** r = table_name##_map_get(table->map, row); \
  if (r == NULL) return; \
  table_name##_row_delete(*r, col); \
  if ((*r)->size == 0) \
    table_name##_delete_row(table, row); \
} \
\
static inline bool table_name##_next_row(table_name##_t* table, int* pos, int* row, table_name##_row_t** row_data) \
{ \
  return table_name##_map_next(table->map, pos, row, row_data); \
} \
\
static inline bool table_name##_next_cell(table_name##_t* table, table_name##_cell_pos_t* pos, int* row, int* col, table_name##_value_t* value) \
{ \
  table_name##_row_t* row_data; \
  if (pos->row_pos == 0) \
  { \
    if (!table_name##_map_next(table->map, &pos->row_pos, row, &row_data)) \
      return false; \
    pos->current_row = *row; \
  } \
  else \
  { \
    table_name##_row_t** r = table_name##_get_row(table, pos->current_row); \
    if (r == NULL) return false; \
    row_data = *r; \
  } \
  if (table_name##_row_next(row_data, &pos->col_pos, col, value)) \
    return true; \
  else \
  { \
    if (!table_name##_map_next(table->map, &pos->row_pos, row, &row_data)) \
      return false; \
    pos->current_row = *row; \
    pos->col_pos = 0; \
    bool result = table_name##_row_next(row_data, &pos->col_pos, col, value); \
    ASSERT(result == true); \
    return result; \
  } \
} \
\
static inline table_name##_t* table_name##_clone(table_name##_t* table) \
{ \
  table_name##_t* t = (table_name##_t*)polymec_malloc(sizeof(table_name##_t)); \
  t->map = table_name##_map_clone(table->map, NULL, NULL, NULL, NULL); \
  t->num_rows = table->num_rows; \
  return t; \
} \
static inline table_name##_cell_pos_t table_name##_start(table_name##_t* table) \
{ \
  table_name##_cell_pos_t pos = {0, 0, 0}; \
  return pos; \
} \

///@}

// Define some tables.
DEFINE_TABLE(int_table, int)
DEFINE_TABLE(real_table, real_t)

#endif
