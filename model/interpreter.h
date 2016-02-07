// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_INTERPRETER_H
#define POLYMEC_INTERPRETER_H

#include "core/polymec.h"
#include "core/mesh.h"
#include "core/st_func.h"
#include "core/unordered_map.h"

// Forward declaration for Lua innards.
struct lua_State;

// This class interprets input files and stores values of variables therein.
typedef struct interpreter_t interpreter_t;

// This defines all valid input variable types.
typedef enum
{
  INTERPRETER_STRING,
  INTERPRETER_NUMBER,
  INTERPRETER_BOOLEAN,
  INTERPRETER_MESH,
  INTERPRETER_STRING_LIST,
  INTERPRETER_POINT,
  INTERPRETER_POINT_LIST,
  INTERPRETER_VECTOR,
  INTERPRETER_BOUNDING_BOX,
  INTERPRETER_SCALAR_FUNCTION,
  INTERPRETER_VECTOR_FUNCTION,
  INTERPRETER_SYM_TENSOR_FUNCTION,
  INTERPRETER_TENSOR_FUNCTION,
  INTERPRETER_SEQUENCE,
  INTERPRETER_TABLE,
  INTERPRETER_USER_DEFINED, // <-- on your head be it!!
  INTERPRETER_TERMINUS // Used only to terminate validation lists.
} interpreter_var_type_t;

// The interpreter makes no distinction between how lists of vectors 
// and points are stored, so this is provided as a convenience.
#define INTERPRETER_VECTOR_LIST INTERPRETER_POINT_LIST

typedef enum
{
  REQUIRED,
  OPTIONAL
} interpreter_required_t;

// This is used to validate variables found within input files. It is a 
// variable/type/required triple that associates variables names with their 
// intended types. A third field indicates whether the variable is 
// REQUIRED or OPTIONAL in the input. Parsing of the input will not succeed 
// if a REQUIRED variable is not found.
typedef struct 
{
  char* variable;                       // Name of the variable.
  interpreter_var_type_t type;          // Type of the variable.
  interpreter_required_t required;      // REQUIRED or OPTIONAL.
} interpreter_validation_t;

// This validation object is used to terminate lists of valid inputs.
static const interpreter_validation_t END_OF_VALID_INPUTS = {(char*)"TERMINUS", INTERPRETER_TERMINUS, REQUIRED};

// The docstring type holds lines of documentation for interpreter functions.
// Docstrings are stored in a static database and destroyed upon exit, so they 
// do not have callable destructors.
typedef struct docstring_t docstring_t;

// Creates a new, empty docstring.
docstring_t* docstring_new();

// Creates a new docstring from the given string s, breaking it into lines 
// at each newline.
docstring_t* docstring_from_string(const char* s);

// Appends a line of text to the given docstring.
void docstring_append(docstring_t* docs, const char* line);

// Allows the retrieval of lines of text from a docstring, returning true if 
// a line was retrieved and false if not. Set pos to 0 to reset the traversal.
bool docstring_next(docstring_t* docs, int* pos, char** line);

// Returns an internal pointer to the first line of the docstring. This is 
// the empty string if the docstring is empty.
char* docstring_first_line(docstring_t* docs);

// Creates a new interpreter for use with a model. The interpreter uses the 
// given set of types to validate variables in input files.
interpreter_t* interpreter_new(interpreter_validation_t* valid_inputs);

// Destroys the given interpreter.
void interpreter_free(interpreter_t* interp);

// Register a user-defined function with the interpreter.
// A docstring may be associated with the function, and is consumed if given.
void interpreter_register_function(interpreter_t* interp, const char* function_name, int (*function)(struct lua_State*), docstring_t* doc);

// Returns true if given user-defined global table has been registered with 
// the interpreter, false if not.
bool interpreter_has_global_table(interpreter_t* interp, const char* table_name);

// Register an empty user-defined global table with the interpreter.
// A docstring may be associated with the function, and is consumed if given.
void interpreter_register_global_table(interpreter_t* interp, const char* table_name, docstring_t* doc);

// Registers a method within the given global table in the interpreter.
// A docstring may be associated with the function, and is consumed if given.
void interpreter_register_global_method(interpreter_t* interp, const char* table_name, const char* method_name, int (*method)(struct lua_State*), docstring_t* doc);

// This function writes documentation for the given registered entity 
// (a function or global symbol) within the interpreter to the given stream, 
// or a generic "not found" string if no such entity exists.
void interpreter_help(interpreter_t* interp, const char* entity, FILE* stream);

// Parses the input file, storing the values in the interpreter.
void interpreter_parse_file(interpreter_t* interp, char* input_file);

// Parses the input string, storing values in the interpreter.
void interpreter_parse_string(interpreter_t* interp, char* input_string);

// Returns true if the given variable exists in the interpreter and 
// matches the given type, false otherwise.
bool interpreter_contains(interpreter_t* interp, const char* variable, interpreter_var_type_t type); 

// Fetches the given string from the interpreter, returning NULL if it 
// is not found or if it is not a string. The string refers to storage 
// within the interpreter and should not be freed.
char* interpreter_get_string(interpreter_t* interp, const char* name);

// Sets the given variable within the interpreter to the given string 
// value. Any existing value of this variable is overwritten, and the 
// string is copied to the interpreter.
void interpreter_set_string(interpreter_t* interp, const char* name, const char* value);

// Fetches the given number from the interpreter, returning -FLT_MAX if it 
// is not found or if it is not a number.
real_t interpreter_get_number(interpreter_t* interp, const char* name);

// Sets the given variable within the interpreter to the given numeric 
// value. Any existing value of this variable is overwritten.
void interpreter_set_number(interpreter_t* interp, const char* name, real_t value);

// Fetches the given boolean value from the interpreter, returning false if it 
// is not found or if it is not a boolean.
bool interpreter_get_boolean(interpreter_t* interp, const char* name);

// Sets the given variable within the interpreter to the given boolean 
// value. Any existing value of this variable is overwritten.
void interpreter_set_boolean(interpreter_t* interp, const char* name, bool value);

// Fetches the given point from the interpreter, returning NULL if it 
// is not found or if it is not a point.
point_t* interpreter_get_point(interpreter_t* interp, const char* name);

// Sets the given variable within the interpreter to the given point
// value. Any existing value of this variable is overwritten.
void interpreter_set_point(interpreter_t* interp, const char* name, point_t* value);

// Fetches the given list of points from the interpreter, returning NULL if it 
// is not found or if it is not a point list. num_points will store the number 
// of points in the list upon successful completion.
point_t* interpreter_get_pointlist(interpreter_t* interp, const char* name, int* num_points);

// Sets the given variable within the interpreter to the given list of points.
// Any existing value of this variable is overwritten.
void interpreter_set_pointlist(interpreter_t* interp, const char* name, point_t* points, int num_points);

// Fetches the given vector from the interpreter, returning NULL if it 
// is not found or if it is not a vector.
vector_t* interpreter_get_vector(interpreter_t* interp, const char* name);

// Sets the given variable within the interpreter to the given point
// value. Any existing value of this variable is overwritten.
void interpreter_set_vector(interpreter_t* interp, const char* name, vector_t* value);

// Fetches the given list of vectors from the interpreter, returning NULL if it 
// is not found or if it is not a vector list. num_vectors will store the number 
// of vectors in the list upon successful completion.
vector_t* interpreter_get_vectorlist(interpreter_t* interp, const char* name, int* num_vectors);

// Sets the given variable within the interpreter to the given list of vectors.
// Any existing value of this variable is overwritten.
void interpreter_set_vectorlist(interpreter_t* interp, const char* name, vector_t* vectors, int num_vectors);

// Fetches the given bounding box from the interpreter, returning NULL if it 
// is not found or if it is not a bounding box.
bbox_t* interpreter_get_bbox(interpreter_t* interp, const char* name);

// Sets the given variable within the interpreter to the given bounding box.
// Any existing value of this variable is overwritten.
void interpreter_set_bbox(interpreter_t* interp, const char* name, bbox_t* value);

// Fetches the given mesh from the interpreter, returning NULL if it 
// is not found or if it is not a mesh. The caller assumes responsibility
// for destroying the mesh after this call.
mesh_t* interpreter_get_mesh(interpreter_t* interp, const char* name);

// Sets the given variable within the interpreter to the given mesh 
// object. Any existing value of this variable is overwritten.
// Control of the mesh is assumed by the interpreter.
void interpreter_set_mesh(interpreter_t* interp, const char* name, mesh_t* value);

// Fetches the given scalar function from the interpreter, returning NULL 
// if it is not found or if it is not such a function.
st_func_t* interpreter_get_scalar_function(interpreter_t* interp, const char* name);

// Sets the given variable within the interpreter to the given scalar function 
// object. Any existing value of this variable is overwritten.
void interpreter_set_scalar_function(interpreter_t* interp, const char* name, st_func_t* value);

// Fetches the given vector function from the interpreter, returning NULL 
// if it is not found or if it is not such a function.
st_func_t* interpreter_get_vector_function(interpreter_t* interp, const char* name);

// Sets the given variable within the interpreter to the given vector function 
// object. Any existing value of this variable is overwritten.
void interpreter_set_vector_function(interpreter_t* interp, const char* name, st_func_t* value);

// Fetches the given (6-component) symmetric tensor function from the 
// interpreter, returning NULL if it is not found or if it is not such a function.
// A symmetric tensor function returns components (Sxx, Sxy, Sxz, Syy, Syz, Szz).
st_func_t* interpreter_get_sym_tensor_function(interpreter_t* interp, const char* name);

// Sets the given variable within the interpreter to the given (6-component)
// symmetric tensor function object. Any existing value of this variable is overwritten.
// A symmetric tensor function returns components (Sxx, Sxy, Sxz, Syy, Syz, Szz).
void interpreter_set_sym_tensor_function(interpreter_t* interp, const char* name, st_func_t* value);

// Fetches the given (9-component) tensor function from the 
// interpreter, returning NULL if it is not found or if it is such a function.
// A tensor function returns components (Sxx, Sxy, Sxz, Syx, Syy, Syz, Szx, Szy, Szz).
st_func_t* interpreter_get_tensor_function(interpreter_t* interp, const char* name);

// Sets the given variable within the interpreter to the given (6-component)
// symmetric tensor function object. Any existing value of this variable is overwritten.
// A tensor function returns components (Sxx, Sxy, Sxz, Syx, Syy, Syz, Szx, Szy, Szz).
void interpreter_set_tensor_function(interpreter_t* interp, const char* name, st_func_t* value);

// Fetches the given list of strings from the interpreter, returning NULL 
// if it is not found or if it is not a list of strings. The size of the list is 
// stored in *len. The caller assumes responsibility for destroying the 
// list of strings after this call.
char** interpreter_get_stringlist(interpreter_t* interp, const char* name, int* len);

// Sets the given variable within the interpreter to the given list of
// strings. Any existing value of this variable is overwritten.
void interpreter_set_stringlist(interpreter_t* interp, const char* name, char** list, int len);

// Fetches the given sequence of numbers from the interpreter, returning NULL 
// if it is not found or if it is not a sequence. The size of the sequence is 
// stored in *len. The caller assumes responsibility for destroying the 
// sequence after this call.
real_t* interpreter_get_sequence(interpreter_t* interp, const char* name, int* len);

// Sets the given variable within the interpreter to the given sequence of
// numbers. Any existing value of this variable is overwritten.
void interpreter_set_sequence(interpreter_t* interp, const char* name, real_t* sequence, int len);

// Fetches the given table from the interpreter, returning NULL if it 
// is not found or if it is not a table. The caller assumes responsibility
// for destroying the table after this call.
string_ptr_unordered_map_t* interpreter_get_table(interpreter_t* interp, const char* name);

// Sets the given variable within the interpreter to the given table of
// objects. Any existing value of this variable is overwritten.
void interpreter_set_table(interpreter_t* interp, const char* name, string_ptr_unordered_map_t* value);

// Fetches the given user-defined object from the interpreter, returning NULL 
// if it is not found or if it is not a user-defined object. The caller 
// assumes responsibility for destroying the user-defined object after 
// this call.
void* interpreter_get_user_defined(interpreter_t* interp, const char* name);

// Sets the given variable within the interpreter to the given pointer.
// Any existing value of this variable is overwritten.
// The interpreter assumes responsibility for destroying the pointer using 
// the supplied destructor (dtor).
void interpreter_set_user_defined(interpreter_t* interp, const char* name, void* value, void (*dtor)(void*));

//------------------------------------------------------------------------
// These methods can be used to create functions that extend an 
// interpreter.
//------------------------------------------------------------------------

// This helper returns true if the object at the given index is a sequence
// of numbers, false if not.
bool lua_issequence(struct lua_State* lua, int index);

// This helper retrieves a sequence from the given index on an active lua 
// interpreter, or returns NULL if the index does not point to a sequence.
real_t* lua_tosequence(struct lua_State* lua, int index, int* len);

// Pushes a sequence onto the interpreter's stack (as a return value for a 
// function).
void lua_pushsequence(struct lua_State* lua, real_t* sequence, int len);

// This helper returns true if the object at the given index is a list of
// strings, false if not.
bool lua_isstringlist(struct lua_State* lua, int index);

// This helper retrieves a list of strings from the given index on an active 
// lua interpreter, or returns NULL if the index does not point to a string list.
char** lua_tostringlist(struct lua_State* lua, int index, int* len);

// Pushes a string list onto the interpreter's stack (as a return value for a 
// function).
void lua_pushstringlist(struct lua_State* lua, char** list, int len);

// This helper returns true if the object at the given index is a point
// in 3-dimensional space, false if not.
bool lua_ispoint(struct lua_State* lua, int index);

// This helper retrieves a point from the given index on an active lua 
// interpreter, or returns NULL if the index does not point to a point.
point_t* lua_topoint(struct lua_State* lua, int index);

// Pushes a point onto the interpreter's stack (as a return value for a 
// function).
void lua_pushpoint(struct lua_State* lua, point_t* point);

// This helper returns true if the object at the given index is a list 
// (array) of points in 3-dimensional space, false if not.
bool lua_ispointlist(struct lua_State* lua, int index);

// This helper retrieves a list (array) of points from the given index on an 
// active lua interpreter, or returns NULL if the index does not point to a 
// point list. Here, the size argument stores the number of points in the 
// list.
point_t* lua_topointlist(struct lua_State* lua, int index, int* size);

// Pushes a list (array) of points onto the interpreter's stack (as a return 
// value for a function).
void lua_pushpointlist(struct lua_State* lua, point_t* points, int size);

// This helper returns true if the object at the given index is a vector
// in 3-dimensional space, false if not.
bool lua_isvector(struct lua_State* lua, int index);

// This helper retrieves a vector from the given index on an active lua 
// interpreter, or returns NULL if the index does not point to a point.
vector_t* lua_tovector(struct lua_State* lua, int index);

// Pushes a vector onto the interpreter's stack (as a return value for a 
// function).
void lua_pushvector(struct lua_State* lua, vector_t* vec);

// This helper returns true if the object at the given index is a list 
// (array) of vectors in 3-dimensional space, false if not.
bool lua_isvectorlist(struct lua_State* lua, int index);

// This helper retrieves a list (array) of vectors from the given index on an 
// active lua interpreter, or returns NULL if the index does not point to a 
// point list. Here, the size argument stores the number of vectors in the 
// list.
vector_t* lua_tovectorlist(struct lua_State* lua, int index, int* size);

// Pushes a list (array) of vectors onto the interpreter's stack (as a return 
// value for a function).
void lua_pushvectorlist(struct lua_State* lua, vector_t* vectors, int size);

// This helper returns true if the object at the given index is a bounding
// box, false if not.
bool lua_isboundingbox(struct lua_State* lua, int index);

// This helper retrieves a bounding box from the given index on an 
// active lua interpreter, or returns NULL if the index does not point to 
// such a thing.
bbox_t* lua_toboundingbox(struct lua_State* lua, int index);

// Pushes a bounding box onto the interpreter's stack (as a 
// return value for a function).
void lua_pushboundingbox(struct lua_State* lua, bbox_t* bbox);

// This helper returns true if the object at the given index is or can be 
// interpreted as a scalar-valued function, false if not.
bool lua_isscalarfunction(struct lua_State* lua, int index);

// This helper retrieves a scalar-valued function from the given index on an 
// active lua interpreter, or returns NULL if the object at the index cannot 
// be interpreted as one.
st_func_t* lua_toscalarfunction(struct lua_State* lua, int index);

// Pushes a scalar-valued function onto the interpreter's stack (as a 
// return value for a function).
void lua_pushscalarfunction(struct lua_State* lua, st_func_t* func);

// This helper returns true if the object at the given index is a vector-valued
// function, false if not.
bool lua_isvectorfunction(struct lua_State* lua, int index);

// This helper retrieves a vector-valued function from the given index on an 
// active lua interpreter, or returns NULL if the index does not point to an st_func.
st_func_t* lua_tovectorfunction(struct lua_State* lua, int index);

// Pushes a vector-valued function onto the interpreter's stack (as a 
// return value for a function).
void lua_pushvectorfunction(struct lua_State* lua, st_func_t* func);

// This helper returns true if the object at the given index is or can be 
// interpreted as a symmetric-tensor-valued function, false if not.
bool lua_issymtensorfunction(struct lua_State* lua, int index);

// This helper retrieves a symmetric-tensor-valued function from the given index on an 
// active lua interpreter, or returns NULL if the index does not point to an st_func.
st_func_t* lua_tosymtensorfunction(struct lua_State* lua, int index);

// Pushes a symmetric-tensor-valued function onto the interpreter's stack (as a 
// return value for a function).
void lua_pushsymtensorfunction(struct lua_State* lua, st_func_t* func);

// This helper returns true if the object at the given index is or can be 
// interpreted as a tensor-valued function, false if not.
bool lua_istensorfunction(struct lua_State* lua, int index);

// This helper retrieves a tensor-valued function from the given index on an 
// active lua interpreter, or returns NULL if the index does not point to an st_func.
st_func_t* lua_totensorfunction(struct lua_State* lua, int index);

// Pushes a tensor-valued function onto the interpreter's stack (as a 
// return value for a function).
void lua_pushtensorfunction(struct lua_State* lua, st_func_t* func);

// This helper returns true if the object at the given index is a mesh,
// false if not.
bool lua_ismesh(struct lua_State* lua, int index);

// This helper retrieves a mesh from the given index on an active lua 
// interpreter, or returns NULL if the index does not point to a mesh.
mesh_t* lua_tomesh(struct lua_State* lua, int index);

// Pushes a mesh onto the interpreter's stack (as a return value for a 
// function).
void lua_pushmesh(struct lua_State* lua, mesh_t* mesh);

// This helper returns true if the object at the given index is a user-defined
// object, false if not.
bool lua_isuserdefined(struct lua_State* lua, int index);

// This helper retrieves a user-defined type from the given index on an active lua 
// interpreter, or returns NULL if the index does not point to such a thing.
void* lua_touserdefined(struct lua_State* lua, int index);

// Pushes a user-defined object onto the interpreter's stack (as a 
// return value for a function). dtor is a destructor for the user-defined object.
void lua_pushuserdefined(struct lua_State* lua, void* userdefined, void (*dtor)(void*));

// Converts the lua table at the given index to an unordered map that associates 
// string keys with user-defined (pointer) objects. This function returns a newly-
// allocated unordered map if the table at the given index has string or integer keys,
// and NULL if no such table can be constructed.
string_ptr_unordered_map_t* lua_tounorderedmap(struct lua_State* lua, int index);

#endif

