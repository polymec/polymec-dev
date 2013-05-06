#ifndef POLYMEC_INTERPRETER_H
#define POLYMEC_INTERPRETER_H

#include "core/polymec.h"
#include "core/mesh.h"
#include "core/st_func.h"
#include "core/unordered_map.h"

#ifdef __cplusplus
extern "C" {
#endif

// Forward declaration for Lua innards.
struct lua_State;

// This class interprets input files and stores values of variables therein.
typedef struct interpreter_t interpreter_t;

// This defines all valid input variable types.
typedef enum
{
  INTERPRETER_STRING,
  INTERPRETER_NUMBER,
  INTERPRETER_MESH,
  INTERPRETER_POINT,
  INTERPRETER_POINT_LIST,
  INTERPRETER_VECTOR,
  INTERPRETER_VECTOR_LIST,
  INTERPRETER_BOUNDING_BOX,
  INTERPRETER_SCALAR_FUNCTION,
  INTERPRETER_VECTOR_FUNCTION,
  INTERPRETER_SEQUENCE,
  INTERPRETER_TABLE,
  INTERPRETER_USER_DEFINED, // <-- on your head be it!!
  INTERPRETER_TERMINUS // Used only to terminate validation lists.
} interpreter_var_type_t;

// This is used to validate variables found within input files. It is a 
// variable/type pair that associates variables names with their intended
// types.
typedef struct 
{
  char* variable;               // Name of the variable.
  interpreter_var_type_t type;  // Type of the variable.
} interpreter_validation_t;

// This validation object is used to terminate lists of valid inputs.
static const interpreter_validation_t END_OF_VALID_INPUTS = {.variable = "TERMINUS", .type = INTERPRETER_TERMINUS};

// Creates a new interpreter for use with a model. The interpreter uses the 
// given set of types to validate variables in input files.
interpreter_t* interpreter_new(interpreter_validation_t* valid_inputs);

// Destroys the given interpreter.
void interpreter_free(interpreter_t* interp);

// Register a user-defined function with the interpreter.
void interpreter_register_function(interpreter_t* interp, const char* function_name, int (*function)(struct lua_State*));

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
double interpreter_get_number(interpreter_t* interp, const char* name);

// Sets the given variable within the interpreter to the given numeric 
// value. Any existing value of this variable is overwritten.
void interpreter_set_number(interpreter_t* interp, const char* name, double value);

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
bbox_t* interpreter_get_boundingbox(interpreter_t* interp, const char* name);

// Sets the given variable within the interpreter to the given bounding box.
// Any existing value of this variable is overwritten.
void interpreter_set_boundingbox(interpreter_t* interp, const char* name, bbox_t* value);

// Fetches the given mesh from the interpreter, returning NULL if it 
// is not found or if it is not a table. The caller assumes responsibility
// for destroying the mesh after this call.
mesh_t* interpreter_get_mesh(interpreter_t* interp, const char* name);

// Sets the given variable within the interpreter to the given mesh 
// object. Any existing value of this variable is overwritten.
// Control of the mesh is assumed by the interpreter.
void interpreter_set_mesh(interpreter_t* interp, const char* name, mesh_t* value);

// Fetches the given scalar function from the interpreter, returning NULL 
// if it is not found or if it is not a table. 
st_func_t* interpreter_get_scalar_function(interpreter_t* interp, const char* name);

// Sets the given variable within the interpreter to the given scalar function 
// object. Any existing value of this variable is overwritten.
void interpreter_set_scalar_function(interpreter_t* interp, const char* name, st_func_t* value);

// Fetches the given vector function from the interpreter, returning NULL 
// if it is not found or if it is not a table. 
st_func_t* interpreter_get_vector_function(interpreter_t* interp, const char* name);

// Sets the given variable within the interpreter to the given vector function 
// object. Any existing value of this variable is overwritten.
void interpreter_set_vector_function(interpreter_t* interp, const char* name, st_func_t* value);

// Fetches the given table from the interpreter, returning NULL if it 
// is not found or if it is not a table. The caller assumes responsibility
// for destroying the table after this call.
str_ptr_unordered_map_t* interpreter_get_table(interpreter_t* interp, const char* name);

// Sets the given variable within the interpreter to the given sequence of
// nubmers.. Any existing value of this variable is overwritten.
void interpreter_set_sequence(interpreter_t* interp, const char* name, double* sequence, int len);

// Sets the given variable within the interpreter to the given table of
// objects. Any existing value of this variable is overwritten.
void interpreter_set_table(interpreter_t* interp, const char* name, str_ptr_unordered_map_t* value);

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
double* lua_tosequence(struct lua_State* lua, int index, int* len);

// Pushes a sequence onto the interpreter's stack (as a return value for a 
// function).
void lua_pushsequence(struct lua_State* lua, double* sequence, int len);

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

// This helper returns true if the object at the given index is a scalar-valued
// function, false if not.
bool lua_isscalarfunction(struct lua_State* lua, int index);

// This helper retrieves a scalar-valued function from the given index on an 
// active lua interpreter, or returns NULL if the index does not point to an st_func.
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

#ifdef __cplusplus
}
#endif

#endif

