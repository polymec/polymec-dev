#ifndef ARBI_INTERPRETER_H
#define ARBI_INTERPRETER_H

#include "core/arbi.h"
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
  INTERPRETER_FUNCTION,
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

// Fetches the given number from the interpreter, returning -FLT_MAX if it 
// is not found or if it is not a number.
double interpreter_get_number(interpreter_t* interp, const char* name);

// Fetches the given mesh from the interpreter, returning NULL if it 
// is not found or if it is not a table. The caller assumes responsibility
// for destroying the mesh after this call.
mesh_t* interpreter_get_mesh(interpreter_t* interp, const char* name);

// Fetches the given space-time function from the interpreter, returning NULL 
// if it is not found or if it is not a table. 
st_func_t* interpreter_get_function(interpreter_t* interp, const char* name);

// Fetches the given table from the interpreter, returning NULL if it 
// is not found or if it is not a table. The caller assumes responsibility
// for destroying the table after this call.
str_ptr_unordered_map_t* interpreter_get_table(interpreter_t* interp, const char* name);

// Fetches the given user-defined object from the interpreter, returning NULL 
// if it is not found or if it is not a user-defined object. The caller 
// assumes responsibility for destroying the user-defined object after 
// this call.
void* interpreter_get_user_defined(interpreter_t* interp, const char* name);

//------------------------------------------------------------------------
// These methods can be used to create functions that extend an 
// interpreter.
//------------------------------------------------------------------------

// This helper returns true if the object at the given index is a space-time
// function, false if not.
bool lua_isstfunc(struct lua_State* lua, int index);

// This helper retrieves an st_func from the given index on an active lua 
// interpreter, or returns NULL if the index does not point to an st_func.
st_func_t* lua_tostfunc(struct lua_State* lua, int index);

// Pushes a space-time function onto the interpreter's stack (as a 
// return value for a function).
void lua_pushstfunc(struct lua_State* lua, st_func_t* func);

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

