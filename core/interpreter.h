#ifndef ARBI_INTERPRETER_H
#define ARBI_INTERPRETER_H

#include "core/arbi.h"
#include "core/unordered_map.h"

#ifdef __cplusplus
extern "C" {
#endif

// This class interprets input files and stores values of variables therein.
typedef struct interpreter_t interpreter_t;

// This defines all valid input variable types.
typedef enum
{
  INTERPRETER_STRING,
  INTERPRETER_NUMBER,
  INTERPRETER_MESH,
  INTERPRETER_FUNCTION,
  INTERPRETER_TABLE
} interpreter_var_type_t;

// This is used to validate variables found within input files. It is a 
// variable/type pair that associates variables names with their intended
// types.
typedef struct 
{
  char* variable;               // Name of the variable.
  interpreter_var_type_t type;  // Type of the variable.
} interpreter_validation_t;

// Creates a new interpreter for use with a model. The interpreter uses the 
// given set of types to validate variables in input files.
interpreter_t* interpreter_new(interpreter_validation_t* valid_inputs, int num_valid_inputs);

// Destroys the given interpreter.
void interpreter_free(interpreter_t* interp);

// Parses the input string, storing values in the interpreter.
void interpreter_parse(interpreter_t* interp, char* input);

// Fetches the given string from the interpreter, returning NULL if it 
// is not found or if it is not a string. The string refers to storage 
// within the interpreter and should not be freed.
char* interpreter_get_string(interpreter_t* interp, const char* name);

// Fetches the given number from the interpreter, returning -FLT_MAX if it 
// is not found or if it is not a number.
double interpreter_get_number(interpreter_t* interp, const char* name);

// Fetches the given table from the interpreter, returning NULL if it 
// is not found or if it is not a table. The caller assumes responsibility
// for destroying the table after this call.
str_ptr_unordered_map_t* interpreter_get_table(interpreter_t* interp, const char* name);

#ifdef __cplusplus
}
#endif

#endif

