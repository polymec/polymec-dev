#ifndef ARBI_INTERPRETER_H
#define ARBI_INTERPRETER_H

#include "core/arbi.h"

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
  INTERPRETER_SP_FUNC,
  INTERPRETER_ST_FUNC,
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

// Fetches the given named value from the interpreter, returning NULL if it 
// is not found.
void* interpreter_get(interpreter_t* interp, const char* name);

#ifdef __cplusplus
}
#endif

#endif

