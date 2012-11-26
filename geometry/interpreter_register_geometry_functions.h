#ifndef ARBI_INTERPRETER_REGISTER_GEOMETRY_FUNCTIONS_H
#define ARBI_INTERPRETER_REGISTER_GEOMETRY_FUNCTIONS_H

#include "core/interpreter.h"

#ifdef __cplusplus
extern "C" {
#endif

// This function registers basic geometry-related functions for the
// interpreter.
void interpreter_register_geometry_functions(interpreter_t* interp);

#ifdef __cplusplus
}
#endif

#endif

