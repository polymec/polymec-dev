#ifndef POLYMEC_INTERPRETER_REGISTER_ADVECT_FUNCTIONS_H
#define POLYMEC_INTERPRETER_REGISTER_ADVECT_FUNCTIONS_H

#include "core/interpreter.h"
#include "advect/advect_bc.h"

// Equips the given interpreter with functions specific to the advect model.
void interpreter_register_advect_functions(interpreter_t* interpreter);

#endif
