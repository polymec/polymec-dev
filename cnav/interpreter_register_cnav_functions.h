#ifndef POLYMEC_INTERPRETER_REGISTER_CNAV_FUNCTIONS_H
#define POLYMEC_INTERPRETER_REGISTER_CNAV_FUNCTIONS_H

#include "core/interpreter.h"
#include "cnav/cnav_bc.h"

// Equips the given interpreter with functions specific to the cnav model.
void interpreter_register_cnav_functions(interpreter_t* interpreter);

#endif
