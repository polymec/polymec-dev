#ifndef POLYMEC_CNAV_EOS_H
#define POLYMEC_CNAV_EOS_H

#include "core/model.h"

#ifdef __cplusplus
extern "C" {
#endif

// This type represents a multi-component equation of state for the 
// compressible Navier-Stokes solver. Objects of this type are 
// garbage-collected.
typedef struct cnav_eos_t cnav_eos_t;

// Returns the name of the equation of state.
char* cnav_eos_name(cnav_eos_t* eos);

// Returns the number of species in the material described by this
// equation of state.
int cnav_eos_num_species(cnav_eos_t* eos);

#ifdef __cplusplus
}
#endif

#endif

