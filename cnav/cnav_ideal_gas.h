#ifndef POLYMEC_CNAV_IDEAL_GAS_H
#define POLYMEC_CNAV_IDEAL_GAS_H

#include "cnav/cnav_eos.h"

// This constructs an equation of state representing a single ideal gas 
// with the given molecular mass mu (in AU) and polytropic index gamma.
cnav_eos_t* cnav_ideal_gas_new(double mu, double gamma);

#endif

