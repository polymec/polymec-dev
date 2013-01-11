#ifndef POLYMEC_CNAV_IDEAL_GAS_H
#define POLYMEC_CNAV_IDEAL_GAS_H

#include "cnav/cnav_eos.h"

#ifdef __cplusplus
extern "C" {
#endif

// This constructs an equation of state representing a single ideal gas 
// with the given molecular mass mu (in AU) and polytropic index gamma.
cnav_eos_t* cnav_ideal_gas_new(double mu, double gamma);

#ifdef __cplusplus
}
#endif

#endif

