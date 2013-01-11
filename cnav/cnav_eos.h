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

// A function pointer type for evaluating the pressure.
typedef double (*cnav_eos_pressure_func)(void*, double, double);

// A function pointer type for evaluating the temperature.
typedef double (*cnav_eos_temperature_func)(void*, double, double);

// A destructor for any given context object.
typedef void (*cnav_eos_dtor)(void*);

// This virtual table must be implemented by any equation of state.
typedef struct 
{
  cnav_eos_pressure_func    pressure;
  cnav_eos_temperature_func temperature;
  cnav_eos_dtor             dtor;
} cnav_eos_vtable;

// Constructs a new equation of state object from the given parameters.
cnav_eos_t* cnav_eos_new(const char* name, void* context, 
                         cnav_eos_vtable vtable, int num_species,
                         double* masses);

// Returns the name of the equation of state.
char* cnav_eos_name(cnav_eos_t* eos);

// Returns the number of species in the material described by this
// equation of state.
int cnav_eos_num_species(cnav_eos_t* eos);

// Retrieves the masses for the species in this equation of state.
void cnav_eos_get_masses(cnav_eos_t* eos, double* masses);

// Returns the temperature given mass density and specific internal energy.
double cnav_eos_temperature(cnav_eos_t* eos, double rho, double eps);

// Returns the pressure given mass density and temperature.
double cnav_eos_pressure(cnav_eos_t* eos, double rho, double T);

#ifdef __cplusplus
}
#endif

#endif

