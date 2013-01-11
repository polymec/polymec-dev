#include "cnav/cnav_ideal_gas.h"

#ifdef __cplusplus
extern "C" {
#endif

// Proton mass and kB in SI units.
static const double proton_mass = 1.67262158e-27;
static const double kB = 1.381e-23; 

typedef struct 
{
  double mu, gamma;
} cnav_ideal_gas_t;

static double ideal_gas_pressure(void* context, double rho, double T)
{
  cnav_ideal_gas_t* ideal_gas = context;
  return rho * kB * T/(proton_mass*ideal_gas->mu);
}

static double ideal_gas_temperature(void* context, double rho, double eps)
{
  cnav_ideal_gas_t* ideal_gas = context;
  return (ideal_gas->gamma - 1.0) * ideal_gas->mu * proton_mass * eps / kB;
}

cnav_eos_t* cnav_ideal_gas_new(double mu, double gamma)
{
  ASSERT(mu > 0.0);
  ASSERT(gamma > 1.0);

  char name[1024];
  sprintf(name, "Ideal gas(mu = %g AU, gamma = %g)", mu, gamma);
  cnav_ideal_gas_t* ideal_gas = malloc(sizeof(cnav_ideal_gas_t));
  ideal_gas->mu = mu;
  ideal_gas->gamma = gamma;
  cnav_eos_vtable vtable = {.pressure = ideal_gas_pressure,
                            .temperature = ideal_gas_temperature,
                            .dtor = free };
  double mass = mu * proton_mass;
  return cnav_eos_new(name, ideal_gas, vtable, 1, &mass);
}

#ifdef __cplusplus
}
#endif

