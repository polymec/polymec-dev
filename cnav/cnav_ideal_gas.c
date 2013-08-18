// Copyright 2012-2013 Jeffrey Johnson.
// 
// This file is part of Polymec, and is licensed under the Apache License, 
// Version 2.0 (the "License"); you may not use this file except in 
// compliance with the License. You may may find the text of the license in 
// the LICENSE file at the top-level source directory, or obtain a copy of 
// it at
// 
// http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "cnav/cnav_ideal_gas.h"

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

static double ideal_gas_sound_speed(void* context, double rho, double eps)
{
  cnav_ideal_gas_t* ideal_gas = context;
  double gamma = ideal_gas->gamma;
  return sqrt(gamma * pow(rho, gamma - 1.0));
}

static double ideal_gas_specific_internal_energy(void* context, double rho, double p)
{
  cnav_ideal_gas_t* ideal_gas = context;
  return p / (rho * (ideal_gas->gamma - 1.0));
}

static double ideal_gas_effective_gamma(void* context, double rho, double T)
{
  cnav_ideal_gas_t* ideal_gas = context;
  return ideal_gas->gamma;
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
                            .sound_speed = ideal_gas_sound_speed,
                            .specific_internal_energy = ideal_gas_specific_internal_energy,
                            .effective_gamma = ideal_gas_effective_gamma,
                            .dtor = free };
  double mass = mu * proton_mass;
  return cnav_eos_new(name, ideal_gas, vtable, 1, &mass);
}

