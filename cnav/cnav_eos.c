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

#include <gc/gc.h>
#include "cnav/cnav_eos.h"

struct cnav_eos_t 
{
  char* name;   // Name of the equation of state
  void* context; // Context for EOS data
  int num_species; // Number of species
  double *masses; // Masses of species
  cnav_eos_vtable vtable; // Virtual table.
};

static void cnav_eos_free(void* ctx, void* dummy)
{
  cnav_eos_t* eos = ctx;
  if ((eos->context != NULL) && (eos->vtable.dtor != NULL))
    eos->vtable.dtor(eos->context);
  free(eos->name);
}

cnav_eos_t* cnav_eos_new(const char* name, void* context, 
                         cnav_eos_vtable vtable, int num_species,
                         double* masses)
{
  ASSERT(vtable.pressure != NULL);
  ASSERT(vtable.temperature != NULL);
  ASSERT(vtable.sound_speed != NULL);
  ASSERT(num_species > 0);
  ASSERT(masses != NULL);
  cnav_eos_t* eos = GC_MALLOC(sizeof(cnav_eos_t));
  eos->name = strdup(name);
  eos->context = context;
  eos->vtable = vtable;
  eos->num_species = num_species;
  eos->masses = malloc(sizeof(double)*num_species);
  for (int s = 0; s < num_species; ++s)
    eos->masses[s] = masses[s];
  GC_register_finalizer(eos, &cnav_eos_free, eos, NULL, NULL);
  return eos;
}

char* cnav_eos_name(cnav_eos_t* eos)
{
  return eos->name;
}

int cnav_eos_num_species(cnav_eos_t* eos)
{
  return eos->num_species;
}

void cnav_eos_get_masses(cnav_eos_t* eos, double* masses)
{
  for (int s = 0; s < eos->num_species; ++s)
    masses[s] = eos->masses[s];
}

double cnav_eos_temperature(cnav_eos_t* eos, double rho, double E)
{
  ASSERT(rho > 0.0);
  ASSERT(E > 0.0);
  return eos->vtable.temperature(eos->context, rho, E);
}

double cnav_eos_pressure(cnav_eos_t* eos, double rho, double T)
{
  ASSERT(rho > 0.0);
  ASSERT(T > 0.0);
  return eos->vtable.pressure(eos->context, rho, T);
}

double cnav_eos_sound_speed(cnav_eos_t* eos, double rho, double T)
{
  ASSERT(rho > 0.0);
  ASSERT(T > 0.0);
  return eos->vtable.sound_speed(eos->context, rho, T);
}

double cnav_eos_specific_internal_energy(cnav_eos_t* eos, double rho, double p)
{
  ASSERT(rho > 0.0);
  return eos->vtable.specific_internal_energy(eos->context, rho, p);
}

double cnav_eos_effective_gamma(cnav_eos_t* eos, double rho, double T)
{
  ASSERT(rho > 0.0);
  ASSERT(T > 0.0);
  return eos->vtable.effective_gamma(eos->context, rho, T);
}

