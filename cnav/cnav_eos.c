#include "cnav/cnav_eos.h"

#ifdef __cplusplus
extern "C" {
#endif

struct cnav_eos_t 
{
  char* name;   // Name of the equation of state
  int num_comp; // Number of components
};

char* cnav_eos_name(cnav_eos_t* eos)
{
  return eos->name;
}

int cnav_eos_num_species(cnav_eos_t* eos)
{
  return eos->num_comp;
}

#ifdef __cplusplus
}
#endif

