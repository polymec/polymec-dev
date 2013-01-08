#include "cnav/cnav_eos.h"

#ifdef __cplusplus
extern "C" {
#endif

struct cnav_eos_t 
{
  char* name;   // Name of the equation of state
  int num_comp; // Number of components
};

#ifdef __cplusplus
}
#endif

