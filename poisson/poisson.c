#include "poisson/poisson_model.h"

#ifdef __cplusplus
extern "C" {
#endif

int main(int argc, char** argv)
{
  // Just use the generic model driver for the Poisson model.
  return model_main("poisson", &poisson_model_new, argc, argv);
}

#ifdef __cplusplus
}
#endif
