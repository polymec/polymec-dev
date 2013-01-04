#include "cnav/cnav_model.h"

#ifdef __cplusplus
extern "C" {
#endif

int main(int argc, char** argv)
{
  // Just use the generic model driver for the cnav model.
  return model_main("cnav", &cnav_model_new, argc, argv);
}

#ifdef __cplusplus
}
#endif
