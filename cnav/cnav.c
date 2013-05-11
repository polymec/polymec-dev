#include "cnav/cnav_model.h"

int main(int argc, char** argv)
{
  // Just use the generic model driver for the cnav model.
  return model_main("cnav", &cnav_model_new, argc, argv);
}

