#include "advect/advect_model.h"

int main(int argc, char** argv)
{
  // Just use the generic model driver for the Advect model.
  return model_main("advect", &advect_model_new, argc, argv);
}

