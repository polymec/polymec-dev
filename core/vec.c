#include "core/vec.h"

#ifdef __cplusplus
extern "C" {
#endif

// We package a context object with a destructor to represent a vector.
struct vec_t 
{
  void* context;
  void (*dtor)(void*); // Destructor.
};

#ifdef __cplusplus
}
#endif

