#include "core/n_vector.h"

struct n_vector_t
{
  char* impl_name;
  void* context;
  n_vector_dtor dtor;
};

n_vector_t* n_vector_new(const char* impl_name,
                     void* context,
                     n_vector_dtor dtor)
{
  ASSERT(context != NULL);
  n_vector_t* vec = malloc(sizeof(n_vector_t));
  vec->impl_name = strdup(impl_name);
  vec->context = context;
  vec->dtor = dtor;
  return vec;
}

char* n_vector_impl_name(n_vector_t* vector)
{
  return vector->impl_name;
}

void* n_vector_context(n_vector_t* vector)
{
  return vector->context;
}

void n_vector_free(n_vector_t* vector)
{
  free(vector->impl_name);
  if (vector->dtor)
    vector->dtor(vector->context);
  free(vector);
}

#ifdef __cplusplus
}
#endif

