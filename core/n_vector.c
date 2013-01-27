#include "core/n_vector.h"

struct n_vector_t
{
  char* impl_name;
  void* context;
  n_vector_vtable vtable;
};

n_vector_t* n_vector_new(const char* impl_name,
                     void* context,
                     n_vector_vtable vtable)
{
  ASSERT(context != NULL);
  ASSERT(vtable.clone != NULL);
  ASSERT(vtable.data != NULL);
  ASSERT(vtable.local_size != NULL);
  ASSERT(vtable.global_size != NULL);
  n_vector_t* vec = malloc(sizeof(n_vector_t));
  vec->impl_name = strdup(impl_name);
  vec->context = context;
  vec->vtable = vtable;
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

n_vector_t* n_vector_clone(n_vector_t* vector)
{
  n_vector_t* vec = malloc(sizeof(n_vector_t));
  vec->impl_name = strdup(vector->impl_name);
  vec->context = vector->vtable.clone(vector->context);
  vec->vtable = vector->vtable;
  return vec;
}

double* n_vector_data(n_vector_t* vector)
{
  return vector->vtable.data(vector->context);
}

int n_vector_local_size(n_vector_t* vector)
{
  return vector->vtable.local_size(vector->context);
}

int n_vector_global_size(n_vector_t* vector)
{
  return vector->vtable.global_size(vector->context);
}

void n_vector_free(n_vector_t* vector)
{
  free(vector->impl_name);
  if (vector->vtable.dtor)
    vector->vtable.dtor(vector->context);
  free(vector);
}

#ifdef __cplusplus
}
#endif

