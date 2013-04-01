#include <gc/gc.h>
#include "geometry/mapping.h"

#ifdef __cplusplus
extern "C" {
#endif

struct mapping_t 
{
  char* name;
  void* context;
  mapping_vtable vtable;
};

static void mapping_free(void* ctx, void* dummy)
{
  mapping_t* mapping = (mapping_t*)ctx;
  if (mapping->vtable.dtor)
    free(mapping->context);
  free(mapping->name);
}

mapping_t* mapping_new(const char* name, void* context, mapping_vtable vtable)
{
  ASSERT(context != NULL);
  ASSERT(vtable.map != NULL);
  ASSERT(vtable.jacobian != NULL);
  mapping_t* m = GC_MALLOC(sizeof(mapping_t));
  m->name = strdup(name);
  m->context = context;
  m->vtable = vtable;
  GC_register_finalizer(m, mapping_free, m, NULL, NULL);
  return m;
}

const char* mapping_name(mapping_t* mapping)
{
  return mapping->name;
}

void* mapping_context(mapping_t* mapping)
{
  return mapping->context;
}

void mapping_map(mapping_t* mapping, point_t* x, point_t* y)
{
  mapping->vtable.map(mapping->context, x, y);
}

void mapping_compute_jacobian(mapping_t* mapping, point_t* x, double* J)
{
  mapping->vtable.jacobian(mapping->context, x, J);
}

#ifdef __cplusplus
}
#endif

