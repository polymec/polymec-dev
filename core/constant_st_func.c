#include "core/constant_st_func.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct
{
  int num_comp;
  double *comp;
} const_st_func_t;

static constant_eval(void* ctx, point_t* x, double t, double* res)
{
  const_st_func_t* f = (const_st_func_t*)ctx;
  for (int i = 0; i < f->num_comp; ++i)
    res[i] = f->comp[i];
}

static constant_dtor(void* ctx)
{
  const_st_func_t* f = (const_st_func_t*)ctx;
  free(f->comp);
  free(f);
}

st_func_t* constant_st_func_new(int num_comp, double comp[])
{
  st_vtable vtable = {.eval = &constant_eval, .dtor = &constant_dtor};
  char name[1024];
  snprintf(name, 1024, "constant space-time function"); // FIXME
  const_st_func_t* f = malloc(sizeof(const_st_func_t));
  f->num_comp = num_comp;
  f->comp = malloc(num_comp*sizeof(double));
  for (int i = 0; i < num_comp; ++i)
    f->comp[i] = comp[i];
  return st_func_new(name, (void*)f, vtable, ST_HOMOGENEOUS, ST_CONSTANT, num_comp);
}

#ifdef __cplusplus
}
#endif

