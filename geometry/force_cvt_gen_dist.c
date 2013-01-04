#include "core/slist.h"
#include "core/point_set.h"
#include "geometry/force_cvt_gen_dist.h"

#ifdef __cplusplus
extern "C" {
#endif

struct cvt_gen_force_t 
{
};

cvt_gen_force_t* linear_spring_force_new(double k)
{
  return NULL;
}

cvt_gen_dist_t* force_cvt_gen_dist_new(cvt_gen_force_t* force,
                                       double tolerance)
{
  return NULL;
}

#ifdef __cplusplus
}
#endif

