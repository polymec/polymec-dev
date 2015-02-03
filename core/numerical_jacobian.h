#ifndef POLYMEC_NUMERICAL_JACOBIAN_H
#define POLYMEC_NUMERICAL_JACOBIAN_H

#include "core/polymec.h"
// compute_numerical_jacobian : a function for assisting in debugging
// 
// computes the numerical jacobian of a rhs function for a serial run
// requires the standard inputs for the rhs function, as well as
// 
// n (int)              size of x array == rhs array for serial process
// dx (real_t)          increment in x to use for the jacobian
// jacobian (real_t**)  the 2d jacobian matrix which must be sized prior to call 
static inline int compute_numerical_jacobian(int (*rhs)(void* context, real_t t, real_t* x, real_t* xdot), 
					     void* context,
					     real_t t, 
					     real_t* xin, 
					     int n, 
					     real_t dx,
					     real_t** jacobian){
  //ASSERT SERIAL PROCESS 
  int status;
  real_t dx_1 = 1./dx; //(one over) the step size for the numerical difference
  real_t xold;
  real_t rhs1[n]; //rhs at x + dx
  real_t rhs0[n]; //rhs at x
  
  status = rhs(context, t, xin, &rhs0[0]);
  if (status != 0)
    {
      log_debug("numerical_jacobian: call to RHS at initial x, t failed.");
      return status;
    }

  
  for(int j = 0; j < n; ++j){
    //store initial value for this dof
    xold = xin[j];
    
    //increment x and evaluate new rhs value
    xin[j] += dx;
    status = rhs(context, t, xin, &rhs1[0]);
    if (status != 0)
      {
	log_debug("numerical_jacobian: call to RHS at x + dx[%d], t failed.", j);
	break; 
      }
    
    for(int i = 0; i < n; ++i){
      (jacobian)[i][j] = (rhs1[i] - rhs0[i]) * dx_1;
    }
    //restore array to previous state
    xin[j] = xold;
  }
  return status;
}

  


  

#endif
