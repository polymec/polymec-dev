// Copyright (c) 2012-2014, Jeffrey N. Johnson
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without 
// modification, are permitted provided that the following conditions are met:
// 
// 1. Redistributions of source code must retain the above copyright notice, this 
// list of conditions and the following disclaimer.
// 
// 2. Redistributions in binary form must reproduce the above copyright notice, 
// this list of conditions and the following disclaimer in the documentation 
// and/or other materials provided with the distribution.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE 
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <string.h>

#include "cmockery.h"
#include "core/polymec.h"
#include "integrators/nonlinear_integrator.h"

// This test solves a nonlinear system that arises from a system
// of partial differential equations. The PDE system is a food web
// population model, with predator-prey interaction and diffusion
// on the unit square in two dimensions. The dependent variable
// vector is the following:
// 
//       1   2         ns
// c = (c , c ,  ..., c  )     (denoted by the variable cc)
// 
// and the PDE's are as follows:
//
//                    i       i
//         0 = d(i)*(c     + c    )  +  f  (x,y,c)   (i=1,...,ns)
//                    xx      yy         i
//
//   where
//
//                   i             ns         j
//   f  (x,y,c)  =  c  * (b(i)  + sum a(i,j)*c )
//    i                           j=1
//
// The number of species is ns = 2 * np, with the first np being
// prey and the last np being predators. The number np is both the
// number of prey and predator species. The coefficients a(i,j),
// b(i), d(i) are:
//
//   a(i,i) = -AA   (all i)
//   a(i,j) = -GG   (i <= np , j >  np)
//   a(i,j) =  EE   (i >  np,  j <= np)
//   b(i) = BB * (1 + alpha * x * y)   (i <= np)
//   b(i) =-BB * (1 + alpha * x * y)   (i >  np)
//   d(i) = DPREY   (i <= np)
//   d(i) = DPRED   ( i > np)
//
// The various scalar parameters are set using define's or in
// routine InitUserData.
//
// The boundary conditions are: normal derivative = 0, and the
// initial guess is constant in x and y, but the final solution
// is not.
//
// The PDEs are discretized by central differencing on an MX by
// MY mesh.
// 
// The nonlinear system is solved by KINSOL using the method
// specified in local variable globalstrat.
//
// The preconditioner matrix is a block-diagonal matrix based on
// the partial derivatives of the interaction terms f only.
//
// Constraints are imposed to make all components of the solution
// positive.
// -----------------------------------------------------------------
// References:
//
// 1. Peter N. Brown and Youcef Saad,
//    Hybrid Krylov Methods for Nonlinear Systems of Equations
//    LLNL report UCRL-97645, November 1987.
//
// 2. Peter N. Brown and Alan C. Hindmarsh,
//    Reduced Storage Matrix Methods in Stiff ODE systems,
//    Lawrence Livermore National Laboratory Report  UCRL-95088,
//    Rev. 1, June 1987, and  Journal of Applied Mathematics and
//    Computation, Vol. 31 (May 1989), pp. 40-91. (Presents a
//    description of the time-dependent version of this test
//    problem.)
// -----------------------------------------------------------------

typedef struct 
{
  realtype **P[MX][MY];
  int *pivot[MX][MY];
  realtype **acoef, *bcoef;
  N_Vector rates;
  realtype *cox, *coy;
  realtype ax, ay, dx, dy;
  realtype uround, sqruround;
  long int mx, my, ns, np;

} user_data_t;

// Interaction rate function routine 
static void web_rate(void* context, real_t xx, real_t yy, real_t *cxy, real_t *ratesxy)
{
  long int i;
  real_t fac;
  user_data_t* data = context;
  
  for (i = 0; i<NUM_SPECIES; i++)
    ratesxy[i] = DotProd(NUM_SPECIES, cxy, acoef[i]);
  
  fac = ONE + ALPHA * xx * yy;
  
  for (i = 0; i < NUM_SPECIES; i++)
    ratesxy[i] = cxy[i] * ( bcoef[i] * fac + ratesxy[i] );  
}

static int foodweb_func(void* context, real_t t, real_t* cc, real_t* fval)
{
  real_t xx, yy, delx, dely, *cxy, *rxy, *fxy, dcyli, dcyui, dcxli, dcxri;
  long int jx, jy, is, idyu, idyl, idxr, idxl;
  user_data_t* data = context;
  
  delx = data->dx;
  dely = data->dy;
  
  // Loop over all mesh points, evaluating rate array at each point
  for (jy = 0; jy < MY; jy++) 
  {
    
    yy = dely*jy;

    // Set lower/upper index shifts, special at boundaries. 
    idyl = (jy != 0   ) ? NSMX : -NSMX;
    idyu = (jy != MY-1) ? NSMX : -NSMX;
    
    for (jx = 0; jx < MX; jx++) {

      xx = delx*jx;

      // Set left/right index shifts, special at boundaries. 
      idxl = (jx !=  0  ) ?  NUM_SPECIES : -NUM_SPECIES;
      idxr = (jx != MX-1) ?  NUM_SPECIES : -NUM_SPECIES;

      cxy = IJ_Vptr(cc,jx,jy);
      rxy = IJ_Vptr(data->rates,jx,jy);
      fxy = IJ_Vptr(fval,jx,jy);

      // Get species interaction rate array at (xx,yy) 
      web_rate(user_data, xx, yy, cxy, rxy);
      
      for(is = 0; is < NUM_SPECIES; is++) {
        
        // Differencing in x direction 
        dcyli = *(cxy+is) - *(cxy - idyl + is) ;
        dcyui = *(cxy + idyu + is) - *(cxy+is);
        
        // Differencing in y direction 
        dcxli = *(cxy+is) - *(cxy - idxl + is);
        dcxri = *(cxy + idxr +is) - *(cxy+is);
        
        // Compute the total rate value at (xx,yy) 
        fxy[is] = (coy)[is] * (dcyui - dcyli) +
          (cox)[is] * (dcxri - dcxli) + rxy[is];
      }
    }
  } 

  return(0);
}

void test_foodweb_integrator(void** state)
{
  user_data_t* user_data = user_data_new();
  nonlinear_integrator_vtable vtable; // FIXME
  nonlinear_integrator_t* integ = nonlinear_integrator_new("Food web",
                                                           user_data,
                                                           MPI_COMM_SELF,
                                                           vtable, GMRES, 5);
  assert_true(strcmp(nonlinear_integrator_name(integ), "Food web") == 0);
  assert_true((user_data_t*)nonlinear_integrator_context(integ) == user_data);

  real_t X[NUM_SPECIES];
  int num_iters;
  bool solved = nonlinear_integrator_solve(integ, 0.0, X, &num_iters);
  assert_true(solved);
  nonlinear_integrator_free(integ);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const UnitTest tests[] = 
  {
    unit_test(test_foodweb_integrator),
  };
  return run_tests(tests);
}
