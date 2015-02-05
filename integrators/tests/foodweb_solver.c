// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/polymec.h"
#include "integrators/newton_solver.h"
#include "integrators/cpr_pc.h"

// We use this for some of the underlying data structures.
#include "sundials/sundials_direct.h"

//------------------------------------------------------------------------
//                          FOOD WEB PROBLEM
//------------------------------------------------------------------------
// This was adapted from KINSol's food web example problem, programmed by 
// Allan Taylor, Alan Hindmarsh, and Radu Serban @ LLNL.
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

/* Problem Constants */

#define NUM_SPECIES     6  /* must equal 2*(number of prey or predators)
                              number of prey = number of predators       */ 

#define PI       RCONST(3.1415926535898)   /* pi */ 

#define MX          8              // MX = number of x mesh points */
#define MY          8              // MY = number of y mesh points */
#define NSMX        (NUM_SPECIES * MX)
#define NEQ         (NSMX * MY)    // number of equations in the system 
#define AA          RCONST(1.0)    // value of coefficient AA in above eqns 
#define EE          RCONST(10000.) // value of coefficient EE in above eqns 
#define GG          RCONST(0.5e-6) // value of coefficient GG in above eqns 
#define BB          RCONST(1.0)    // value of coefficient BB in above eqns 
#define DPREY       RCONST(1.0)    // value of coefficient dprey above 
#define DPRED       RCONST(0.5)    // value of coefficient dpred above 
#define ALPHA       RCONST(1.0)    // value of coefficient alpha above 
#define AX          RCONST(1.0)    // total range of x variable 
#define AY          RCONST(1.0)    // total range of y variable 
#define FTOL        RCONST(1.e-7)  // ftol tolerance 
#define STOL        RCONST(1.e-13) // stol tolerance 
#define THOUSAND    RCONST(1000.0) // one thousand 
#define ZERO        RCONST(0.)     // 0. 
#define ONE         RCONST(1.0)    // 1. 
#define TWO         RCONST(2.0)    // 2. 
#define PREYIN      RCONST(1.0)    // initial guess for prey concentrations. 
#define PREDIN      RCONST(30000.0)// initial guess for predator concs.  

// User-defined vector access macro: IJ_Vptr 

// IJ_Vptr is defined in order to translate from the underlying 3D structure
// of the dependent variable vector to the 1D storage scheme for an N-vector.
// IJ_Vptr(vv,i,j) returns a pointer to the location in vv corresponding to 
// indices is = 0, jx = i, jy = j.    

#define IJ_Vptr(vv,i,j)   (&vv[i*NUM_SPECIES + j*NSMX])

typedef struct 
{
  real_t **acoef, *bcoef;
  real_t *rates;
  real_t *cox, *coy;
  real_t ax, ay, dx, dy;
  real_t uround, sqruround;
  long int mx, my, ns, np;

  // Sparsity graph for preconditioner matrix.
  adj_graph_t* sparsity;

} foodweb_t;

// Readability definitions used in other routines below.
#define acoef  (data->acoef)
#define bcoef  (data->bcoef)
#define cox    (data->cox)
#define coy    (data->coy)

// Newly initialized food web data context.
static foodweb_t* foodweb_new()
{
  real_t *a1,*a2, *a3, *a4;

  foodweb_t* data = malloc(sizeof(foodweb_t));
  acoef = newDenseMat(NUM_SPECIES, NUM_SPECIES);
  bcoef = malloc(NUM_SPECIES * sizeof(real_t));
  cox   = malloc(NUM_SPECIES * sizeof(real_t));
  coy   = malloc(NUM_SPECIES * sizeof(real_t));
  
  data->mx = MX;
  data->my = MY;
  data->ns = NUM_SPECIES;
  data->np = NUM_SPECIES/2;
  data->ax = AX;
  data->ay = AY;
  data->dx = (data->ax)/(MX-1);
  data->dy = (data->ay)/(MY-1);
  data->uround = UNIT_ROUNDOFF;
  data->sqruround = sqrt(data->uround);
  data->rates = malloc(sizeof(real_t) * NEQ);

  // Set up the coefficients a and b plus others found in the equations.
  long np = data->np;

  real_t dx2=(data->dx)*(data->dx), dy2=(data->dy)*(data->dy);

  for (long i = 0; i < np; i++) 
  {
    a1= &(acoef[i][np]);
    a2= &(acoef[i+np][0]);
    a3= &(acoef[i][0]);
    a4= &(acoef[i+np][np]);

    // Fill in the portion of acoef in the four quadrants, row by row...
    for (long j = 0; j < np; j++) 
    {
      *a1++ =  -GG;
      *a2++ =   EE;
      *a3++ = ZERO;
      *a4++ = ZERO;
    }

    // ...and then change the diagonal elements of acoef to -AA.
    acoef[i][i]=-AA;
    acoef[i+np][i+np] = -AA;

    bcoef[i] = BB;
    bcoef[i+np] = -BB;

    cox[i]=DPREY/dx2;
    cox[i+np]=DPRED/dx2;

    coy[i]=DPREY/dy2;
    coy[i+np]=DPRED/dy2;
  }  

  // Construct a sparsity graph.
  adj_graph_t* sparsity = adj_graph_new(MPI_COMM_SELF, MX*MY);
  for (int jy = 0; jy < MY; jy++) 
  {
    // Set lower/upper index shifts, special at boundaries. 
    int idyl = (jy != 0   ) ? MX : -MX;
    int idyu = (jy != MY-1) ? MX : -MX;
    
    for (int jx = 0; jx < MX; jx++) 
    {
      // Set left/right index shifts, special at boundaries. 
      int idxl = (jx !=  0  ) ?  1 : -1;
      int idxr = (jx != MX-1) ?  1 : -1;

      // Find the edges for the vertex corresponding to (jx, jy).
      int idx = jx + MX*jy;
      int num_edges = 0;
      int edges[4];
      if (jy > 0)
        edges[num_edges++] = idx - idyl; // lower
      if (jy < (MY-1))
        edges[num_edges++] = idx + idyu; // upper
      if (jx > 0)
        edges[num_edges++] = idx - idxl; // left
      if (jx < (MX-1))
        edges[num_edges++] = idx + idxr; // right

      // Set the edges within the sparsity graph.
      adj_graph_set_num_edges(sparsity, idx, num_edges);
      memcpy(adj_graph_edges(sparsity, idx), edges, sizeof(int) * num_edges);
    }
  }

  data->sparsity = adj_graph_new_with_block_size(sparsity, NUM_SPECIES);
  adj_graph_free(sparsity);

  return data;
}

// Food web data context destructor.
static void foodweb_dtor(void* context)
{
  foodweb_t* data = context;
  
  destroyMat(acoef);
  free(bcoef);
  free(cox);
  free(coy);
  free(data->rates);
  adj_graph_free(data->sparsity);
  free(data);
}

// Dot product routine for real_t arrays 
static real_t dot_prod(long int size, real_t *x1, real_t *x2)
{
  long int i;
  real_t *xx1, *xx2, temp = ZERO;
  
  xx1 = x1; xx2 = x2;
  for (i = 0; i < size; i++) 
    temp += (*xx1++) * (*xx2++);

  return temp;
}

// Interaction rate function routine 
static void web_rate(void* context, real_t xx, real_t yy, real_t *cxy, real_t *ratesxy)
{
  long int i;
  real_t fac;
  foodweb_t* data = context;
  
  for (i = 0; i<NUM_SPECIES; i++)
    ratesxy[i] = dot_prod(NUM_SPECIES, cxy, acoef[i]);
  
  fac = ONE + ALPHA * xx * yy;
  
  for (i = 0; i < NUM_SPECIES; i++)
    ratesxy[i] = cxy[i] * ( bcoef[i] * fac + ratesxy[i] );  
}

static int foodweb_func(void* context, real_t t, real_t* cc, real_t* fval)
{
  real_t xx, yy, delx, dely, *cxy, *rxy, *fxy, dcyli, dcyui, dcxli, dcxri;
  long int jx, jy, is, idyu, idyl, idxr, idxl;
  foodweb_t* data = context;
  
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
      web_rate(data, xx, yy, cxy, rxy);
      
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

  return 0;
}

static void foodweb_set_x_scale(void* context, real_t* x_scale)
{
  int i, jx, jy;
  real_t *sloc;
  real_t stemp[NUM_SPECIES];
  
  // Initialize the stemp array used in the loading process.
  for (i = 0; i < NUM_SPECIES/2; i++) 
    stemp[i] = ONE;
  for (i = NUM_SPECIES/2; i < NUM_SPECIES; i++)
    stemp[i] = RCONST(0.00001);

  for (jy = 0; jy < MY; jy++) 
  {
    for (jx = 0; jx < MX; jx++) 
    {
      sloc = IJ_Vptr(x_scale,jx,jy);
      for (i = 0; i < NUM_SPECIES; i++) 
        sloc[i] = stemp[i];
    }
  }
}

static void foodweb_set_F_scale(void* context, real_t* F_scale)
{
  // We scale F the same way we scale x.
  foodweb_set_x_scale(context, F_scale);
}

static void foodweb_set_constraints(void* context, real_t* constraints)
{
  // Enforce positivity on all components.
  for (int i = 0; i < NEQ; ++i)
    constraints[i] = 2.0;
}

// Constructor for food web integrator with no preconditioner.
static newton_solver_t* foodweb_solver_new()
{
  // Set up a nonlinear integrator using GMRES with no globalization 
  // strategy.
  foodweb_t* data = foodweb_new();
  newton_solver_vtable vtable = {.eval = foodweb_func,
                                        .set_x_scale = foodweb_set_x_scale,
                                        .set_F_scale = foodweb_set_F_scale,
                                        .set_constraints = foodweb_set_constraints,
                                        .dtor = foodweb_dtor};
  newton_solver_t* integ = gmres_newton_solver_new("Food web",
                                                                 data,
                                                                 MPI_COMM_SELF,
                                                                 NEQ,
                                                                 vtable, 
                                                                 NO_GLOBAL_STRATEGY, 
                                                                 15, 2);
  return integ;
}

// Constructor for block-Jacobi-preconditioned food web integrator.
newton_solver_t* block_jacobi_precond_foodweb_solver_new()
{
  newton_solver_t* integ = foodweb_solver_new();
  foodweb_t* data = newton_solver_context(integ);
  int block_size = NUM_SPECIES;
  newton_pc_t* precond = block_jacobi_cpr_pc_from_function(MPI_COMM_WORLD, data, foodweb_func, NULL, data->sparsity, NEQ/block_size, 0, block_size);
  newton_solver_set_preconditioner(integ, precond);
  return integ;
}

// Constructor for LU-preconditioned food web integrator.
newton_solver_t* lu_precond_foodweb_solver_new()
{
  newton_solver_t* integ = foodweb_solver_new();
  foodweb_t* data = newton_solver_context(integ);
  newton_pc_t* precond = lu_cpr_pc_from_function(MPI_COMM_WORLD, data, foodweb_func, NULL, data->sparsity, NEQ, 0);
  newton_solver_set_preconditioner(integ, precond);
  return integ;
}

// Constructor for ILU-preconditioned food web integrator.
newton_solver_t* ilu_precond_foodweb_solver_new()
{
  newton_solver_t* integ = foodweb_solver_new();
  foodweb_t* data = newton_solver_context(integ);
  ilu_params_t* ilu_params = ilu_params_new();
  newton_pc_t* precond = ilu_cpr_pc_from_function(MPI_COMM_WORLD, data, foodweb_func, NULL, data->sparsity, NEQ, 0, ilu_params);
  newton_solver_set_preconditioner(integ, precond);
  return integ;
}

// Returns initial conditions.
real_t* foodweb_initial_conditions()
{
  real_t* cc = malloc(sizeof(real_t) * NEQ);
  int i, jx, jy;
  real_t *cloc;
  real_t  ctemp[NUM_SPECIES];
  
  for (i = 0; i < NUM_SPECIES/2; i++)
    ctemp[i] = PREYIN;
  for (i = NUM_SPECIES/2; i < NUM_SPECIES; i++) 
    ctemp[i] = PREDIN;

  for (jy = 0; jy < MY; jy++) 
  {
    for (jx = 0; jx < MX; jx++) 
    {
      cloc = IJ_Vptr(cc,jx,jy);
      for (i = 0; i < NUM_SPECIES; i++) 
        cloc[i] = ctemp[i];
    }
  }
  return cc;
}

