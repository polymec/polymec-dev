// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/polymec.h"
#include "integrators/bdf_ode_integrator.h"
#include "integrators/ark_ode_integrator.h"
#include "integrators/cpr_newton_pc.h"

//------------------------------------------------------------------------
//               Diurnal kinetic advection-diffusion problem
//------------------------------------------------------------------------
// This is adapted from the CVODE Diurnal example problem programmed by 
// Scott D. Cohen, Alan C. Hindmarsh and Radu Serban @ LLNL
// An ODE system is generated from the following 2-species diurnal
// kinetics advection-diffusion PDE system in 2 space dimensions:
//
// dc(i)/dt = Kh*(d/dx)^2 c(i) + V*dc(i)/dx + (d/dy)(Kv(y)*dc(i)/dy)
//                 + Ri(c1,c2,t)      for i = 1,2,   where
//   R1(c1,c2,t) = -q1*c1*c3 - q2*c1*c2 + 2*q3(t)*c3 + q4(t)*c2 ,
//   R2(c1,c2,t) =  q1*c1*c3 - q2*c1*c2 - q4(t)*c2 ,
//   Kv(y) = Kv0*exp(y/5) ,
// Kh, V, Kv0, q1, q2, and c3 are constants, and q3(t) and q4(t)
// vary diurnally. The problem is posed on the square
//   0 <= x <= 20,    30 <= y <= 50   (all in km),
// with homogeneous Neumann boundary conditions, and for time t in
//   0 <= t <= 86400 sec (1 day).
// The PDE system is treated by central differences on a uniform
// 10 x 10 mesh, with simple polynomial initial profiles.
// The problem is solved with CVODE, with the BDF/GMRES
// method (i.e. using the CVSPGMR linear solver) and the
// block-diagonal part of the Newton matrix as a left
// preconditioner. A copy of the block-diagonal part of the
// Jacobian is saved and conditionally reused within the Precond
// routine.

#include <cvode/cvode_spgmr.h>        /* prototypes & constants for CVSPGMR */
#include <nvector/nvector_serial.h>   /* serial N_Vector types, fct., macros */
#include <sundials/sundials_dense.h>  /* use generic dense solver in precond. */
#include <sundials/sundials_types.h>  /* definition of real_t */
#include <sundials/sundials_math.h>   /* contains the macros SUNABS, SUNSQR, SUNRexp*/

/* Problem Constants */

#define ZERO RCONST(0.0)
#define ONE  RCONST(1.0)
#define TWO  RCONST(2.0)

#define NUM_SPECIES  2                 /* number of species         */
#define KH           RCONST(4.0e-6)    /* horizontal diffusivity Kh */
#define VEL          RCONST(0.001)     /* advection velocity V      */
#define KV0          RCONST(1.0e-8)    /* coefficient in Kv(y)      */
#define Q1           RCONST(1.63e-16)  /* coefficients q1, q2, c3   */ 
#define Q2           RCONST(4.66e-16)
#define C3           RCONST(3.7e16)
#define A3           RCONST(22.62)     /* coefficient in expression for q3(t) */
#define A4           RCONST(7.601)     /* coefficient in expression for q4(t) */
#define C1_SCALE     RCONST(1.0e6)     /* coefficients in initial profiles    */
#define C2_SCALE     RCONST(1.0e12)

#define T0           ZERO                 /* initial time */
#define NOUT         12                   /* number of output times */
#define TWOHR        RCONST(7200.0)       /* number of seconds in two hours  */
#define HALFDAY      RCONST(4.32e4)       /* number of seconds in a half day */
#define PI       RCONST(3.1415926535898)  /* pi */ 

#define XMIN         ZERO                 /* grid boundaries in x  */
#define XMAX         RCONST(20.0)           
#define YMIN         RCONST(30.0)         /* grid boundaries in y  */
#define YMAX         RCONST(50.0)
#define XMID         RCONST(10.0)         /* grid midpoints in x,y */          
#define YMID         RCONST(40.0)

#define MX           10             /* MX = number of x mesh points */
#define MY           10             /* MY = number of y mesh points */
#define NSMX         20             /* NSMX = NUM_SPECIES*MX */
#define MM           (MX*MY)        /* MM = MX*MY */

/* CVodeInit Constants */

#define RTOL    RCONST(1.0e-5)    /* scalar relative tolerance */
#define FLOOR   RCONST(100.0)     /* value of C1 or C2 at which tolerances */
                                  /* change from relative to absolute      */
#define ATOL    (RTOL*FLOOR)      /* scalar absolute tolerance */
#define NEQ     (NUM_SPECIES*MM)  /* NEQ = number of equations */

/* User-defined vector and matrix accessor macros: IJKth, IJth */

/* IJKth is defined in order to isolate the translation from the
   mathematical 3-dimensional structure of the dependent variable vector
   to the underlying 1-dimensional storage. IJth is defined in order to
   write code which indexes into small dense matrices with a (row,column)
   pair, where 1 <= row, column <= NUM_SPECIES.   
   
   IJKth(vdata,i,j,k) references the element in the vdata array for
   species i at mesh point (j,k), where 1 <= i <= NUM_SPECIES,
   0 <= j <= MX-1, 0 <= k <= MY-1. The vdata array is obtained via
   the macro call vdata = NV_DATA_S(v), where v is an N_Vector. 
   For each mesh point (j,k), the elements for species i and i+1 are
   contiguous within vdata.

   IJth(a,i,j) references the (i,j)th entry of the small matrix real_t **a,
   where 1 <= i,j <= NUM_SPECIES. The small matrix routines in sundials_dense.h
   work with matrices stored by column in a 2-dimensional array. In C,
   arrays are indexed starting at 0, not 1. */

#define IJKth(vdata,i,j,k) (vdata[i-1 + (j)*NUM_SPECIES + (k)*NSMX])
#define IJth(a,i,j)        (a[j-1][i-1])

typedef struct 
{
  real_t q4, om, dx, dy, hdco, haco, vdco;

  // Sparsity graph.
  adj_graph_t* sparsity;
} diurnal_t;

// Newly initialized data context.
static diurnal_t* diurnal_new()
{
  diurnal_t* data = malloc(sizeof(diurnal_t));

  data->om = PI/HALFDAY;
  data->dx = (XMAX-XMIN)/(MX-1);
  data->dy = (YMAX-YMIN)/(MY-1);
  data->hdco = KH/SUNSQR(data->dx);
  data->haco = VEL/(TWO*data->dx);
  data->vdco = (ONE/SUNSQR(data->dy))*KV0;

  // Construct a sparsity graph.
  adj_graph_t* sparsity = adj_graph_new(MPI_COMM_SELF, MX*MY);
  for (int jy = 0; jy < MY; jy++) 
  {
    // Set lower/upper index shifts, special at boundaries. 
    int idyd = (jy == 0   ) ?  MX : -MX;
    int idyu = (jy == MY-1) ?  -MX : MX;
    
    for (int jx = 0; jx < MX; jx++) 
    {
      // Set left/right index shifts, special at boundaries. 
      int idxl = (jx ==  0  ) ?  1 : -1;
      int idxr = (jx == MX-1) ?  -1 : 1;

      // Find the edges for the vertex corresponding to (jx, jy).
      int idx = jx + MX*jy;
      int num_edges = 0;
      int edges[4];
      if (jy > 0)
        edges[num_edges++] = idx + idyd; // lower
      if (jy < (MY-1))
        edges[num_edges++] = idx + idyu; // upper
      if (jx > 0)
        edges[num_edges++] = idx + idxl; // left
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

static int diurnal_rhs(void* context, real_t t, real_t* u, real_t* udot)
{
  real_t q3, c1, c2, c1dn, c2dn, c1up, c2up, c1lt, c2lt;
  real_t c1rt, c2rt, cydn, cyup, hord1, hord2, horad1, horad2;
  real_t qq1, qq2, qq3, qq4, rkin1, rkin2, s, vertd1, vertd2, ydn, yup;
  real_t q4coef, dely, verdco, hordco, horaco;
  int jx, jy, idn, iup, ileft, iright;
  diurnal_t* data = context;

  // We don't bother with FPE for now.
  polymec_suspend_fpe();

  // Set diurnal rate coefficients. 

  s = sin(data->om*t);
  if (s > ZERO) {
    q3 = SUNRexp(-A3/s);
    data->q4 = SUNRexp(-A4/s);
  } else {
      q3 = ZERO;
      data->q4 = ZERO;
  }

  /* Make local copies of problem variables, for efficiency. */

  q4coef = data->q4;
  dely = data->dy;
  verdco = data->vdco;
  hordco  = data->hdco;
  horaco  = data->haco;

  /* Loop over all grid points. */

  for (jy=0; jy < MY; jy++) {

    /* Set vertical diffusion coefficients at jy +- 1/2 */

    ydn = YMIN + (jy - RCONST(0.5))*dely;
    yup = ydn + dely;
    cydn = verdco*SUNRexp(RCONST(0.2)*ydn);
    cyup = verdco*SUNRexp(RCONST(0.2)*yup);
    idn = (jy == 0) ? 1 : -1;
    iup = (jy == MY-1) ? -1 : 1;
    for (jx=0; jx < MX; jx++) {

      /* Extract c1 and c2, and set kinetic rate terms. */

      c1 = IJKth(u,1,jx,jy); 
      c2 = IJKth(u,2,jx,jy);
      qq1 = Q1*c1*C3;
      qq2 = Q2*c1*c2;
      qq3 = q3*C3;
      qq4 = q4coef*c2;
      rkin1 = -qq1 - qq2 + TWO*qq3 + qq4;
      rkin2 = qq1 - qq2 - qq4;

      /* Set vertical diffusion terms. */

      c1dn = IJKth(u,1,jx,jy+idn);
      c2dn = IJKth(u,2,jx,jy+idn);
      c1up = IJKth(u,1,jx,jy+iup);
      c2up = IJKth(u,2,jx,jy+iup);
      vertd1 = cyup*(c1up - c1) - cydn*(c1 - c1dn);
      vertd2 = cyup*(c2up - c2) - cydn*(c2 - c2dn);

      /* Set horizontal diffusion and advection terms. */

      ileft = (jx == 0) ? 1 : -1;
      iright =(jx == MX-1) ? -1 : 1;
      c1lt = IJKth(u,1,jx+ileft,jy); 
      c2lt = IJKth(u,2,jx+ileft,jy);
      c1rt = IJKth(u,1,jx+iright,jy);
      c2rt = IJKth(u,2,jx+iright,jy);
      hord1 = hordco*(c1rt - TWO*c1 + c1lt);
      hord2 = hordco*(c2rt - TWO*c2 + c2lt);
      horad1 = horaco*(c1rt - c1lt);
      horad2 = horaco*(c2rt - c2lt);

      /* Load all terms into udot. */

      IJKth(udot, 1, jx, jy) = vertd1 + hord1 + horad1 + rkin1; 
      IJKth(udot, 2, jx, jy) = vertd2 + hord2 + horad2 + rkin2;
    }
  }

  // Resume FPE.
  polymec_restore_fpe();

  return 0;
}

// Data context destructor.
static void diurnal_dtor(void* context)
{
  diurnal_t* data = context;
  adj_graph_free(data->sparsity);
  free(data);
}

// Returns initial conditions.
real_t* diurnal_initial_conditions(ode_integrator_t* integ)
{
  diurnal_t* data = bdf_ode_integrator_context(integ);
  real_t dx = data->dx, dy = data->dy;
  int jx, jy;
  real_t x, y, cx, cy;
  real_t* udata = malloc(sizeof(real_t) * NEQ);

  // Load initial profiles of c1 and c2 into u vector 

  for (jy=0; jy < MY; jy++) {
    y = YMIN + jy*dy;
    cy = SUNSQR(RCONST(0.1)*(y - YMID));
    cy = ONE - cy + RCONST(0.5)*SUNSQR(cy);
    for (jx=0; jx < MX; jx++) {
      x = XMIN + jx*dx;
      cx = SUNSQR(RCONST(0.1)*(x - XMID));
      cx = ONE - cx + RCONST(0.5)*SUNSQR(cx);
      IJKth(udata,1,jx,jy) = C1_SCALE*cx*cy; 
      IJKth(udata,2,jx,jy) = C2_SCALE*cx*cy;
    }
  }
  return udata; 
}

// Constructor for a BDF diurnal integrator with the given preconditioner.
static ode_integrator_t* bdf_diurnal_integrator_new(diurnal_t* data, newton_pc_t* precond)
{
  // Set up a time integrator using GMRES with a maximum order of 2 and 
  // a Krylov space of maximum dimension 5.
  ode_integrator_t* integ = jfnk_bdf_ode_integrator_new(2, MPI_COMM_SELF, NEQ, 0,
                                                        data, diurnal_rhs, NULL, 
                                                        diurnal_dtor, precond, 
                                                        JFNK_BDF_GMRES, 5);

  return integ;
}

// Constructor for block-Jacobi-preconditioned BDF diurnal integrator.
ode_integrator_t* block_jacobi_precond_bdf_diurnal_integrator_new()
{
  diurnal_t* data = diurnal_new();
  newton_pc_t* precond = block_jacobi_cpr_newton_pc_from_function(MPI_COMM_WORLD, data, diurnal_rhs, NULL, data->sparsity, NEQ/NUM_SPECIES, 0, NUM_SPECIES);
  ode_integrator_t* integ = bdf_diurnal_integrator_new(data, precond);
  return integ;
}

// Constructor for LU-preconditioned BDF diurnal integrator.
ode_integrator_t* lu_precond_bdf_diurnal_integrator_new()
{
  diurnal_t* data = diurnal_new();
  newton_pc_t* precond = lu_cpr_newton_pc_from_function(MPI_COMM_WORLD, data, diurnal_rhs, NULL, data->sparsity, NEQ, 0);
  ode_integrator_t* integ = bdf_diurnal_integrator_new(data, precond);
  return integ;
}

// Constructor for ILU-preconditioned BDF diurnal integrator.
ode_integrator_t* ilu_precond_bdf_diurnal_integrator_new()
{
  diurnal_t* data = diurnal_new();
  ilu_params_t* ilu_params = ilu_params_new();
  newton_pc_t* precond = ilu_cpr_newton_pc_from_function(MPI_COMM_WORLD, data, diurnal_rhs, NULL, data->sparsity, NEQ, 0, ilu_params);
  ode_integrator_t* integ = bdf_diurnal_integrator_new(data, precond);
  return integ;
}

// Constructor for an ARK diurnal integrator with the given preconditioner.
static ode_integrator_t* ark_diurnal_integrator_new(diurnal_t* data, newton_pc_t* precond)
{
  // Set up a time integrator using GMRES with a maximum order of 2 and 
  // a Krylov space of maximum dimension 5.
  ode_integrator_t* integ;
  if (precond != NULL)
  {
    integ = jfnk_ark_ode_integrator_new(3, MPI_COMM_SELF, NEQ, 0,
                                        data, NULL, diurnal_rhs, false, false, 
                                        NULL, NULL, diurnal_dtor, precond, 
                                        JFNK_ARK_GMRES, 5);
  }
  else
  {
    integ = functional_ark_ode_integrator_new(3, MPI_COMM_SELF, NEQ, 0,
                                              data, NULL, diurnal_rhs, NULL, 
                                              diurnal_dtor, 15);
  }

  return integ;
}

// Constructor for functional ARK diurnal integrator.
ode_integrator_t* functional_ark_diurnal_integrator_new()
{
  diurnal_t* data = diurnal_new();
  ode_integrator_t* integ = ark_diurnal_integrator_new(data, NULL);
  return integ;
}

// Constructor for block-Jacobi-preconditioned ARK diurnal integrator.
ode_integrator_t* block_jacobi_precond_ark_diurnal_integrator_new()
{
  diurnal_t* data = diurnal_new();
  newton_pc_t* precond = block_jacobi_cpr_newton_pc_from_function(MPI_COMM_WORLD, data, diurnal_rhs, NULL, data->sparsity, NEQ/NUM_SPECIES, 0, NUM_SPECIES);
  ode_integrator_t* integ = ark_diurnal_integrator_new(data, precond);
  return integ;
}

// Constructor for LU-preconditioned ARK diurnal integrator.
ode_integrator_t* lu_precond_ark_diurnal_integrator_new()
{
  diurnal_t* data = diurnal_new();
  newton_pc_t* precond = lu_cpr_newton_pc_from_function(MPI_COMM_WORLD, data, diurnal_rhs, NULL, data->sparsity, NEQ, 0);
  ode_integrator_t* integ = ark_diurnal_integrator_new(data, precond);
  return integ;
}

// Constructor for ILU-preconditioned ARK diurnal integrator.
ode_integrator_t* ilu_precond_ark_diurnal_integrator_new()
{
  diurnal_t* data = diurnal_new();
  ilu_params_t* ilu_params = ilu_params_new();
  newton_pc_t* precond = ilu_cpr_newton_pc_from_function(MPI_COMM_WORLD, data, diurnal_rhs, NULL, data->sparsity, NEQ, 0, ilu_params);
  ode_integrator_t* integ = ark_diurnal_integrator_new(data, precond);
  return integ;
}

