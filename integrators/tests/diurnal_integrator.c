// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/polymec.h"
#include "core/options.h"
#include "core/declare_nd_array.h"
#include "integrators/bdf_ode_integrator.h"
#include "integrators/ark_ode_integrator.h"
#include "integrators/bj_newton_pc.h"

#include <setjmp.h>
#include "cmocka.h"

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

#include <cvode/cvode_spgmr.h>        // prototypes & constants for CVSPGMR 
#include <nvector/nvector_serial.h>   // serial N_Vector types, fct., macros 
#include <sundials/sundials_dense.h>  // use generic dense solver in precond. 
#include <sundials/sundials_types.h>  // definition of real_t 
#include <sundials/sundials_math.h>   // contains the macros SUNABS, SUNSQR, SUNRexp

// Problem Constants 

#define ZERO RCONST(0.0)
#define ONE  RCONST(1.0)
#define TWO  RCONST(2.0)

#define NUM_SPECIES  2                 // number of species         
#define KH           RCONST(4.0e-6)    // horizontal diffusivity Kh 
#define VEL          RCONST(0.001)     // advection velocity V      
#define KV0          RCONST(1.0e-8)    // coefficient in Kv(y)      
#define Q1           RCONST(1.63e-16)  // coefficients q1, q2, c3   
#define Q2           RCONST(4.66e-16)
#define C3           RCONST(3.7e16)
#define A3           RCONST(22.62)     // coefficient in expression for q3(t) 
#define A4           RCONST(7.601)     // coefficient in expression for q4(t) 
#define C1_SCALE     RCONST(1.0e6)     // coefficients in initial profiles    
#define C2_SCALE     RCONST(1.0e12)

#define T0           ZERO                 // initial time 
#define NOUT         12                   // number of output times 
#define TWOHR        RCONST(7200.0)       // number of seconds in two hours  
#define HALFDAY      RCONST(4.32e4)       // number of seconds in a half day 
#define PI       RCONST(3.1415926535898)  // pi 

#define XMIN         ZERO                 // grid boundaries in x  
#define XMAX         RCONST(20.0)           
#define YMIN         RCONST(30.0)         // grid boundaries in y  
#define YMAX         RCONST(50.0)
#define XMID         RCONST(10.0)         // grid midpoints in x,y 
#define YMID         RCONST(40.0)

#define MX           10             // MX = number of x mesh points 
#define MY           10             // MY = number of y mesh points 
#define NSMX         20             // NSMX = NUM_SPECIES*MX 
#define MM           (MX*MY)        // MM = MX*MY 

// CVodeInit Constants 

#define RTOL    RCONST(1.0e-5)    // scalar relative tolerance 
#define FLOOR   RCONST(100.0)     // value of C1 or C2 at which tolerances 
                                  // change from relative to absolute      
#define ATOL    (RTOL*FLOOR)      // scalar absolute tolerance 
#define NEQ     (NUM_SPECIES*MM)  // NEQ = number of equations 

typedef struct 
{
  real_t q4, om, dx, dy, hdco, haco, vdco;

  // Sparsity graph.
  adj_graph_t* sparsity;
} diurnal_t;

// Newly initialized data context.
static diurnal_t* diurnal_new()
{
  diurnal_t* data = polymec_malloc(sizeof(diurnal_t));

  data->om = PI/HALFDAY;
  data->dx = (XMAX-XMIN)/(MX-1);
  data->dy = (YMAX-YMIN)/(MY-1);
  data->hdco = KH/SUNSQR(data->dx);
  data->haco = VEL/(TWO*data->dx);
  data->vdco = (ONE/SUNSQR(data->dy))*KV0;

  // Construct a sparsity graph.
  data->sparsity = adj_graph_new(MPI_COMM_SELF, MX*MY*NUM_SPECIES);
  for (int i = 0; i < MX; ++i) 
  {
    int i_left = (i == 0) ? 1 : -1;
    int i_right =(i == MX-1) ? -1 : 1;
    for (int j = 0; j < MY; ++j) 
    {
      int j_down = (j == 0) ? 1 : -1;
      int j_up = (j == MY-1) ? -1 : 1;
      for (int s = 0; s < 2; ++s)
      {
        int i_self = ARRAY_INDEX_3D(MX, MY, NUM_SPECIES, i, j, s);
        // Set the edges within the sparsity graph.
        int num_edges = 0;
        int edges[5];
        edges[num_edges++] = ARRAY_INDEX_3D(MX, MY, NUM_SPECIES, i, j, (s+1) % NUM_SPECIES);
        if (j > 0)
          edges[num_edges++] = ARRAY_INDEX_3D(MX, MY, NUM_SPECIES, i, j+j_down, s);
        if (j < (MY-1))
          edges[num_edges++] = ARRAY_INDEX_3D(MX, MY, NUM_SPECIES, i, j+j_up, s);
        if (i > 0)
          edges[num_edges++] = ARRAY_INDEX_3D(MX, MY, NUM_SPECIES, i+i_left, j, s);
        if (i < (MX-1))
          edges[num_edges++] = ARRAY_INDEX_3D(MX, MY, NUM_SPECIES, i+i_right, j, s);
        adj_graph_set_num_edges(data->sparsity, i_self, num_edges);
        memcpy(adj_graph_edges(data->sparsity, i_self), edges, sizeof(int) * num_edges);
      }
    }
  }

  return data;
}

static int diurnal_rhs(void* context, real_t t, real_t* U, real_t* U_dot)
{
  diurnal_t* data = context;
  DECLARE_3D_ARRAY(real_t, U_ijk, U, MX, MY, NUM_SPECIES);
  DECLARE_3D_ARRAY(real_t, U_dot_ijk, U_dot, MX, MY, NUM_SPECIES);

  // We don't bother with FPE for now.
  polymec_suspend_fpe();

  // Set diurnal rate coefficients. 
  real_t s = sin(data->om*t);
  real_t q3;
  if (s > ZERO) 
  {
    q3 = SUNRexp(-A3/s);
    data->q4 = SUNRexp(-A4/s);
  } 
  else 
  {
    q3 = ZERO;
    data->q4 = ZERO;
  }

  // Make local copies of problem variables, for efficiency. 
  real_t q4coef = data->q4;
  real_t dely = data->dy;
  real_t verdco = data->vdco;
  real_t hordco  = data->hdco;
  real_t horaco  = data->haco;

  // Loop over all grid points. 
  for (int i=0; i < MX; i++) 
  {
    int ileft = (i == 0) ? 1 : -1;
    int iright =(i == MX-1) ? -1 : 1;
    for (int j=0; j < MY; ++j) 
    {
      // Set vertical diffusion coefficients at j +- 1/2 
      real_t ydn = YMIN + (j - RCONST(0.5))*dely;
      real_t yup = ydn + dely;
      real_t cydn = verdco*SUNRexp(RCONST(0.2)*ydn);
      real_t cyup = verdco*SUNRexp(RCONST(0.2)*yup);
      int idn = (j == 0) ? 1 : -1;
      int iup = (j == MY-1) ? -1 : 1;

      // Extract c1 and c2, and set kinetic rate terms. 
      real_t c1 = U_ijk[i][j][0]; 
      real_t c2 = U_ijk[i][j][1];
      real_t qq1 = Q1*c1*C3;
      real_t qq2 = Q2*c1*c2;
      real_t qq3 = q3*C3;
      real_t qq4 = q4coef*c2;
      real_t rkin1 = -qq1 - qq2 + TWO*qq3 + qq4;
      real_t rkin2 = qq1 - qq2 - qq4;

      // Set vertical diffusion terms. 
      real_t c1dn = U_ijk[i][j+idn][0];
      real_t c2dn = U_ijk[i][j+idn][1];
      real_t c1up = U_ijk[i][j+iup][0];
      real_t c2up = U_ijk[i][j+iup][1];
      real_t vertd1 = cyup*(c1up - c1) - cydn*(c1 - c1dn);
      real_t vertd2 = cyup*(c2up - c2) - cydn*(c2 - c2dn);

      // Set horizontal diffusion and advection terms. 
      real_t c1lt = U_ijk[i+ileft][j][0]; 
      real_t c2lt = U_ijk[i+ileft][j][1];
      real_t c1rt = U_ijk[i+iright][j][0];
      real_t c2rt = U_ijk[i+iright][j][1];
      real_t hord1 = hordco*(c1rt - TWO*c1 + c1lt);
      real_t hord2 = hordco*(c2rt - TWO*c2 + c2lt);
      real_t horad1 = horaco*(c1rt - c1lt);
      real_t horad2 = horaco*(c2rt - c2lt);

      // Load all terms into udot. 
      U_dot_ijk[i][j][0] = vertd1 + hord1 + horad1 + rkin1; 
      U_dot_ijk[i][j][1] = vertd2 + hord2 + horad2 + rkin2;
    }
  }

  // Resume FPE.
  polymec_restore_fpe();

  return 0;
}

// Function for accumulating a Jacobian column value into a column map.
static void accumulate_J_value(index_real_unordered_map_t* col_map, 
                               index_t column, 
                               real_t value)
{
  real_t* val_ptr = index_real_unordered_map_get(col_map, column);
  if (val_ptr == NULL) // New value is inserted.
    index_real_unordered_map_insert(col_map, column, value);
  else // accumulate the given value into the existing one.
    index_real_unordered_map_insert(col_map, column, value + *val_ptr);
}

// Function for inserting a row of values into the Jacobian matrix.
static void insert_J_values(index_t row, 
                            index_real_unordered_map_t* col_map, 
                            krylov_matrix_t* J)
{
  index_t num_cols = (index_t)col_map->size;
  index_t indices[num_cols];
  real_t values[num_cols];
  int pos = 0, k = 0;
  while (index_real_unordered_map_next(col_map, &pos, &indices[k], &values[k])) ++k;
  krylov_matrix_set_values(J, 1, &num_cols, &row, indices, values);
  index_real_unordered_map_clear(col_map);
}

// Function for constructing the Jacobian matrix for the diurnal system.
static int diurnal_J(void* context, real_t t, real_t* U, real_t* U_dot, krylov_matrix_t* J)
{
  diurnal_t* data = context;
  DECLARE_3D_ARRAY(real_t, U_ijk, U, MX, MY, NUM_SPECIES);

  // We don't bother with FPE for now.
  polymec_suspend_fpe();

  // Set diurnal rate coefficients. 
  real_t s = sin(data->om*t);
  if (s > ZERO) 
    data->q4 = SUNRexp(-A4/s);
  else 
    data->q4 = ZERO;

  // Make local copies of problem variables, for efficiency. 
  real_t q4coef = data->q4;
  real_t dely = data->dy;
  real_t verdco = data->vdco;
  real_t hordco  = data->hdco;
  real_t horaco  = data->haco;

  // Maps from indices to their values in the Jacobian.
  index_real_unordered_map_t* I1_map = index_real_unordered_map_new();
  index_real_unordered_map_t* I2_map = index_real_unordered_map_new();

  // Loop over all grid points. 
  for (int i = 0; i < MX; ++i) 
  {
    int i_left = (i == 0) ? 1 : -1;
    int i_right =(i == MX-1) ? -1 : 1;

    for (int j = 0; j < MY; ++j) 
    {
      // Set vertical diffusion coefficients at j +- 1/2 
      real_t ydn = YMIN + (j - RCONST(0.5))*dely;
      real_t yup = ydn + dely;
      real_t cydn = verdco*SUNRexp(RCONST(0.2)*ydn);
      real_t cyup = verdco*SUNRexp(RCONST(0.2)*yup);

      int i_down = (j == 0) ? 1 : -1;
      int i_up = (j == MY-1) ? -1 : 1;

      // There are up to 12 Jacobian contributions: 6 for each of 2 species: 
      // 1 in each of the 5 stencil points, plus 1 reaction term in which the 
      // species interact directly with one another.
      real_t J1_self = 0.0, J1_rxn = 0.0, J1_left = 0.0, J1_right = 0.0, J1_up = 0.0, J1_down = 0.0;
      real_t J2_self = 0.0, J2_rxn = 0.0, J2_left = 0.0, J2_right = 0.0, J2_up = 0.0, J2_down = 0.0;
      index_t I1_self  = ARRAY_INDEX_3D(MX, MY, NUM_SPECIES, i, j, 0), 
              I1_left  = ARRAY_INDEX_3D(MX, MY, NUM_SPECIES, i+i_left, j, 0),
              I1_right = ARRAY_INDEX_3D(MX, MY, NUM_SPECIES, i+i_right, j, 0),
              I1_up    = ARRAY_INDEX_3D(MX, MY, NUM_SPECIES, i, j+i_up, 0),
              I1_down  = ARRAY_INDEX_3D(MX, MY, NUM_SPECIES, i, j+i_down, 0),
              I2_self  = ARRAY_INDEX_3D(MX, MY, NUM_SPECIES, i, j, 1), 
              I2_left  = ARRAY_INDEX_3D(MX, MY, NUM_SPECIES, i+i_left, j, 1),
              I2_right = ARRAY_INDEX_3D(MX, MY, NUM_SPECIES, i+i_right, j, 1),
              I2_up    = ARRAY_INDEX_3D(MX, MY, NUM_SPECIES, i, j+i_up, 1),
              I2_down  = ARRAY_INDEX_3D(MX, MY, NUM_SPECIES, i, j+i_down, 1);

      // Extract c1 and c2 at the current location.
      real_t c1 = U_ijk[i][j][0];
      real_t c2 = U_ijk[i][j][1];

      // Add in kinetic rate terms. 
      J1_self += -(Q1*C3 + Q2*c2);
      J1_rxn  +=  (q4coef - Q2*c1);
      J2_rxn  +=  (Q1*C3 - Q2*c2);
      J2_self += -(q4coef + Q2*c1);

      // Add in vertical diffusion terms.
      J1_self += -(cyup+cydn);
      J1_up   += cyup;
      J1_down += cydn;
      J2_self += -(cyup+cydn);
      J2_up   += cyup;
      J2_down += cydn;

      // Add in horizontal diffusion terms.
      J1_self  -= TWO*hordco;
      J1_left  += hordco;
      J1_right += hordco;
      J2_self  -= TWO*hordco;
      J2_left  += hordco;
      J2_right += hordco;

      // Add in horizontal advection terms.
      J1_left  += -horaco;
      J1_right +=  horaco;
      J2_left  += -horaco;
      J2_right +=  horaco;

      // Aggregate the values, since some of them are aliased on top of each 
      // other (owing to periodic boundary conditions).
      accumulate_J_value(I1_map, I1_self, J1_self);
      accumulate_J_value(I1_map, I2_self, J1_rxn);
      accumulate_J_value(I1_map, I1_left, J1_left);
      accumulate_J_value(I1_map, I1_right, J1_right);
      accumulate_J_value(I1_map, I1_up, J1_up);
      accumulate_J_value(I1_map, I1_down, J1_down);

      accumulate_J_value(I2_map, I2_self, J2_self);
      accumulate_J_value(I2_map, I1_self, J2_rxn);
      accumulate_J_value(I2_map, I2_left, J2_left);
      accumulate_J_value(I2_map, I2_right, J2_right);
      accumulate_J_value(I2_map, I2_up, J2_up);
      accumulate_J_value(I2_map, I2_down, J2_down);

      // Stick the data into the Jacobian matrix.
      insert_J_values(I1_self, I1_map, J);
      insert_J_values(I2_self, I2_map, J);
    }
  }

  // Resume FPE.
  polymec_restore_fpe();

  // Assemble the Jacobian matrix.
  krylov_matrix_assemble(J);

  // Clean up.
  index_real_unordered_map_free(I1_map);
  index_real_unordered_map_free(I2_map);

  return 0;
}

// Data context destructor.
static void diurnal_dtor(void* context)
{
  diurnal_t* data = context;
  adj_graph_free(data->sparsity);
  polymec_free(data);
}

// Returns initial conditions.
real_t* diurnal_initial_conditions(ode_integrator_t* integ);
real_t* diurnal_initial_conditions(ode_integrator_t* integ)
{
  // Get the data for the diurnal system.
  diurnal_t* data;
  const char* integ_name = ode_integrator_name(integ);
  if (string_contains(integ_name, "INK Backwards-Difference"))
    data = ink_bdf_ode_integrator_context(integ);
  else
    data = bdf_ode_integrator_context(integ);
  real_t dx = data->dx, dy = data->dy;
  real_t* u_data = polymec_malloc(sizeof(real_t) * NEQ);
  DECLARE_3D_ARRAY(real_t, u_ijk, u_data, MX, MY, NUM_SPECIES);

  // Load initial profiles of c1 and c2 into u vector 
  for (int j=0; j < MY; j++) 
  {
    real_t y = YMIN + j*dy;
    real_t cy = SUNSQR(RCONST(0.1)*(y - YMID));
    cy = ONE - cy + RCONST(0.5)*SUNSQR(cy);
    for (int i=0; i < MX; i++) 
    {
      real_t x = XMIN + i*dx;
      real_t cx = SUNSQR(RCONST(0.1)*(x - XMID));
      cx = ONE - cx + RCONST(0.5)*SUNSQR(cx);
      u_ijk[i][j][0] = C1_SCALE*cx*cy; 
      u_ijk[i][j][1] = C2_SCALE*cx*cy;
    }
  }
  return u_data; 
}

// Constructor for a BDF diurnal integrator with the given preconditioner.
static ode_integrator_t* jfnk_bdf_diurnal_integrator_new(diurnal_t* data, newton_pc_t* precond)
{
  // Set up a time integrator using GMRES with a maximum order of 5 and 
  // a Krylov space of maximum dimension 15.
  ode_integrator_t* integ = jfnk_bdf_ode_integrator_new(5, MPI_COMM_SELF, NEQ, 0,
                                                        data, diurnal_rhs, NULL, 
                                                        diurnal_dtor, precond, 
                                                        JFNK_BDF_GMRES, 15);

  return integ;
}

// Constructor for block-Jacobi-preconditioned BDF diurnal integrator.
ode_integrator_t* bj_jfnk_bdf_diurnal_integrator_new(newton_pc_side_t side);
ode_integrator_t* bj_jfnk_bdf_diurnal_integrator_new(newton_pc_side_t side)
{
  diurnal_t* data = diurnal_new();
  newton_pc_t* precond = cpr_bj_newton_pc_new(MPI_COMM_WORLD, data, diurnal_rhs, NULL, side, data->sparsity, NEQ/NUM_SPECIES, 0, NUM_SPECIES);
  ode_integrator_t* integ = jfnk_bdf_diurnal_integrator_new(data, precond);
  return integ;
}

// Constructor for an Inexact Newton-Krylov BDF diurnal integrator.
ode_integrator_t* ink_bdf_diurnal_integrator_new(krylov_factory_t* factory);
ode_integrator_t* ink_bdf_diurnal_integrator_new(krylov_factory_t* factory)
{
  diurnal_t* data = diurnal_new();
  // Set up a time integrator using GMRES with a maximum order of 5 and 
  // a Krylov space of maximum dimension 15.
  matrix_sparsity_t* J_sparsity = matrix_sparsity_from_graph(data->sparsity, NULL);
  ode_integrator_t* integ = ink_bdf_ode_integrator_new(5, MPI_COMM_SELF, factory,
                                                       J_sparsity, data, diurnal_rhs, 
                                                       diurnal_J, diurnal_dtor);
  ink_bdf_ode_integrator_use_gmres(integ, 15);
  return integ;
}

// Constructor for an ARK diurnal integrator with the given preconditioner.
static ode_integrator_t* jfnk_ark_diurnal_integrator_new(diurnal_t* data, newton_pc_t* precond)
{
  // Set up a time integrator using GMRES with a maximum order of 3 and 
  // a Krylov space of maximum dimension 15.
  ode_integrator_t* integ;
  if (precond != NULL)
  {
    integ = jfnk_ark_ode_integrator_new(3, MPI_COMM_SELF, NEQ, 0,
                                        data, NULL, diurnal_rhs, false, false, 
                                        NULL, NULL, diurnal_dtor, precond, 
                                        JFNK_ARK_GMRES, 15);
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
ode_integrator_t* functional_ark_diurnal_integrator_new(void);
ode_integrator_t* functional_ark_diurnal_integrator_new()
{
  diurnal_t* data = diurnal_new();
  ode_integrator_t* integ = jfnk_ark_diurnal_integrator_new(data, NULL);
  return integ;
}

// Constructor for block-Jacobi-preconditioned ARK diurnal integrator.
ode_integrator_t* bj_jfnk_ark_diurnal_integrator_new(newton_pc_side_t side);
ode_integrator_t* bj_jfnk_ark_diurnal_integrator_new(newton_pc_side_t side)
{
  diurnal_t* data = diurnal_new();
  newton_pc_t* precond = cpr_bj_newton_pc_new(MPI_COMM_WORLD, data, diurnal_rhs, NULL, side, data->sparsity, NEQ/NUM_SPECIES, 0, NUM_SPECIES);
  ode_integrator_t* integ = jfnk_ark_diurnal_integrator_new(data, precond);
  return integ;
}

// Test harness for ODE integrator.
int test_diurnal_step(void** state, ode_integrator_t* integ, int max_steps);
int test_diurnal_step(void** state, ode_integrator_t* integ, int max_steps)
{
  // Set up the problem.
#if POLYMEC_HAVE_DOUBLE_PRECISION
  bdf_ode_integrator_set_tolerances(integ, 1e-5, 1e-3);
#else
  bdf_ode_integrator_set_tolerances(integ, 1e-4, 1e-2);
#endif
  real_t* U = diurnal_initial_conditions(integ);
  DECLARE_3D_ARRAY(real_t, U_ijk, U, 10, 10, 2);

  // Are we asked to plot the result?
  options_t* opts = options_argv();
  char* plotfile = options_value(opts, "plotfile");
  int plot_I = -1, plot_J = -1;
  {
    char* I_str = options_value(opts, "I");
    if (string_is_integer(I_str))
      plot_I = atoi(I_str);
    char* J_str = options_value(opts, "J");
    if (string_is_integer(J_str))
      plot_J = atoi(J_str);
  }
  FILE* plot = NULL;
  if ((plotfile != NULL) && (plot_I >= 0) && (plot_J >= 0))
    plot = fopen(plotfile, "w");

  // Integrate it out to t = 86400 s (24 hours).
  real_t t = 0.0;
  int step = 0;
  while (t < 86400.0)
  {
    real_t t_old = t;
    bool integrated = ode_integrator_step(integ, 7200.0, &t, U);
    assert_true(integrated);
    log_detail("Step %d: t = %g, dt = %g", step, t, t - t_old);
    ++step;

    // Plot if requested.
    if (plot != NULL)
      fprintf(plot, "%g %g %g\n", t, U_ijk[plot_I][plot_J][0], U_ijk[plot_I][plot_J][1]);

    if (step >= max_steps)
      break;
  }
  printf("Final time: %g\n", t);
  bdf_ode_integrator_diagnostics_t diags;
  bdf_ode_integrator_get_diagnostics(integ, &diags);
  bdf_ode_integrator_diagnostics_fprintf(&diags, stdout);

  // If we've opened a plot file, close it.
  if (plot != NULL)
    fclose(plot);

  // Did the integrator finish before taking the max number of steps?
  assert_true(step < max_steps);

  ode_integrator_free(integ);
  free(U);
  return step;
}

