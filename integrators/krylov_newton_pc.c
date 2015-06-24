// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "sundials/sundials_spgmr.h"
#include "sundials/sundials_spbcgs.h"
#include "sundials/sundials_sptfqmr.h"
#include "core/sundials_helpers.h"
#include "integrators/krylov_newton_pc.h"

typedef struct
{
  MPI_Comm comm;
  void* context;
  SpgmrMemRec* gmr;
  SpbcgMemRec* bcg;
  SptfqmrMemRec* tfqmr;
  int N_local, N_remote, N_total;
  int (*F)(void* context, real_t t, real_t* x, real_t* F);
  int (*dae_F)(void* context, real_t t, real_t* x, real_t* x_dot, real_t* F);
  int (*Jv)(void* context, 
            real_t alpha, real_t beta, real_t gamma, 
            real_t t, real_t* x, real_t* v, real_t* Jv);
  int (*dae_Jv)(void* context, 
                real_t alpha, real_t beta, real_t gamma, 
                real_t t, real_t* x, real_t* x_dot, 
                real_t* v, real_t* Jv);
  int (*Atimes)(void* A_data, N_Vector v, N_Vector Av);
  void (*dtor)(void* context);
  int krylov_dim;

  // State information.
  real_t t;
  N_Vector x, x_dot;

  // Preconditioning parameters.
  real_t alpha, beta, gamma;

  // Bookkeeping.
  N_Vector z, r, z_scaling, r_scaling, temp1, temp2, temp3, temp4;
  real_t delta;
  int num_iters;

  // GMRES stuff.
  int max_restarts;
  gmres_gram_schmidt_t gs;
} krylov_pc_t;

static void krylov_compute_p(void* context, 
                             real_t alpha, real_t beta, real_t gamma, 
                             real_t t, real_t* x, real_t* xdot)
{
  krylov_pc_t* krylov = context;
  krylov->alpha = alpha;
  krylov->beta = beta;
  krylov->gamma = gamma;
  if ((alpha != 0.0) && (beta != 0.0) && (gamma != 0.0))
    log_debug("krylov_newton_pc: approximating J = %g * I + %g * dF/dx + %g * dF/d(xdot)...", alpha, beta, gamma);
  else if ((alpha == 0.0) && (beta != 0.0) && (gamma != 0.0))
    log_debug("krylov_newton_pc: approximating J = %g * dF/dx + %g * dF/d(xdot)...", beta, gamma);
  else if ((alpha == 0.0) && (beta == 0.0) && (gamma != 0.0))
    log_debug("krylov_newton_pc: approximating J = %g * dF/d(xdot)...", gamma);
  else if ((alpha != 0.0) && (beta != 0.0))
    log_debug("krylov_newton_pc: approximating J = %g * I + %g * dF/dx...", alpha, beta);
  else if ((alpha == 0.0) && (beta != 0.0))
    log_debug("krylov_newton_pc: approximating J = %g * dF/dx...", beta);
}

static void krylov_pc_dtor(void* context)
{
  krylov_pc_t* krylov = context;
  if ((krylov->dtor != NULL) && (krylov->context != NULL))
    krylov->dtor(krylov->context);
  if (krylov->gmr != NULL)
    SpgmrFree(krylov->gmr);
  if (krylov->bcg != NULL)
    SpbcgFree(krylov->bcg);
  if (krylov->tfqmr != NULL)
    SptfqmrFree(krylov->tfqmr);
  N_VDestroy(krylov->x);
  N_VDestroy(krylov->x_dot);
  N_VDestroy(krylov->z);
  N_VDestroy(krylov->r);
  if (krylov->z_scaling != NULL)
  {
    N_VDestroy(krylov->z_scaling);
  }
  if (krylov->r_scaling != NULL)
  {
    N_VDestroy(krylov->r_scaling);
  }
  if (krylov->temp1 != NULL)
  {
    N_VDestroy(krylov->temp1);
    N_VDestroy(krylov->temp2);
    N_VDestroy(krylov->temp3);
    N_VDestroy(krylov->temp4);
  }
  polymec_free(krylov);
}

static int dae_F_adaptor(void* context, real_t t, real_t* x, real_t* x_dot, real_t* F)
{
  // We are passed the actual preconditioner as our context pointer, so get the 
  // "real" one here.
  krylov_pc_t* krylov = context;
  ASSERT(krylov->F != NULL);
  return krylov->F(krylov->context, t, x, F);
}

static int dae_Jv_adaptor(void* context, 
                          real_t alpha, real_t beta, real_t gamma,
                          real_t t, real_t* x, real_t* x_dot, 
                          real_t* v, real_t* Jv)
{
  // We are passed the actual preconditioner as our context pointer, so get the 
  // "real" one here.
  krylov_pc_t* krylov = context;
  ASSERT(krylov->Jv != NULL);
  return krylov->Jv(krylov->context, alpha, beta, gamma, t, x, v, Jv);
}

// This approximates J*v using a difference quotient applied to a function F.
static int Atimes_DQ_adaptor(void* A_data, N_Vector v, N_Vector Av)
{
  krylov_pc_t* krylov = A_data;
  void* F_context = NULL;
  int (*F)(void* context, real_t t, real_t* x, real_t* x_dot, real_t* F) = NULL;
  if (krylov->F == NULL)
  {
    F_context = krylov->context;
    F = krylov->dae_F;
  }
  else
  {
    F_context = krylov;
    F = dae_F_adaptor;
  }

  // Allocate workspaces if needed.
  if (krylov->temp1 == NULL)
  {
    krylov->temp1 = N_VNew(krylov->comm, krylov->N_local);
    krylov->temp2 = N_VNew(krylov->comm, krylov->N_local);
    krylov->temp3 = N_VNew(krylov->comm, krylov->N_local);
    krylov->temp4 = N_VNew(krylov->comm, krylov->N_local);
  }

  // Evaluate the function F at its present value here, storing it in temp4.
  int status;
  if (krylov->x_dot != NULL)
  {
    status = F(F_context, 
               krylov->t, NV_DATA(krylov->x), NV_DATA(krylov->x_dot), 
               NV_DATA(krylov->temp4));
  }
  else
  {
    status = F(F_context, 
               krylov->t, NV_DATA(krylov->x), NULL,
               NV_DATA(krylov->temp4));
  }
  if (status != 0)
    return status;

  // We estimate sigma here, following Brown and Saad (1990), p 469.
  // See KINSpilsDQJTimes in kinsol_spils.c for a cheat sheet.
  if (krylov->z_scaling != NULL)
  {
    N_VProd(v, krylov->z_scaling, krylov->temp1);
    N_VProd(krylov->x, krylov->z_scaling, Av);
  }
  else
  {
    N_VScale(1.0, v, krylov->temp1);
    N_VScale(1.0, krylov->x, Av);
  }
  real_t sutsv = N_VDotProd(Av, krylov->temp1);
  real_t vtv = N_VDotProd(krylov->temp1, krylov->temp1);
  real_t sq1norm = N_VL1Norm(krylov->temp1);
  real_t sign = SIGN(sutsv);
  real_t sqrt_relfunc = sqrt(UNIT_ROUNDOFF);
  real_t sigma = sign * sqrt_relfunc * MAX(fabs(sutsv), sq1norm) / vtv;
  real_t sigma_inv = 1.0 / sigma;

  // Add the identity term into Av.
  N_VScale(krylov->alpha, v, Av);

  if (krylov->beta != 0.0)
  {
    // Compute the value of x at which to evaluate F -> temp1.
    N_VLinearSum(1.0, krylov->x, sigma, v, krylov->temp1);
    
    // Call F -> temp2.
    int status = F(F_context, 
                   krylov->t, NV_DATA(krylov->temp1), NV_DATA(krylov->x_dot), 
                   NV_DATA(krylov->temp2));
    if (status != 0)
      return status;

    // Store the difference quotient in temp3 and add it into Av.
    N_VLinearSum(sigma_inv, krylov->temp2, -sigma_inv, krylov->temp4, krylov->temp3);
    N_VLinearSum(1.0, Av, krylov->beta, krylov->temp3, Av);
  }

  if (krylov->gamma != 0.0)
  {
    ASSERT(krylov->x_dot != NULL);
    // Compute the value of x_dot at which to evaluate F -> temp1.
    N_VLinearSum(1.0, krylov->x_dot, sigma, v, krylov->temp1);
    
    // Call F -> temp2.
    int status = F(F_context, 
                   krylov->t, NV_DATA(krylov->x), NV_DATA(krylov->temp1), 
                   NV_DATA(krylov->temp2));
    if (status != 0)
      return status;

    // Store the difference quotient in temp3 and add it into Av.
    N_VLinearSum(sigma_inv, krylov->temp2, -sigma_inv, krylov->temp4, krylov->temp3);
    N_VLinearSum(1.0, Av, krylov->gamma, krylov->temp3, Av);
  }

  return 0; 
}

static int Atimes_Jv_adaptor(void* A_data, N_Vector v, N_Vector Av)
{
  krylov_pc_t* krylov = A_data;
  void* Jv_context = NULL;
  int (*Jv)(void* context, 
            real_t alpha, real_t beta, real_t gamma, 
            real_t t, real_t* x, real_t* x_dot, 
            real_t* v, real_t* Jv) = NULL;
  if (krylov->dae_Jv != NULL)
  {
    Jv_context = krylov->context;
    Jv = krylov->dae_Jv;
  }
  else
  {
    Jv_context = krylov;
    Jv = dae_Jv_adaptor;
  }
  return Jv(Jv_context, krylov->alpha, krylov->beta, krylov->gamma, 
            krylov->t, NV_DATA(krylov->x), NV_DATA(krylov->x_dot), 
            NV_DATA(v), NV_DATA(Av));
}

static krylov_pc_t* krylov_pc_new(MPI_Comm comm,
                                  void* context,
                                  int (*F)(void* context, real_t t, real_t* x, real_t* F),
                                  int (*Jv)(void* context, 
                                           real_t alpha, real_t beta, real_t gamma, 
                                           real_t t, real_t* x, real_t* v, real_t* Jv),
                                  void (*dtor)(void* context),
                                  int num_local_values, 
                                  int num_remote_values,
                                  int krylov_dim)
{
  ASSERT((F != NULL) || (Jv != NULL));
  ASSERT((F == NULL) || (Jv == NULL));
  ASSERT(num_local_values > 0);
  ASSERT(krylov_dim >= 3);
  krylov_pc_t* krylov = polymec_malloc(sizeof(krylov_pc_t));
  krylov->comm = comm;
  krylov->context = context;
  krylov->gmr = NULL;
  krylov->bcg = NULL;
  krylov->tfqmr = NULL;
  krylov->F = F;
  krylov->Jv = Jv;
  krylov->dae_F = NULL;
  krylov->dae_Jv = NULL;
  krylov->Atimes = (F != NULL) ? Atimes_DQ_adaptor : Atimes_Jv_adaptor;
  krylov->dtor = dtor;
  krylov->N_local = num_local_values;
  krylov->N_remote = num_remote_values;
  krylov->N_total = krylov->N_local + krylov->N_remote;
  krylov->t = 0.0;
  krylov->x = N_VNew(comm, krylov->N_local);
  krylov->x_dot = N_VNew(comm, krylov->N_local);
  krylov->z = N_VNew(comm, krylov->N_local);
  krylov->r = N_VNew(comm, krylov->N_local);
  krylov->z_scaling = NULL;
  krylov->r_scaling = NULL;
  krylov->temp1 = NULL;
  krylov->temp2 = NULL;
  krylov->temp3 = NULL;
  krylov->temp4 = NULL;
  krylov->delta = 1e-4;
  return krylov;
}

static krylov_pc_t* dae_krylov_pc_new(MPI_Comm comm,
                                      void* context,
                                      int (*F)(void* context, real_t t, real_t* x, real_t* x_dot, real_t* F),
                                      int (*Jv)(void* context, 
                                                real_t alpha, real_t beta, real_t gamma, 
                                                real_t t, real_t* x, real_t* x_dot, 
                                                real_t* v, real_t* Jv),
                                      void (*dtor)(void* context),
                                      int num_local_values, 
                                      int num_remote_values,
                                      int krylov_dim)
{
  ASSERT((F != NULL) || (Jv != NULL));
  ASSERT((F == NULL) || (Jv == NULL));
  ASSERT(num_local_values > 0);
  ASSERT(krylov_dim >= 3);
  krylov_pc_t* krylov = polymec_malloc(sizeof(krylov_pc_t));
  krylov->comm = comm;
  krylov->context = context;
  krylov->gmr = NULL;
  krylov->bcg = NULL;
  krylov->tfqmr = NULL;
  krylov->F = NULL;
  krylov->Jv = NULL;
  krylov->dae_F = F;
  krylov->dae_Jv = Jv;
  krylov->Atimes = (F != NULL) ? Atimes_DQ_adaptor : Atimes_Jv_adaptor;
  krylov->dtor = dtor;
  krylov->N_local = num_local_values;
  krylov->N_remote = num_remote_values;
  krylov->N_total = krylov->N_local + krylov->N_remote;
  krylov->t = 0.0;
  krylov->x = N_VNew(comm, krylov->N_local);
  krylov->x_dot = N_VNew(comm, krylov->N_local);
  krylov->z = N_VNew(comm, krylov->N_local);
  krylov->r = N_VNew(comm, krylov->N_local);
  krylov->z_scaling = NULL;
  krylov->r_scaling = NULL;
  krylov->temp1 = NULL;
  krylov->temp2 = NULL;
  krylov->temp3 = NULL;
  krylov->temp4 = NULL;
  krylov->delta = 1e-4;
  return krylov;
}

static bool gmres_solve(void* context, 
                        real_t t, real_t* x, real_t* x_dot,
                        real_t* r, real_t* z)
{
  krylov_pc_t* krylov = context;
  krylov->t = t;
  memcpy(NV_DATA(krylov->r), r, sizeof(real_t) * krylov->N_local);
  memcpy(NV_DATA(krylov->z), z, sizeof(real_t) * krylov->N_local);
  memcpy(NV_DATA(krylov->x), x, sizeof(real_t) * krylov->N_local);
  if (x_dot != NULL)
    memcpy(NV_DATA(krylov->x_dot), x, sizeof(real_t) * krylov->N_local);
  real_t res_norm;
  int nli, nps;
  int status = SpgmrSolve(krylov->gmr, krylov, krylov->z, krylov->r, PREC_NONE, krylov->gs,
                          krylov->delta, krylov->max_restarts, NULL, krylov->r_scaling, krylov->z_scaling,
                          krylov->Atimes, NULL, &res_norm, &nli, &nps);
  switch(status)
  {
    case SPGMR_RES_REDUCED:
      log_debug("  GMRES: Residual norm was reduced (%g), but not below tolerance (%g)", res_norm, krylov->delta);
      break;
    case SPGMR_CONV_FAIL:
      log_debug("  GMRES: Failed to converge.");
      break;
    case SPGMR_QRFACT_FAIL:
      log_debug("  GMRES: Singular matrix in QR factorization.");
      break;
    case SPGMR_ATIMES_FAIL_REC:
      log_debug("  GMRES: Jacobian-vector product failed recoverably.");
      break;
    case SPGMR_ATIMES_FAIL_UNREC:
      log_debug("  GMRES: Jacobian-vector product failed unrecoverably.");
      break;
    case SPGMR_GS_FAIL:
      log_debug("  GMRES: Gram-Schmidt orthogonalization failed.");
      break;
    case SPGMR_QRSOL_FAIL:
      log_debug("  GMRES: Singular matrix in QR solve.");
    default:
      break;
  }
  return (status == SPGMR_SUCCESS);
}

newton_pc_t* gmres_newton_pc_new(MPI_Comm comm,
                                 void* context,
                                 int (*F)(void* context, real_t t, real_t* x, real_t* F),
                                 int (*Jv)(void* context, 
                                           real_t alpha, real_t beta, real_t gamma, 
                                           real_t t, real_t* x, real_t* v, real_t* Jv),
                                 void (*dtor)(void* context),
                                 int num_local_values, 
                                 int num_remote_values,
                                 int krylov_dim,
                                 int max_restarts,
                                 gmres_gram_schmidt_t gs)
{
  ASSERT(max_restarts >= 0);
  krylov_pc_t* krylov = krylov_pc_new(comm, context, F, Jv, dtor, num_local_values, num_remote_values, krylov_dim);
  krylov->gmr = SpgmrMalloc(krylov_dim, krylov->z);
  krylov->max_restarts = max_restarts;
  krylov->gs = gs;
  newton_pc_vtable vtable = {.solve = gmres_solve,
                             .compute_p = krylov_compute_p,
                             .dtor = krylov_pc_dtor};
  return newton_pc_new("Krylov preconditioner (GMRES)", krylov, vtable);
}
                                        
newton_pc_t* dae_gmres_newton_pc_new(MPI_Comm comm,
                                     void* context,
                                     int (*F)(void* context, real_t t, real_t* x, real_t* x_dot, real_t* F),
                                     int (*Jv)(void* context, 
                                               real_t alpha, real_t beta, real_t gamma, 
                                               real_t t, real_t* x, real_t* x_dot, 
                                               real_t* v, real_t* Jv),
                                     void (*dtor)(void* context),
                                     int num_local_values, 
                                     int num_remote_values,
                                     int krylov_dim,
                                     int max_restarts,
                                     gmres_gram_schmidt_t gs)
{
  ASSERT(max_restarts >= 0);
  krylov_pc_t* krylov = dae_krylov_pc_new(comm, context, F, Jv, dtor, num_local_values, num_remote_values, krylov_dim);
  krylov->gmr = SpgmrMalloc(krylov_dim, krylov->z);
  krylov->max_restarts = max_restarts;
  krylov->gs = gs;
  newton_pc_vtable vtable = {.solve = gmres_solve,
                             .compute_p = krylov_compute_p,
                             .dtor = krylov_pc_dtor};
  return newton_pc_new("Krylov DAE preconditioner (GMRES)", krylov, vtable);
}

static bool bicgstab_solve(void* context, 
                           real_t t, real_t* x, real_t* x_dot,
                           real_t* r, real_t* z)
{
  krylov_pc_t* krylov = context;
  krylov->t = t;
  memcpy(NV_DATA(krylov->r), r, sizeof(real_t) * krylov->N_local);
  memcpy(NV_DATA(krylov->z), z, sizeof(real_t) * krylov->N_local);
  memcpy(NV_DATA(krylov->x), x, sizeof(real_t) * krylov->N_local);
  if (x_dot != NULL)
    memcpy(NV_DATA(krylov->x_dot), x, sizeof(real_t) * krylov->N_local);
  real_t res_norm;
  int nli, nps;
  int status = SpbcgSolve(krylov->bcg, krylov, krylov->z, krylov->r, PREC_NONE, krylov->delta, 
                          NULL, krylov->z_scaling, krylov->r_scaling, krylov->Atimes, NULL, 
                          &res_norm, &nli, &nps);
  switch(status)
  {
    case SPBCG_RES_REDUCED:
      log_debug("  Bi-CGSTAB: Residual norm was reduced (%g), but not below tolerance (%g)", res_norm, krylov->delta);
      break;
    case SPBCG_CONV_FAIL:
      log_debug("  Bi-CGSTAB: Failed to converge.");
      break;
    case SPBCG_ATIMES_FAIL_REC:
      log_debug("  Bi-CGSTAB: Jacobian-vector product failed recoverably.");
      break;
    case SPBCG_ATIMES_FAIL_UNREC:
      log_debug("  Bi-CGSTAB: Jacobian-vector product failed unrecoverably.");
    default:
      break;
  }
  return (status == SPBCG_SUCCESS);
}

newton_pc_t* bicgstab_newton_pc_new(MPI_Comm comm,
                                    void* context,
                                    int (*F)(void* context, real_t t, real_t* x, real_t* F),
                                    int (*Jv)(void* context, 
                                              real_t alpha, real_t beta, real_t gamma, 
                                              real_t t, real_t* x, real_t* v, real_t* Jv),
                                    void (*dtor)(void* context),
                                    int num_local_values, 
                                    int num_remote_values,
                                    int krylov_dim)
{
  krylov_pc_t* krylov = krylov_pc_new(comm, context, F, Jv, dtor, num_local_values, num_remote_values, krylov_dim);
  krylov->bcg = SpbcgMalloc(krylov_dim, krylov->z);
  newton_pc_vtable vtable = {.solve = bicgstab_solve,
                             .compute_p = krylov_compute_p,
                             .dtor = krylov_pc_dtor};
  return newton_pc_new("Krylov preconditioner (Bi-CGSTAB)", krylov, vtable);
}

newton_pc_t* dae_bicgstab_newton_pc_new(MPI_Comm comm,
                                        void* context,
                                        int (*F)(void* context, real_t t, real_t* x, real_t* x_dot, real_t* F),
                                        int (*Jv)(void* context, 
                                                  real_t alpha, real_t beta, real_t gamma, 
                                                  real_t t, real_t* x, real_t* x_dot, 
                                                  real_t* v, real_t* Jv),
                                        void (*dtor)(void* context),
                                        int num_local_values, 
                                        int num_remote_values,
                                        int krylov_dim)
{
  krylov_pc_t* krylov = dae_krylov_pc_new(comm, context, F, Jv, dtor, num_local_values, num_remote_values, krylov_dim);
  krylov->bcg = SpbcgMalloc(krylov_dim, krylov->z);
  newton_pc_vtable vtable = {.solve = bicgstab_solve,
                             .compute_p = krylov_compute_p,
                             .dtor = krylov_pc_dtor};
  return newton_pc_new("Krylov DAE preconditioner (Bi-CGSTAB)", krylov, vtable);
}

static bool tfqmr_solve(void* context, 
                        real_t t, real_t* x, real_t* x_dot,
                        real_t* r, real_t* z)
{
  krylov_pc_t* krylov = context;
  memcpy(NV_DATA(krylov->r), r, sizeof(real_t) * krylov->N_local);
  memcpy(NV_DATA(krylov->z), z, sizeof(real_t) * krylov->N_local);
  memcpy(NV_DATA(krylov->x), x, sizeof(real_t) * krylov->N_local);
  if (x_dot != NULL)
    memcpy(NV_DATA(krylov->x_dot), x, sizeof(real_t) * krylov->N_local);
  real_t res_norm;
  int nli, nps;
  int status = SptfqmrSolve(krylov->tfqmr, krylov, krylov->z, krylov->r, PREC_NONE, krylov->delta, 
                            NULL, krylov->z_scaling, krylov->r_scaling, krylov->Atimes, NULL, 
                            &res_norm, &nli, &nps);
  switch(status)
  {
    case SPTFQMR_RES_REDUCED:
      log_debug("  TFQMR: Residual norm was reduced (%g), but not below tolerance (%g)", res_norm, krylov->delta);
      break;
    case SPTFQMR_CONV_FAIL:
      log_debug("  TFQMR: Failed to converge.");
      break;
    case SPTFQMR_ATIMES_FAIL_REC:
      log_debug("  TFQMR: Jacobian-vector product failed recoverably.");
      break;
    case SPTFQMR_ATIMES_FAIL_UNREC:
      log_debug("  TFQMR: Jacobian-vector product failed unrecoverably.");
    default:
      break;
  }
  return (status == SPTFQMR_SUCCESS);
}

newton_pc_t* tfqmr_newton_pc_new(MPI_Comm comm, 
                                 void* context,
                                 int (*F)(void* context, real_t t, real_t* x, real_t* F),
                                 int (*Jv)(void* context, 
                                           real_t alpha, real_t beta, real_t gamma, 
                                           real_t t, real_t* x, real_t* v, real_t* Jv),
                                 void (*dtor)(void* context),
                                 int num_local_values, 
                                 int num_remote_values,
                                 int krylov_dim)
{
  krylov_pc_t* krylov = krylov_pc_new(comm, context, F, Jv, dtor, num_local_values, num_remote_values, krylov_dim);
  krylov->tfqmr = SptfqmrMalloc(krylov_dim, krylov->z);
  newton_pc_vtable vtable = {.solve = tfqmr_solve,
                             .compute_p = krylov_compute_p,
                             .dtor = krylov_pc_dtor};
  return newton_pc_new("Krylov preconditioner (TFQMR)", krylov, vtable);
}

newton_pc_t* dae_tfqmr_newton_pc_new(MPI_Comm comm,
                                     void* context,
                                     int (*F)(void* context, real_t t, real_t* x, real_t* x_dot, real_t* F),
                                     int (*Jv)(void* context, 
                                               real_t alpha, real_t beta, real_t gamma, 
                                               real_t t, real_t* x, real_t* x_dot, 
                                               real_t* v, real_t* Jv),
                                     void (*dtor)(void* context),
                                     int num_local_values, 
                                     int num_remote_values,
                                     int krylov_dim)
{
  krylov_pc_t* krylov = dae_krylov_pc_new(comm, context, F, Jv, dtor, num_local_values, num_remote_values, krylov_dim);
  krylov->tfqmr = SptfqmrMalloc(krylov_dim, krylov->z);
  newton_pc_vtable vtable = {.solve = tfqmr_solve,
                             .compute_p = krylov_compute_p,
                             .dtor = krylov_pc_dtor};
  return newton_pc_new("Krylov DAE preconditioner (TFQMR)", krylov, vtable);
}

bool newton_pc_is_krylov_newton_pc(newton_pc_t* pc)
{
  return string_contains(newton_pc_name(pc), "Krylov");
}

void krylov_newton_pc_set_scaling(newton_pc_t* pc, 
                                  real_t* z_scaling, 
                                  real_t* r_scaling)
{
  krylov_pc_t* krylov = newton_pc_context(pc);
  if (z_scaling != NULL)
  {
    if (krylov->z_scaling == NULL)
      krylov->z_scaling = N_VNew(krylov->comm, krylov->N_local);
    memcpy(NV_DATA(krylov->z_scaling), z_scaling, sizeof(krylov->N_total));
  }
  if (r_scaling != NULL)
  {
    if (krylov->r_scaling == NULL)
      krylov->r_scaling = N_VNew(krylov->comm, krylov->N_local);
    memcpy(NV_DATA(krylov->r_scaling), r_scaling, sizeof(krylov->N_total));
  }
}

void krylov_newton_pc_set_tolerance(newton_pc_t* pc,
                                    real_t res_norm)
{
  ASSERT(res_norm > 0.0);
  krylov_pc_t* krylov = newton_pc_context(pc);
  krylov->delta = res_norm;
}

int krylov_newton_pc_get_iterations(newton_pc_t* pc)
{
  krylov_pc_t* krylov = newton_pc_context(pc);
  return krylov->num_iters;
}

