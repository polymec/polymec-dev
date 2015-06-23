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
  void* F_context;
  SpgmrMemRec* gmr;
  SpbcgMemRec* bcg;
  SptfqmrMemRec* tfqmr;
  int N_local, N_remote, N_total;
  int (*F)(void* context, real_t t, real_t* x, real_t* F);
  int (*dae_F)(void* context, real_t t, real_t* x, real_t* xdot, real_t* F);
  int (*Jv)(void* context, real_t t, real_t* x, real_t* v, real_t* Jv);
  int (*dae_Jv)(void* context, real_t t, real_t* x, real_t* xdot, real_t* v, real_t* Jv);
  int (*Atimes)(void* A_data, N_Vector v, N_Vector Av);
  void (*dtor)(void* context);
  int krylov_dim;
  N_Vector z, r, z_scaling, r_scaling;
  int max_restarts;
  gmres_gram_schmidt_t gs;
  real_t delta;
  int num_iters;
} krylov_pc_t;

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

static int dae_Jv_adaptor(void* context, real_t t, real_t* x, real_t* x_dot, real_t* v, real_t* Jv)
{
  // We are passed the actual preconditioner as our context pointer, so get the 
  // "real" one here.
  krylov_pc_t* krylov = context;
  ASSERT(krylov->Jv != NULL);
  return krylov->Jv(krylov->context, t, x, v, Jv);
}

static int Atimes_DQ_adaptor(void* A_data, N_Vector v, N_Vector Av)
{
  return 0; // FIXME
}

static int Atimes_Jv_adaptor(void* A_data, N_Vector v, N_Vector Av)
{
  return 0; // FIXME
}

static krylov_pc_t* krylov_pc_new(MPI_Comm comm,
                                  void* context,
                                  int (*F)(void* context, real_t t, real_t* x, real_t* F),
                                  int (*Jv)(void* context, real_t t, real_t* x, real_t* v, real_t* Jv),
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
  krylov->F_context = context;
  krylov->gmr = NULL;
  krylov->bcg = NULL;
  krylov->tfqmr = NULL;
  krylov->F = F;
  krylov->Jv = Jv;
  krylov->dae_F = (F != NULL) ? dae_F_adaptor : NULL;
  krylov->dae_Jv = (Jv != NULL) ? dae_Jv_adaptor : NULL;
  krylov->Atimes = (F != NULL) ? Atimes_DQ_adaptor : Atimes_Jv_adaptor;
  krylov->dtor = dtor;
  krylov->N_local = num_local_values;
  krylov->N_remote = num_remote_values;
  krylov->N_total = krylov->N_local + krylov->N_remote;
  krylov->z = N_VNew(comm, krylov->N_local);
  krylov->r = N_VNew(comm, krylov->N_local);
  krylov->z_scaling = NULL;
  krylov->r_scaling = NULL;
  return krylov;
}

static krylov_pc_t* dae_krylov_pc_new(MPI_Comm comm,
                                      void* context,
                                      int (*F)(void* context, real_t t, real_t* x, real_t* x_dot, real_t* F),
                                      int (*Jv)(void* context, real_t t, real_t* x, real_t* x_dot, real_t* v, real_t* Jv),
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
  krylov->F_context = krylov;
  krylov->gmr = NULL;
  krylov->bcg = NULL;
  krylov->tfqmr = NULL;
  krylov->F = NULL;
  krylov->Jv = NULL;
  krylov->dae_F = F;
  krylov->dae_Jv = Jv;
  krylov->Atimes = (F != NULL) ? Atimes_DQ_adaptor : Atimes_Jv_adaptor;
  krylov->dtor = dtor;
  krylov->N_remote = num_remote_values;
  krylov->N_total = krylov->N_local + krylov->N_remote;
  krylov->z = N_VNew(comm, krylov->N_local);
  krylov->r = N_VNew(comm, krylov->N_local);
  krylov->z_scaling = NULL;
  krylov->r_scaling = NULL;
  return krylov;
}

static bool gmres_solve(void* context, 
                        real_t t, real_t* x, real_t* xdot,
                        real_t* r, real_t* z)
{
  krylov_pc_t* krylov = context;
  real_t res_norm;
  int nli, nps;
  int status = SpgmrSolve(krylov->gmr, krylov->F_context, krylov->z, krylov->r, PREC_NONE, krylov->gs,
                          krylov->delta, krylov->max_restarts, NULL, krylov->r_scaling, krylov->z_scaling,
                          krylov->Atimes, NULL, &res_norm, &nli, &nps);
  return (status == SPGMR_SUCCESS);
}

newton_pc_t* gmres_newton_pc_new(MPI_Comm comm,
                                 void* context,
                                 int (*F)(void* context, real_t t, real_t* x, real_t* F),
                                 int (*Jv)(void* context, real_t t, real_t* x, real_t* v, real_t* Jv),
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
                             .dtor = krylov_pc_dtor};
  return newton_pc_new("Krylov preconditioner (GMRES)", krylov, vtable);
}
                                        
newton_pc_t* dae_gmres_newton_pc_new(MPI_Comm comm,
                                     void* context,
                                     int (*F)(void* context, real_t t, real_t* x, real_t* x_dot, real_t* F),
                                     int (*Jv)(void* context, real_t t, real_t* x, real_t* x_dot, real_t* v, real_t* Jv),
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
                             .dtor = krylov_pc_dtor};
  return newton_pc_new("Krylov DAE preconditioner (GMRES)", krylov, vtable);
}

static bool bicgstab_solve(void* context, 
                           real_t t, real_t* x, real_t* xdot,
                           real_t* r, real_t* z)
{
  krylov_pc_t* krylov = context;
  real_t res_norm;
  int nli, nps;
  int status = SpbcgSolve(krylov->bcg, krylov->F_context, krylov->z, krylov->r, PREC_NONE, krylov->delta, 
                          NULL, krylov->z_scaling, krylov->r_scaling, krylov->Atimes, NULL, 
                          &res_norm, &nli, &nps);
  return (status == SPBCG_SUCCESS);
}

newton_pc_t* bicgstab_newton_pc_new(MPI_Comm comm,
                                    void* context,
                                    int (*F)(void* context, real_t t, real_t* x, real_t* F),
                                    int (*Jv)(void* context, real_t t, real_t* x, real_t* v, real_t* Jv),
                                    void (*dtor)(void* context),
                                    int num_local_values, 
                                    int num_remote_values,
                                    int krylov_dim)
{
  krylov_pc_t* krylov = krylov_pc_new(comm, context, F, Jv, dtor, num_local_values, num_remote_values, krylov_dim);
  krylov->bcg = SpbcgMalloc(krylov_dim, krylov->z);
  newton_pc_vtable vtable = {.solve = bicgstab_solve,
                             .dtor = krylov_pc_dtor};
  return newton_pc_new("Krylov preconditioner (Bi-CGSTAB)", krylov, vtable);
}

newton_pc_t* dae_bicgstab_newton_pc_new(MPI_Comm comm,
                                        void* context,
                                        int (*F)(void* context, real_t t, real_t* x, real_t* x_dot, real_t* F),
                                        int (*Jv)(void* context, real_t t, real_t* x, real_t* x_dot, real_t* v, real_t* Jv),
                                        void (*dtor)(void* context),
                                        int num_local_values, 
                                        int num_remote_values,
                                        int krylov_dim)
{
  krylov_pc_t* krylov = dae_krylov_pc_new(comm, context, F, Jv, dtor, num_local_values, num_remote_values, krylov_dim);
  krylov->bcg = SpbcgMalloc(krylov_dim, krylov->z);
  newton_pc_vtable vtable = {.solve = bicgstab_solve,
                             .dtor = krylov_pc_dtor};
  return newton_pc_new("Krylov DAE preconditioner (Bi-CGSTAB)", krylov, vtable);
}

static bool tfqmr_solve(void* context, 
                        real_t t, real_t* x, real_t* xdot,
                        real_t* r, real_t* z)
{
  krylov_pc_t* krylov = context;
  real_t res_norm;
  int nli, nps;
  int status = SptfqmrSolve(krylov->tfqmr, krylov->F_context, krylov->z, krylov->r, PREC_NONE, krylov->delta, 
                            NULL, krylov->z_scaling, krylov->r_scaling, krylov->Atimes, NULL, 
                            &res_norm, &nli, &nps);
  return (status == SPTFQMR_SUCCESS);
}

newton_pc_t* tfqmr_newton_pc_new(MPI_Comm comm, 
                                 void* context,
                                 int (*F)(void* context, real_t t, real_t* x, real_t* F),
                                 int (*Jv)(void* context, real_t t, real_t* x, real_t* v, real_t* Jv),
                                 void (*dtor)(void* context),
                                 int num_local_values, 
                                 int num_remote_values,
                                 int krylov_dim)
{
  krylov_pc_t* krylov = krylov_pc_new(comm, context, F, Jv, dtor, num_local_values, num_remote_values, krylov_dim);
  krylov->tfqmr = SptfqmrMalloc(krylov_dim, krylov->z);
  newton_pc_vtable vtable = {.solve = tfqmr_solve,
                             .dtor = krylov_pc_dtor};
  return newton_pc_new("Krylov preconditioner (TFQMR)", krylov, vtable);
}

newton_pc_t* dae_tfqmr_newton_pc_new(MPI_Comm comm,
                                     void* context,
                                     int (*F)(void* context, real_t t, real_t* x, real_t* x_dot, real_t* F),
                                     int (*Jv)(void* context, real_t t, real_t* x, real_t* x_dot, real_t* v, real_t* Jv),
                                     void (*dtor)(void* context),
                                     int num_local_values, 
                                     int num_remote_values,
                                     int krylov_dim)
{
  krylov_pc_t* krylov = dae_krylov_pc_new(comm, context, F, Jv, dtor, num_local_values, num_remote_values, krylov_dim);
  krylov->tfqmr = SptfqmrMalloc(krylov_dim, krylov->z);
  newton_pc_vtable vtable = {.solve = tfqmr_solve,
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

