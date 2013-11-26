// Copyright 2012-2013 Jeffrey Johnson.
// 
// This file is part of Polymec, and is licensed under the Apache License, 
// Version 2.0 (the "License"); you may not use this file except in 
// compliance with the License. You may may find the text of the license in 
// the LICENSE file at the top-level source directory, or obtain a copy of 
// it at
// 
// http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include <gc/gc.h>
#include "core/polynomial.h"

struct polynomial_t 
{
  int order;
  double* coeffs;
  point_t x0;
};

static void polynomial_free(void* ctx, void* dummy)
{
  polynomial_t* p = ctx;
  free(p->coeffs);
}

static const int N_coeffs[5] = {1, 4, 10, 20, 35};

polynomial_t* polynomial_new(int order, double* coeffs, point_t* x0)
{
  ASSERT(order >= 0);
  ASSERT(order <= 4);
  polynomial_t* p = GC_MALLOC(sizeof(polynomial_t));
  p->order = order;
  p->coeffs = malloc(sizeof(double) * N_coeffs[order]);
  memset(p->coeffs, 0, sizeof(double) * N_coeffs[order]);
  p->x0 = *x0;
  GC_register_finalizer(p, polynomial_free, p, NULL, NULL);
  return p;
}

int polynomial_order(polynomial_t* p)
{
  return p->order;
}

int polynomial_num_coeffs(polynomial_t* p)
{
  return N_coeffs[p->order];
}

double* polynomial_coeffs(polynomial_t* p)
{
  return p->coeffs;
}

point_t* polynomial_x0(polynomial_t* p)
{
  return &p->x0;
}

static inline double poly_value0(polynomial_t* p, point_t* x)
{
  return p->coeffs[0];
}

static inline double poly1_value(polynomial_t* p, point_t* x)
{
  return poly_value0(p, x) + 
         p->coeffs[1]*x->x + p->coeffs[2]*x->y + p->coeffs[3]*x->z;
}

static inline double poly2_value(polynomial_t* p, point_t* X)
{
  double x = X->x, y = X->y, z = X->z;
  return poly1_value(p, X) + 
         p->coeffs[4]*x*x + p->coeffs[5]*x*y + p->coeffs[6]*x*z + 
         p->coeffs[7]*y*y + p->coeffs[8]*y*z + p->coeffs[9]*z*z;
}

static double poly3_value(polynomial_t* p, point_t* X)
{
  double x = X->x, y = X->y, z = X->z;
  return poly2_value(p, X) + 
         p->coeffs[10]*x*x*x + p->coeffs[11]*x*x*y + p->coeffs[12]*x*x*z + 
         p->coeffs[13]*x*y*y + p->coeffs[14]*x*y*z + p->coeffs[15]*x*z*z + 
         p->coeffs[16]*y*y*y + p->coeffs[17]*y*y*z + p->coeffs[18]*y*z*z + 
         p->coeffs[19]*z*z*z;

}

static double poly4_value(polynomial_t* p, point_t* X)
{
  double x = X->x, y = X->y, z = X->z;
  return poly3_value(p, X) + 
         p->coeffs[20]*x*x*x*x + p->coeffs[21]*x*x*x*y + p->coeffs[22]*x*x*x*z + 
         p->coeffs[23]*x*x*y*y + p->coeffs[24]*x*x*y*z + p->coeffs[25]*x*x*z*z + 
         p->coeffs[26]*x*y*y*y + p->coeffs[27]*x*y*y*z + p->coeffs[28]*x*y*z*z + 
         p->coeffs[29]*x*z*z*z + p->coeffs[30]*y*y*y*y + p->coeffs[31]*y*y*y*z + 
         p->coeffs[32]*y*y*z*z + p->coeffs[33]*y*z*z*z + p->coeffs[34]*z*z*z*z;
}

typedef double (*poly_value_func)(polynomial_t*, point_t*);
static poly_value_func poly_value[] = { poly_value0, poly1_value, poly2_value, poly3_value, poly4_value};

double polynomial_value(polynomial_t* p, point_t* x)
{
  return poly_value[p->order](p, x);
}

// Polynomial derivatives -- kind of a pain.
typedef double (*poly_deriv_func)(polynomial_t*, point_t*);

static double poly_zero(polynomial_t* p, point_t* x)
{
  return 0.0;
}

static double poly1_dx(polynomial_t* p, point_t* x)
{
  return p->coeffs[1];
}

static double poly1_dy(polynomial_t* p, point_t* x)
{
  return p->coeffs[2];
}

static double poly1_dz(polynomial_t* p, point_t* x)
{
  return p->coeffs[3];
}

static poly_deriv_func poly1_deriv[2][2][2] = { 
  { { poly1_value, poly1_dz },
    { poly1_dy, poly_zero } },
  { { poly1_dx, poly_zero },
    { poly_zero, poly_zero} } };

static double poly2_dx(polynomial_t* p, point_t* X)
{
  double x = X->x, y = X->y, z = X->z;
  return p->coeffs[1] + 2.0*p->coeffs[4]*x + p->coeffs[5]*y + p->coeffs[6]*z;
}

static double poly2_dy(polynomial_t* p, point_t* X)
{
  double x = X->x, y = X->y, z = X->z;
  return p->coeffs[2] + p->coeffs[5]*x + 2.0*p->coeffs[7]*y + p->coeffs[8]*z;
}

static double poly2_dz(polynomial_t* p, point_t* X)
{
  double x = X->x, y = X->y, z = X->z;
  return p->coeffs[3] + p->coeffs[6]*x + p->coeffs[8]*y + 2.0*p->coeffs[9]*z;
}

static double poly2_dxx(polynomial_t* p, point_t* X)
{
  return 2.0*p->coeffs[4];
}

static double poly2_dxy(polynomial_t* p, point_t* X)
{
  return p->coeffs[5];
}

static double poly2_dxz(polynomial_t* p, point_t* X)
{
  return p->coeffs[6];
}

static double poly2_dyy(polynomial_t* p, point_t* X)
{
  return 2.0*p->coeffs[7];
}

static double poly2_dyz(polynomial_t* p, point_t* X)
{
  return p->coeffs[8];
}

static double poly2_dzz(polynomial_t* p, point_t* X)
{
  return 2.0*p->coeffs[9];
}

static poly_deriv_func poly2_deriv[3][3][3] = { 
  { { poly2_value, poly2_dz, poly2_dzz },
    { poly2_dy, poly2_dyz, poly_zero },
    { poly2_dyy, poly_zero, poly_zero } },
  { { poly2_dx, poly2_dxz, poly_zero },
    { poly2_dxy, poly_zero, poly_zero},
    { poly_zero, poly_zero, poly_zero} },
  { { poly2_dxx, poly_zero, poly_zero },
    { poly_zero, poly_zero, poly_zero},
    { poly_zero, poly_zero, poly_zero} } };

static double poly3_dx(polynomial_t* p, point_t* X)
{
  double x = X->x, y = X->y, z = X->z;
  return p->coeffs[1] + 2.0*p->coeffs[4]*x + p->coeffs[5]*y + p->coeffs[6]*z +
         3.0*p->coeffs[10]*x*x + 2.0*p->coeffs[11]*x*y + 2.0*p->coeffs[12]*x*z +
         p->coeffs[13]*y*y + p->coeffs[14]*y*z + p->coeffs[15]*z*z;
}

static double poly3_dy(polynomial_t* p, point_t* X)
{
  double x = X->x, y = X->y, z = X->z;
  return p->coeffs[2] + p->coeffs[5]*x + 2.0*p->coeffs[7]*y + p->coeffs[8]*z +
         p->coeffs[11]*x*x + 2.0*p->coeffs[13]*x*y + p->coeffs[14]*x*z + 
         3.0*p->coeffs[16]*y*y + 2.0*p->coeffs[17]*y*z + p->coeffs[18]*z*z;
}

static double poly3_dz(polynomial_t* p, point_t* X)
{
  double x = X->x, y = X->y, z = X->z;
  return p->coeffs[3] + p->coeffs[6]*x + p->coeffs[8]*y + 2.0*p->coeffs[9]*z +
         p->coeffs[12]*x*x + p->coeffs[14]*x*y + 2.0*p->coeffs[15]*x*z + 
         p->coeffs[17]*y*y + 2.0*p->coeffs[18]*y*z + 3.0*p->coeffs[19]*z*z;
}

static double poly3_dxx(polynomial_t* p, point_t* X)
{
  double x = X->x, y = X->y, z = X->z;
  return 2.0*p->coeffs[4] + 6.0*p->coeffs[10]*x + 2.0*p->coeffs[11]*y + 2.0*p->coeffs[12]*z;
}

static double poly3_dxy(polynomial_t* p, point_t* X)
{
  double x = X->x, y = X->y, z = X->z;
  return p->coeffs[5] + 2.0*p->coeffs[11]*x + 
         2.0*p->coeffs[13]*y + p->coeffs[14]*z;
}

static double poly3_dxz(polynomial_t* p, point_t* X)
{
  double x = X->x, y = X->y, z = X->z;
  return p->coeffs[6] + 2.0*p->coeffs[12]*x +
         p->coeffs[14]*y + 2.0*p->coeffs[15]*z;
}

static double poly3_dyy(polynomial_t* p, point_t* X)
{
  double x = X->x, y = X->y, z = X->z;
  return 2.0*p->coeffs[7] + 2.0*p->coeffs[13]*x +
         6.0*p->coeffs[16]*y + 2.0*p->coeffs[17]*z;
}

static double poly3_dyz(polynomial_t* p, point_t* X)
{
  double x = X->x, y = X->y, z = X->z;
  return p->coeffs[8] + p->coeffs[14]*x + 
         2.0*p->coeffs[17]*y + 2.0*p->coeffs[18]*z;
}

static double poly3_dzz(polynomial_t* p, point_t* X)
{
  double y = X->y, z = X->z;
  return 2.0*p->coeffs[9] + 2.0*p->coeffs[15]*z + 
         2.0*p->coeffs[18]*y + 6.0*p->coeffs[19]*z;
}

static double poly3_dxxx(polynomial_t* p, point_t* X)
{
  return 6.0*p->coeffs[10];
}

static double poly3_dxxy(polynomial_t* p, point_t* X)
{
  return 2.0*p->coeffs[11];
}

static double poly3_dxxz(polynomial_t* p, point_t* X)
{
  return 2.0*p->coeffs[12];
}

static double poly3_dxyy(polynomial_t* p, point_t* X)
{
  return 2.0*p->coeffs[13];
}

static double poly3_dxyz(polynomial_t* p, point_t* X)
{
  return p->coeffs[14];
}

static double poly3_dxzz(polynomial_t* p, point_t* X)
{
  return 2.0*p->coeffs[15];
}

static double poly3_dyyy(polynomial_t* p, point_t* X)
{
  return 6.0*p->coeffs[16];
}

static double poly3_dyyz(polynomial_t* p, point_t* X)
{
  return 2.0*p->coeffs[17];
}

static double poly3_dyzz(polynomial_t* p, point_t* X)
{
  return 2.0*p->coeffs[18];
}

static double poly3_dzzz(polynomial_t* p, point_t* X)
{
  return 6.0*p->coeffs[19];
}

static poly_deriv_func poly3_deriv[4][4][4] = { 
  { { poly3_value, poly3_dz, poly3_dzz, poly3_dzzz },
    { poly3_dy, poly3_dyz, poly3_dyzz, poly_zero },
    { poly3_dyy, poly3_dyyz, poly_zero, poly_zero },
    { poly3_dyyy, poly_zero, poly_zero, poly_zero} },
  { { poly3_dx, poly3_dxz, poly3_dxzz, poly_zero },
    { poly3_dxy, poly3_dxyz, poly_zero, poly_zero},
    { poly3_dxyy, poly_zero, poly_zero, poly_zero},
    { poly_zero, poly_zero, poly_zero, poly_zero} },
  { { poly3_dxx, poly3_dxxz, poly_zero, poly_zero},
    { poly3_dxxy, poly_zero, poly_zero, poly_zero},
    { poly_zero, poly_zero, poly_zero, poly_zero},
    { poly_zero, poly_zero, poly_zero, poly_zero} },
  { { poly3_dxxx, poly_zero, poly_zero, poly_zero },
    { poly_zero, poly_zero, poly_zero, poly_zero},
    { poly_zero, poly_zero, poly_zero, poly_zero},
    { poly_zero, poly_zero, poly_zero, poly_zero} } };

static double poly4_dx(polynomial_t* p, point_t* X)
{
  double x = X->x, y = X->y, z = X->z;
  return poly3_dx(p, X) + 
    4.0*p->coeffs[20]*x*x*x + 3.0*p->coeffs[21]*x*x*y + 3.0*p->coeffs[22]*x*x*z +
    2.0*p->coeffs[23]*x*y*y + 2.0*p->coeffs[24]*x*y*z + 2.0*p->coeffs[25]*x*z*z + 
    p->coeffs[26]*y*y*y + p->coeffs[27]*y*y*z + p->coeffs[28]*y*z*z + p->coeffs[29]*z*z*z;
}

static double poly4_dy(polynomial_t* p, point_t* X)
{
  double x = X->x, y = X->y, z = X->z;
  return poly3_dy(p, X) + 
    p->coeffs[21]*x*x*x + 2.0*p->coeffs[23]*x*x*y + p->coeffs[24]*x*x*z + 
    3.0*p->coeffs[26]*x*y*y + 2.0*p->coeffs[27]*x*y*z + p->coeffs[28]*x*z*z + 
    4.0*p->coeffs[30]*y*y*y + 3.0*p->coeffs[31]*y*y*z + 2.0*p->coeffs[32]*y*z*z + 
    p->coeffs[33]*z*z*z;
}

static double poly4_dz(polynomial_t* p, point_t* X)
{
  double x = X->x, y = X->y, z = X->z;
  return poly3_dz(p, X) + 
    p->coeffs[22]*x*x*x + p->coeffs[24]*x*x*y + 2.0*p->coeffs[25]*x*x*z + 
    p->coeffs[27]*x*y*y + 2.0*p->coeffs[28]*x*y*z + 3.0*p->coeffs[29]*x*z*z + 
    p->coeffs[31]*y*y*y + 2.0*p->coeffs[32]*y*y*z + 3.0*p->coeffs[33]*y*z*z + 
    4.0*p->coeffs[34]*z*z*z;
}

static double poly4_dxx(polynomial_t* p, point_t* X)
{
  double x = X->x, y = X->y, z = X->z;
  return poly3_dxx(p, X) + 
    12.0*p->coeffs[20]*x*x + 6.0*p->coeffs[21]*x*y + 6.0*p->coeffs[22]*x*z +
    2.0*p->coeffs[23]*y*y + 2.0*p->coeffs[24]*y*z + 2.0*p->coeffs[25]*z*z;
}

static double poly4_dxy(polynomial_t* p, point_t* X)
{
  double x = X->x, y = X->y, z = X->z;
  return poly3_dxy(p, X) + 
    3.0*p->coeffs[21]*x*x + 4.0*p->coeffs[23]*x*y + 2.0*p->coeffs[24]*x*z + 
    3.0*p->coeffs[26]*y*y + 2.0*p->coeffs[27]*y*z + p->coeffs[28]*z*z;
}

static double poly4_dxz(polynomial_t* p, point_t* X)
{
  double x = X->x, y = X->y, z = X->z;
  return poly3_dxz(p, X) + 
    3.0*p->coeffs[22]*x*x + 2.0*p->coeffs[24]*x*y + 4.0*p->coeffs[25]*x*z + 
    p->coeffs[27]*y*y + 2.0*p->coeffs[28]*y*z + 3.0*p->coeffs[29]*z*z;
}

static double poly4_dyy(polynomial_t* p, point_t* X)
{
  double x = X->x, y = X->y, z = X->z;
  return poly3_dyy(p, X) + 
    2.0*p->coeffs[23]*x*x + 6.0*p->coeffs[26]*x*y + 2.0*p->coeffs[27]*x*z + 
    12.0*p->coeffs[30]*y*y + 6.0*p->coeffs[31]*y*z + 2.0*p->coeffs[32]*z*z;
}

static double poly4_dyz(polynomial_t* p, point_t* X)
{
  double x = X->x, y = X->y, z = X->z;
  return poly3_dyz(p, X) + 
    p->coeffs[24]*x*x + 2.0*p->coeffs[27]*x*y + 2.0*p->coeffs[28]*x*z + 
    3.0*p->coeffs[31]*y*y + 4.0*p->coeffs[32]*y*z + 3.0*p->coeffs[33]*z*z;
}

static double poly4_dzz(polynomial_t* p, point_t* X)
{
  double x = X->x, y = X->y, z = X->z;
  return poly3_dzz(p, X) + 
    2.0*p->coeffs[25]*x*x + 2.0*p->coeffs[28]*x*y + 6.0*p->coeffs[29]*x*z + 
    2.0*p->coeffs[32]*y*y + 6.0*p->coeffs[33]*y*z + 12.0*p->coeffs[34]*z*z;
}

static double poly4_dxxx(polynomial_t* p, point_t* X)
{
  double x = X->x, y = X->y, z = X->z;
  return poly3_dxxx(p, X) + 
    24.0*p->coeffs[20]*x + 6.0*p->coeffs[21]*y + 6.0*p->coeffs[22]*z;
}

static double poly4_dxxy(polynomial_t* p, point_t* X)
{
  double x = X->x, y = X->y, z = X->z;
  return poly3_dxxy(p, X) + 
    6.0*p->coeffs[21]*x + 4.0*p->coeffs[23]*y + 2.0*p->coeffs[24]*z;
}

static double poly4_dxxz(polynomial_t* p, point_t* X)
{
  double x = X->x, y = X->y, z = X->z;
  return poly3_dxxz(p, X) + 
    6.0*p->coeffs[22]*x + 2.0*p->coeffs[24]*y + 4.0*p->coeffs[25]*z;
}

static double poly4_dxyy(polynomial_t* p, point_t* X)
{
  double x = X->x, y = X->y, z = X->z;
  return poly3_dxyy(p, X) + 
    4.0*p->coeffs[23]*x + 6.0*p->coeffs[26]*y + 2.0*p->coeffs[27]*z;
}

static double poly4_dxyz(polynomial_t* p, point_t* X)
{
  double x = X->x, y = X->y, z = X->z;
  return poly3_dxyz(p, X) + 
    2.0*p->coeffs[24]*x + 2.0*p->coeffs[27]*y + 2.0*p->coeffs[28]*z;
}

static double poly4_dxzz(polynomial_t* p, point_t* X)
{
  double x = X->x, y = X->y, z = X->z;
  return poly3_dxzz(p, X) + 
    4.0*p->coeffs[25]*x + 2.0*p->coeffs[28]*y + 6.0*p->coeffs[29]*z;
}

static double poly4_dyyy(polynomial_t* p, point_t* X)
{
  double x = X->x, y = X->y, z = X->z;
  return poly3_dyyy(p, X) + 
    6.0*p->coeffs[26]*x + 24.0*p->coeffs[30]*y + 6.0*p->coeffs[31]*z;
}

static double poly4_dyyz(polynomial_t* p, point_t* X)
{
  double x = X->x, y = X->y, z = X->z;
  return poly3_dyyz(p, X) + 
    2.0*p->coeffs[27]*x + 6.0*p->coeffs[31]*y + 4.0*p->coeffs[32]*z;
}

static double poly4_dyzz(polynomial_t* p, point_t* X)
{
  double x = X->x, y = X->y, z = X->z;
  return poly3_dyzz(p, X) + 
    2.0*p->coeffs[28]*x + 4.0*p->coeffs[32]*y + 6.0*p->coeffs[33]*z;
}

static double poly4_dzzz(polynomial_t* p, point_t* X)
{
  double x = X->x, y = X->y, z = X->z;
  return poly3_dzzz(p, X) + 
    6.0*p->coeffs[29]*x + 6.0*p->coeffs[33]*y + 24.0*p->coeffs[34]*z;
}

static double poly4_dxxxx(polynomial_t* p, point_t* X)
{
  return 24.0*p->coeffs[20];
}

static double poly4_dxxxy(polynomial_t* p, point_t* X)
{
  return 6.0*p->coeffs[21];
}

static double poly4_dxxxz(polynomial_t* p, point_t* X)
{
  return 6.0*p->coeffs[22];
}

static double poly4_dxxyy(polynomial_t* p, point_t* X)
{
  return 4.0*p->coeffs[23];
}

static double poly4_dxxyz(polynomial_t* p, point_t* X)
{
  return 2.0*p->coeffs[24];
}

static double poly4_dxxzz(polynomial_t* p, point_t* X)
{
  return 4.0*p->coeffs[25];
}

static double poly4_dxyyy(polynomial_t* p, point_t* X)
{
  return 6.0*p->coeffs[26];
}

static double poly4_dxyyz(polynomial_t* p, point_t* X)
{
  return 2.0*p->coeffs[27];
}

static double poly4_dxyzz(polynomial_t* p, point_t* X)
{
  return 2.0*p->coeffs[28];
}

static double poly4_dxzzz(polynomial_t* p, point_t* X)
{
  return 6.0*p->coeffs[29];
}

static double poly4_dyyyy(polynomial_t* p, point_t* X)
{
  return 24.0*p->coeffs[30];
}

static double poly4_dyyyz(polynomial_t* p, point_t* X)
{
  return 6.0*p->coeffs[31];
}

static double poly4_dyyzz(polynomial_t* p, point_t* X)
{
  return 4.0*p->coeffs[32];
}

static double poly4_dyzzz(polynomial_t* p, point_t* X)
{
  return 6.0*p->coeffs[33];
}

static double poly4_dzzzz(polynomial_t* p, point_t* X)
{
  return 24.0*p->coeffs[34];
}

static poly_deriv_func poly4_deriv[5][5][5] = { 
  { { poly4_value, poly4_dz, poly4_dzz, poly4_dzzz, poly4_dzzzz },
    { poly4_dy, poly4_dyz, poly4_dyzz, poly4_dyzzz, poly_zero },
    { poly4_dyy, poly4_dyyz, poly4_dyyzz, poly_zero, poly_zero },
    { poly4_dyyy, poly4_dyyyz, poly_zero, poly_zero, poly_zero},
    { poly4_dyyyy, poly_zero, poly_zero, poly_zero, poly_zero} },
  { { poly4_dx, poly4_dxz, poly4_dxzz, poly4_dxzzz, poly_zero },
    { poly4_dxy, poly4_dxyz, poly4_dxyzz, poly_zero, poly_zero},
    { poly4_dxyy, poly4_dxyyz, poly_zero, poly_zero, poly_zero},
    { poly4_dxyyy, poly_zero, poly_zero, poly_zero, poly_zero},
    { poly_zero, poly_zero, poly_zero, poly_zero, poly_zero} },
  { { poly4_dxx, poly4_dxxz, poly4_dxxzz, poly_zero, poly_zero},
    { poly4_dxxy, poly4_dxxyz, poly_zero, poly_zero, poly_zero},
    { poly4_dxxyy, poly_zero, poly_zero, poly_zero, poly_zero},
    { poly_zero, poly_zero, poly_zero, poly_zero, poly_zero},
    { poly_zero, poly_zero, poly_zero, poly_zero, poly_zero} },
  { { poly4_dxxx, poly4_dxxxz, poly_zero, poly_zero, poly_zero },
    { poly4_dxxxy, poly_zero, poly_zero, poly_zero, poly_zero},
    { poly_zero, poly_zero, poly_zero, poly_zero, poly_zero},
    { poly_zero, poly_zero, poly_zero, poly_zero, poly_zero},
    { poly_zero, poly_zero, poly_zero, poly_zero, poly_zero} },
  { { poly4_dxxxx, poly_zero, poly_zero, poly_zero, poly_zero },
    { poly_zero, poly_zero, poly_zero, poly_zero, poly_zero},
    { poly_zero, poly_zero, poly_zero, poly_zero, poly_zero},
    { poly_zero, poly_zero, poly_zero, poly_zero, poly_zero},
    { poly_zero, poly_zero, poly_zero, poly_zero, poly_zero} } };

double polynomial_deriv(polynomial_t* p, int x_deriv, int y_deriv, int z_deriv, point_t* x)
{
  if ((x_deriv + y_deriv + z_deriv) > p->order)
    return 0.0;
  else if (p->order == 1)
    return poly1_deriv[x_deriv][y_deriv][z_deriv](p, x);
  else if (p->order == 2)
    return poly2_deriv[x_deriv][y_deriv][z_deriv](p, x);
  else if (p->order == 3)
    return poly3_deriv[x_deriv][y_deriv][z_deriv](p, x);
  else 
  {
    ASSERT(p->order == 4);
    return poly4_deriv[x_deriv][y_deriv][z_deriv](p, x);
  }
}

static int x_pow[35] = {0, 1, 0, 0, 2, 1, 1, 0, 0, 0, 
                        3, 2, 2, 1, 1, 1, 0, 0, 0, 0, 
                        4, 3, 3, 2, 2, 2, 1, 1, 1, 1, 
                        0, 0, 0, 0, 0};
static int y_pow[35] = {0, 0, 1, 0, 0, 1, 0, 2, 1, 0, 
                        0, 1, 0, 2, 1, 0, 3, 2, 1, 0, 
                        0, 1, 0, 2, 1, 0, 3, 2, 1, 0, 
                        4, 3, 2, 1, 0};
static int z_pow[35] = {0, 0, 0, 1, 0, 0, 1, 0, 1, 2, 
                        0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 
                        0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 
                        0, 1, 2, 3, 4};

bool polynomial_next(polynomial_t* p, int* pos, double* coeff, int* x_power, int* y_power, int* z_power)
{
  if (*pos >= N_coeffs[p->order])
    return false;
  *coeff = p->coeffs[*pos];
  *x_power = x_pow[*pos];
  *y_power = y_pow[*pos];
  *z_power = z_pow[*pos];
  ++(*pos);
  return true;
}

static void wrap_eval(void* context, point_t* x, double* result)
{
  polynomial_t* p = context;
  *result = polynomial_value(p, x);
}

static void wrap_eval_deriv(void* context, int deriv, point_t* x, double* result)
{
  polynomial_t* p = context;
  int result_size = pow(3, deriv);
  if (deriv > p->order)
    memset(result, 0, sizeof(double) * result_size);
  else if (deriv == 0)
    *result = polynomial_value(p, x);
  else if (deriv == 1)
  {
    result[0] = polynomial_deriv(p, 1, 0, 0, x);
    result[1] = polynomial_deriv(p, 0, 1, 0, x);
    result[2] = polynomial_deriv(p, 0, 0, 1, x);
  }
  else 
  {
    // FIXME
    POLYMEC_NOT_IMPLEMENTED
  }
}

static bool wrap_has_deriv(void* context, int deriv)
{
  // Polynomials are analytic.
  return true;
}

sp_func_t* polynomial_sp_func(polynomial_t* p)
{
  sp_vtable vtable = {.eval = wrap_eval,
                      .eval_deriv = wrap_eval_deriv,
                      .has_deriv = wrap_has_deriv,
                      .dtor = NULL};
  char name[128];
  snprintf(name, 128, "polynomial (p = %d)", p->order);
  return sp_func_new(name, p, vtable, SP_INHOMOGENEOUS, 1);
}

