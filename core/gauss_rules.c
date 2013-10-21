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

#include "core/polymec.h"
#include "core/gauss_rules.h"

void get_gauss_points(int n, double* points, double* weights)
{
  ASSERT(n >= 0);
  ASSERT(points != NULL);
  ASSERT(weights != NULL);

  int m = (n+1)/2;
  double p1, p2, p3;
  double pp, z;
  for (int i = 1; i <= m; ++i)
  {
    z = cos(M_PI * (i - 0.25) / (n + 0.5));

    while(1)
    {
      p1 = 1;
      p2 = 0;
      for (int j = 1; j <= n; j++)
      {
        p3 = p2;
        p2 = p1;
        p1 = ((2 * j - 1) * z * p2 - (j - 1) * p3) / j;
      }
      // p1 is Legendre polynomial

      pp = n * (z*p1-p2) / (z*z - 1);

      if (fabs(p1/pp) < 2e-16) break;

      z = z - p1/pp;
    }

    z = ((1 - z) + p1/pp)/2;

    points[i-1] = z;
    points[n-i] = 1 - z;
    weights[i-1] = 1./(4*z*(1 - z)*pp*pp);
  }
}

void get_gauss_legendre_points(int n, double* points, double* weights)
{
  ASSERT(n >= 0);
  ASSERT(points != NULL);
  ASSERT(weights != NULL);
}

void get_gauss_radau_points(int n, double* points, double* weights)
{
  ASSERT(n >= 0);
  ASSERT(points != NULL);
  ASSERT(weights != NULL);
}

void get_gauss_lobatto_points(int n, double* points, double* weights)
{
  ASSERT(n >= 0);
  ASSERT(points != NULL);
  ASSERT(weights != NULL);

  if (n == 0)
  {
    points[0] = 0.5;
    weights[0] = 0.5;
  }
  else
  {
    points[0] = 0.;
    weights[0] = 0.;
    points[n] = 1.;
    weights[n] = 1.;
    if (n == 1) return;

    int m = (n - 1)/2, odd_n = n%2;

    if (!odd_n)
    {
      points[m+1] = 0.5;
      weights[m+1] = 1.0;
    }
    for (int i = 0; i < m; )
    {
      double y, z, d0, s0;
      z = cos(M_PI*(i + 1)/n);

      int k = 0;
      while (true)
      {
        // compute d0, s0 -- P'_p(z), P"_p(z)
        // (n+1)*P_{n+1}(z) = (2*n+1)*z*P_n(z)-n*P_{n-1}(z)
        // P'_{n+1}(z) = (2*n+1)* P_n(z)+P'_{n-1}(z)
        // P"_{n+1}(z) = (2*n+1)*P'_n(z)+P"_{n-1}(z)
        {
          double p0, p1, p2, d1;
          p2 = 1.;
          p1 = z;
          d0 = odd_n;
          d1 = 1 - odd_n;
          s0 = 0;
          for (int nn = 1; true; nn++)
          {
            p0 = ((2*nn+1)*z*p1 - nn*p2)/(nn + 1);
            if (nn%2 != odd_n)
            {
              d0 += (2*nn+1)*p1;
              s0 += (2*nn+1)*d1;
            }
            else
            {
              d1 += (2*nn+1)*p1;
            }
            if (nn == n - 1) break;
            p2 = p1;
            p1 = p0;
          }
        }

        if (fabs(d0/s0) < 2e-16) break;
        ++k;
        ASSERT(k < 6);

        z = z - d0/s0;
      }

      y = ((1 - z) + d0/s0)/2;

      points[++i] = y;
      weights[++i] = y;
      points[n-i] = 1 - y;
      weights[n-i] = y;
    }
  }
}

