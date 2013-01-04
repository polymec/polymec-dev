#include <gc/gc.h>
#include "core/polymec.h"
#include "core/linear_algebra.h"
#include "advect/slope_estimator.h"

#ifdef __cplusplus
extern "C" {
#endif

struct slope_estimator_t 
{
  double (*compute_factor)(double forward_slope, double backward_slope);
};

static double minmod(double forward_slope, double backward_slope)
{
  return 0.5 * (SIGN(forward_slope) + SIGN(backward_slope)) * MIN(fabs(forward_slope), fabs(backward_slope));
}

static double lax_wendroff(double forward_slope, double backward_slope)
{
  return 1.0;
}

static double beam_warming(double forward_slope, double backward_slope)
{
  return backward_slope / (forward_slope + 1e-15);
}

static double superbee(double forward_slope, double backward_slope)
{
  double theta = backward_slope / (forward_slope + 1e-15);
  double phi = MAX(0.0, MAX(MIN(1.0, 2.0*theta), MIN(theta, 2.0)));
  return forward_slope * phi;
}

static double van_leer(double forward_slope, double backward_slope)
{
  double theta = backward_slope / (forward_slope + 1e-15);
  double phi = (fabs(theta) + theta) / (1.0 + fabs(theta));
  return forward_slope * phi;
}

// This is a weighting function for our parabolic least-squares fit.
static inline double weight(double distance, double length_scale)
{
  return 1.0;
  double X = distance/length_scale;
  static double epsilon = 1e-4;
  return 1.0 / (X*X + epsilon*epsilon);
}

static inline void add_to_ls_system(double weight, double x, double val, double* A, double* b)
{
  A[0] += weight*1.0; A[3] += weight*x;     A[6] += weight*x*x;
  A[1] += weight*x;   A[4] += weight*x*x;   A[7] += weight*x*x*x;
  A[2] += weight*x*x; A[5] += weight*x*x*x; A[8] += weight*x*x*x*x;

  b[0] += weight*1.0*val;
  b[1] += weight*x*val;
  b[2] += weight*x*x*val;
}

static void compute_parabolic_fit(point_t* center_point,
                                  double center_value,
                                  point_t* forward_point,
                                  double forward_value,
                                  point_t* other_points,
                                  double* other_values,
                                  int num_other_points,
                                  double* coeffs)
{
  // Our length scale is the distance between the center point and the 
  // forward point.
  vector_t center_to_forward;
  point_displacement(center_point, forward_point, &center_to_forward);
  double h = vector_mag(&center_to_forward);
  ASSERT(h != 0.0);
//printf("xc = %g %g %g, xf = %g %g %g, h = %g\n", center_point->x, center_point->y, center_point->z, forward_point->x, forward_point->y, forward_point->z, h);
  vector_normalize(&center_to_forward); // Now a unit vector!

  // Consider the line on which lie the center point and forward point. 
  // We will weight the contributions of the points based on their 
  // distance from this line.
  double center_weight = weight(0.0, h);
  double forward_weight = center_weight;
  double other_weights[num_other_points];
  double projs[num_other_points];
  for (int i = 0; i < num_other_points; ++i)
  {
    // Find the distance from the other point to the line.
    vector_t other_to_center;
    point_displacement(&other_points[i], center_point, &other_to_center);
    vector_t cross_prod;
    vector_cross(&center_to_forward, &other_to_center, &cross_prod);
    double distance = vector_mag(&cross_prod); //  / vector_mag(&center_to_forward);

    // Project the other point to the line.
    projs[i] = -vector_dot(&other_to_center, &center_to_forward);
//if (projs[i] != 0.0)
//  printf("points[%d] = %g %g %g, proj = %g, val = %g\n", i, other_points[i].x, other_points[i].y, other_points[i].z, projs[i], other_values[i]);

    // Weight it accordingly.
    other_weights[i] = weight(distance, h);
//printf("weights[%d] = %g\n", i, other_weights[i]);
  }

  // Now we assemble the 3x3 moment matrix A and right-hand side vector b.
  double A[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
         b[3] = {0.0, 0.0, 0.0};

  // The center point only contributes to the constant term, since it's at 
  // the center!
  add_to_ls_system(center_weight, 0.0, center_value, A, b);

  // The forward point contributes at x = h.
  add_to_ls_system(forward_weight, h, forward_value, A, b);

  // Other contributions.
  for (int i = 0; i < num_other_points; ++i)
    add_to_ls_system(other_weights[i], projs[i], other_values[i], A, b);

//  printf("A = \n");
//  matrix_fprintf(A, 3, 3, stdout);
//  printf("b = \n");
//  vector_fprintf(b, 3, stdout);

  // Solve the system for the coefficients.
  solve_3x3(A, b, coeffs);
}

slope_estimator_t* slope_estimator_new(slope_limiter_t function)
{
  slope_estimator_t* estimator = GC_MALLOC(sizeof(slope_estimator_t));
  switch (function)
  {
    case (SLOPE_LIMITER_BEAM_WARMING):
      estimator->compute_factor = beam_warming;
      break;
    case (SLOPE_LIMITER_VAN_LEER):
      estimator->compute_factor = van_leer;
      break;
    case (SLOPE_LIMITER_SUPERBEE):
      estimator->compute_factor = superbee;
      break;
    case (SLOPE_LIMITER_MINMOD):
      estimator->compute_factor = minmod;
      break;
    case (SLOPE_LIMITER_LAX_WENDROFF):
      estimator->compute_factor = lax_wendroff;
      break;
  }
  return estimator;
}

double slope_estimator_value(slope_estimator_t* estimator, 
                             point_t* center_point,
                             double center_value,
                             point_t* forward_point,
                             double forward_value,
                             point_t* other_points,
                             double* other_values,
                             int num_other_points)
{
  // There must be at least one "other" point to proceed.
  ASSERT(num_other_points >= 1);

  // If we don't have a compute_factor function, just return 0.
  if (estimator->compute_factor == NULL)
    return 0.0;

  // Compute the coefficients for the parabolic fit.
  double coeffs[3];
  compute_parabolic_fit(center_point, center_value,
                        forward_point, forward_value,
                        other_points, other_values, num_other_points,
                        coeffs);
//printf("%g %g %g\n", coeffs[0], coeffs[1], coeffs[2]);

  // Compute the "backward" value of the solution using the fit.
  double L = point_distance(center_point, forward_point);
  double backward_value = coeffs[0] + coeffs[1]*(-L) + coeffs[2]*L*L;

  // Compute forward and backward slopes.
  double forward_slope = (forward_value - center_value) / L;
  double backward_slope = (center_value - backward_value) / L;
//printf("bval = %g, cval = %g, fval = %g\n", backward_value, center_value, forward_value);
//printf("fslope = %g, bslope = %g\n", forward_slope, backward_slope);

  // Now compute the estimator factor with the forward and backward slopes.
  return estimator->compute_factor(forward_slope, backward_slope);
}

#ifdef __cplusplus
}
#endif

