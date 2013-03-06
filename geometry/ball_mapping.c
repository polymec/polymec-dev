#include "geometry/ball_mapping.h"

#ifdef __cplusplus
extern "C" {
#endif

// We use the Gnomonic coordinates explained in
// Ronchi et al, JCP 124 p. 93 (1996).

typedef struct
{
  point_t x0;
  double r, l;
  int index;
} ball_mapping_t;

static void ball_map0(void* context, point_t* xi, double* x)
{
  // Inner cube.
  ball_mapping_t* ball = context;
  x[0] = ball->x0.x + ball->l * (xi->x - 0.5);
  x[1] = ball->x0.y + ball->l * (xi->y - 0.5);
  x[2] = ball->x0.z + ball->l * (xi->z - 0.5);
}

static void ball_map1(void* context, point_t* xi, double* x)
{
  // -x cube -- Region III in Ronchi et al.
  ball_mapping_t* ball = context;

  // Gnomonic and spherical coordinates for the surface of the ball.
  double Xi = M_PI/2 * (0.5 - xi->y);
  double Eta = M_PI/2 * (xi->z - 0.5);

  double phi = -M_PI + Xi;
  double Y = tan(Eta);
  double theta = -atan2(cos(phi), Y);
  double r0 = 0.5*ball->l * sqrt(1.0 + sin(Xi)*sin(Xi) + sin(Eta)*sin(Eta)); // Closest approach
  double r = ball->r - xi->x * (ball->r - r0);
  double sin_theta = sin(theta);

  x[0] = ball->x0.x + r * cos(phi) * sin_theta;
  x[1] = ball->x0.y + r * sin(phi) * sin_theta;
  x[2] = ball->x0.z + r * cos(theta);
}

static void ball_map2(void* context, point_t* xi, double* x)
{
  // +x cube -- Region I in Ronchi et al.
  ball_mapping_t* ball = context;

  // Gnomonic and spherical coordinates for the surface of the ball.
  double Xi = M_PI/2 * (xi->y - 0.5);
  double Eta = M_PI/2 * (xi->z - 0.5);

  double phi = Xi;
  double Y = tan(Eta);
  double theta = atan2(cos(phi), Y);
  double r0 = 0.5*ball->l * sqrt(1.0 + sin(Xi)*sin(Xi) + sin(Eta)*sin(Eta)); // Closest approach
  double r = r0 + xi->x * (ball->r - r0);
  double sin_theta = sin(theta);

  x[0] = ball->x0.x + r * cos(phi) * sin_theta;
  x[1] = ball->x0.y + r * sin(phi) * sin_theta;
  x[2] = ball->x0.z + r * cos(theta);
}

static void ball_map3(void* context, point_t* xi, double* x)
{
  // -y cube -- Region IV in Ronchi et al.
  ball_mapping_t* ball = context;

  // Gnomonic and spherical coordinates for the surface of the ball.
  double Xi = M_PI/2 * (xi->x - 0.5);
  double Eta = M_PI/2 * (xi->z - 0.5);

  double phi = -M_PI/2 + Xi;
  double Y = tan(Eta);
  double theta = -atan2(sin(phi), Y);
  double r0 = 0.5*ball->l * sqrt(1.0 + sin(Xi)*sin(Xi) + sin(Eta)*sin(Eta)); // Closest approach
  double r = ball->r - xi->y * (ball->r - r0);
  double sin_theta = sin(theta);

  x[0] = ball->x0.x + r * cos(phi) * sin_theta;
  x[1] = ball->x0.y + r * sin(phi) * sin_theta;
  x[2] = ball->x0.z + r * cos(theta);
}

static void ball_map4(void* context, point_t* xi, double* x)
{
  // +y cube -- Region II in Ronchi et al.
  ball_mapping_t* ball = context;

  // Gnomonic and spherical coordinates for the surface of the ball.
  double Xi = M_PI/2 * (0.5 - xi->x);
  double Eta = M_PI/2 * (xi->z - 0.5);

  double phi = M_PI/2 + Xi;
  double Y = tan(Eta);
  double theta = atan2(sin(phi), Y);
  double r0 = 0.5*ball->l * sqrt(1.0 + sin(Xi)*sin(Xi) + sin(Eta)*sin(Eta)); // Closest approach
  double r = r0 + xi->y * (ball->r - r0);
  double sin_theta = sin(theta);

  x[0] = ball->x0.x + r * cos(phi) * sin_theta;
  x[1] = ball->x0.y + r * sin(phi) * sin_theta;
  x[2] = ball->x0.z + r * cos(theta);
}

static void ball_map5(void* context, point_t* xi, double* x)
{
  // -z cube -- Region V in Ronchi et al.
  ball_mapping_t* ball = context;

  // Gnomonic and spherical coordinates for the surface of the ball.
  double Xi = M_PI/2 * (xi->x - 0.5);
  double Eta = M_PI/2 * (xi->y - 0.5);

  double X = tan(Xi);
  double Y = tan(Eta);
  double phi = atan2(Y, X);
  double theta = atan(X/cos(phi));
  double r0 = 0.5*ball->l * sqrt(1.0 + sin(Xi)*sin(Xi) + sin(Eta)*sin(Eta)); // Closest approach
  double r = r0 + xi->z * (ball->r - r0);
  double sin_theta = sin(theta);

  x[0] = ball->x0.x + r * cos(phi) * sin_theta;
  x[1] = ball->x0.y + r * sin(phi) * sin_theta;
  x[2] = ball->x0.z + r * cos(theta);
}

static void ball_map6(void* context, point_t* xi, double* x)
{
  // -z cube -- Region VI in Ronchi et al.
  ball_mapping_t* ball = context;

  // Gnomonic and spherical coordinates for the surface of the ball.
  double Xi = M_PI/2 * (xi->x - 0.5);
  double Eta = M_PI/2 * (xi->y - 0.5);

  double X = tan(Xi);
  double Y = tan(Eta);
  double phi = atan2(Y, X);
  double theta = M_PI - atan(X/cos(phi));
  double r0 = 0.5*ball->l * sqrt(1.0 + sin(Xi)*sin(Xi) + sin(Eta)*sin(Eta)); // Closest approach
  double r = ball->r - xi->z * (ball->r - r0);
  double sin_theta = sin(theta);

  x[0] = ball->x0.x + r * cos(phi) * sin_theta;
  x[1] = ball->x0.y + r * sin(phi) * sin_theta;
  x[2] = ball->x0.z + r * cos(theta);
}

static void ball_dtor(void* context)
{
  ball_mapping_t* ball = context;
  free(ball);
}

sp_func_t* ball_mapping_new(point_t* x0, double r, double l, int block_index)
{
  ASSERT(block_index >= 0);
  ASSERT(block_index <= 6);
  sp_vtable vtable = {.dtor = ball_dtor};
  switch (block_index)
  {
    case 0:
      vtable.eval = ball_map0;
      break;
    case 1:
      vtable.eval = ball_map1;
      break;
    case 2:
      vtable.eval = ball_map2;
      break;
    case 3:
      vtable.eval = ball_map3;
      break;
    case 4:
      vtable.eval = ball_map4;
      break;
    case 5:
      vtable.eval = ball_map5;
      break;
    case 6:
      vtable.eval = ball_map6;
      break;
  }
  char name[1024];
  snprintf(name, 1024, "Multi-block ball mapping (block %d)", block_index);
  ball_mapping_t* m = malloc(sizeof(ball_mapping_t));
  m->r = r;
  m->l = l;
  m->x0.x = x0->x;
  m->x0.y = x0->y;
  m->x0.z = x0->z;
  m->index = block_index;
  return sp_func_new(name, (void*)m, vtable, SP_INHOMOGENEOUS, 3);
}

#ifdef __cplusplus
}
#endif

