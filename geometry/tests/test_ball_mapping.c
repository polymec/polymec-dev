#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <string.h>
#include "cmockery.h"
#include "geometry/ball_mapping.h"

void plot_points(point_t* points, int num_points, const char* filename)
{
  FILE* fd = fopen(filename, "w");
  fprintf(fd, "# x y z\n");
  for (int i = 0; i < num_points; ++i)
    fprintf(fd, "%g %g %g\n", points[i].x, points[i].y, points[i].z);
  fclose(fd);
}

void write_source(int block_index, point_t* points, int num_points)
{
  printf("  // block %d\n"
         "  {", block_index);
  for (int i = 0; i < num_points; ++i)
  {
    printf("   {.x = %.16f, .y = %.16f, .z = %.16f},\n", points[i].x, points[i].y, points[i].z);
  }
  printf("  },\n");
}

// These test points were generated by write_source() above.
static const double L = 1.0;
static const int num_test_points = 27;
static point_t test_xs[7][27] = {
  // block 0
  {
   {.x = -0.2500000000000000, .y = -0.2500000000000000, .z = -0.2500000000000000},
   {.x = -0.2500000000000000, .y = -0.2500000000000000, .z = 0.0000000000000000},
   {.x = -0.2500000000000000, .y = -0.2500000000000000, .z = 0.2500000000000000},
   {.x = -0.2500000000000000, .y = 0.0000000000000000, .z = -0.2500000000000000},
   {.x = -0.2500000000000000, .y = 0.0000000000000000, .z = 0.0000000000000000},
   {.x = -0.2500000000000000, .y = 0.0000000000000000, .z = 0.2500000000000000},
   {.x = -0.2500000000000000, .y = 0.2500000000000000, .z = -0.2500000000000000},
   {.x = -0.2500000000000000, .y = 0.2500000000000000, .z = 0.0000000000000000},
   {.x = -0.2500000000000000, .y = 0.2500000000000000, .z = 0.2500000000000000},
   {.x = 0.0000000000000000, .y = -0.2500000000000000, .z = -0.2500000000000000},
   {.x = 0.0000000000000000, .y = -0.2500000000000000, .z = 0.0000000000000000},
   {.x = 0.0000000000000000, .y = -0.2500000000000000, .z = 0.2500000000000000},
   {.x = 0.0000000000000000, .y = 0.0000000000000000, .z = -0.2500000000000000},
   {.x = 0.0000000000000000, .y = 0.0000000000000000, .z = 0.0000000000000000},
   {.x = 0.0000000000000000, .y = 0.0000000000000000, .z = 0.2500000000000000},
   {.x = 0.0000000000000000, .y = 0.2500000000000000, .z = -0.2500000000000000},
   {.x = 0.0000000000000000, .y = 0.2500000000000000, .z = 0.0000000000000000},
   {.x = 0.0000000000000000, .y = 0.2500000000000000, .z = 0.2500000000000000},
   {.x = 0.2500000000000000, .y = -0.2500000000000000, .z = -0.2500000000000000},
   {.x = 0.2500000000000000, .y = -0.2500000000000000, .z = 0.0000000000000000},
   {.x = 0.2500000000000000, .y = -0.2500000000000000, .z = 0.2500000000000000},
   {.x = 0.2500000000000000, .y = 0.0000000000000000, .z = -0.2500000000000000},
   {.x = 0.2500000000000000, .y = 0.0000000000000000, .z = 0.0000000000000000},
   {.x = 0.2500000000000000, .y = 0.0000000000000000, .z = 0.2500000000000000},
   {.x = 0.2500000000000000, .y = 0.2500000000000000, .z = -0.2500000000000000},
   {.x = 0.2500000000000000, .y = 0.2500000000000000, .z = 0.0000000000000000},
   {.x = 0.2500000000000000, .y = 0.2500000000000000, .z = 0.2500000000000000}
  },
  // block 1
  {
   {.x = -0.4082482904638630, .y = -0.4082482904638631, .z = -0.8164965809277259},
   {.x = -0.7071067811865475, .y = -0.7071067811865476, .z = 0.0000000000000001},
   {.x = -0.4082482904638630, .y = -0.4082482904638630, .z = 0.8164965809277260},
   {.x = -0.7071067811865476, .y = -0.0000000000000001, .z = -0.7071067811865475},
   {.x = -1.0000000000000000, .y = -0.0000000000000001, .z = 0.0000000000000001},
   {.x = -0.7071067811865476, .y = -0.0000000000000001, .z = 0.7071067811865475},
   {.x = -0.4082482904638632, .y = 0.4082482904638630, .z = -0.8164965809277259},
   {.x = -0.7071067811865477, .y = 0.7071067811865475, .z = 0.0000000000000001},
   {.x = -0.4082482904638632, .y = 0.4082482904638630, .z = 0.8164965809277259},
   {.x = -0.2762929288806347, .y = -0.2762929288806348, .z = -0.5525858577612693},
   {.x = -0.4618065660663286, .y = -0.4618065660663286, .z = 0.0000000000000000},
   {.x = -0.2762929288806347, .y = -0.2762929288806347, .z = 0.5525858577612694},
   {.x = -0.4618065660663286, .y = -0.0000000000000001, .z = -0.4618065660663286},
   {.x = -0.6250000000000000, .y = -0.0000000000000001, .z = 0.0000000000000000},
   {.x = -0.4618065660663286, .y = -0.0000000000000001, .z = 0.4618065660663286},
   {.x = -0.2762929288806348, .y = 0.2762929288806347, .z = -0.5525858577612693},
   {.x = -0.4618065660663287, .y = 0.4618065660663286, .z = 0.0000000000000000},
   {.x = -0.2762929288806348, .y = 0.2762929288806347, .z = 0.5525858577612693},
   {.x = -0.1443375672974064, .y = -0.1443375672974065, .z = -0.2886751345948128},
   {.x = -0.2165063509461096, .y = -0.2165063509461096, .z = 0.0000000000000000},
   {.x = -0.1443375672974064, .y = -0.1443375672974064, .z = 0.2886751345948129},
   {.x = -0.2165063509461096, .y = -0.0000000000000000, .z = -0.2165063509461096},
   {.x = -0.2500000000000000, .y = -0.0000000000000000, .z = 0.0000000000000000},
   {.x = -0.2165063509461096, .y = -0.0000000000000000, .z = 0.2165063509461096},
   {.x = -0.1443375672974065, .y = 0.1443375672974064, .z = -0.2886751345948128},
   {.x = -0.2165063509461097, .y = 0.2165063509461096, .z = 0.0000000000000000},
   {.x = -0.1443375672974065, .y = 0.1443375672974064, .z = 0.2886751345948128}
  },
  // block 2
  {
   {.x = 0.1443375672974065, .y = -0.1443375672974065, .z = -0.2886751345948129},
   {.x = 0.2165063509461096, .y = -0.2165063509461096, .z = 0.0000000000000000},
   {.x = 0.1443375672974065, .y = -0.1443375672974064, .z = 0.2886751345948129},
   {.x = 0.2165063509461096, .y = 0.0000000000000000, .z = -0.2165063509461096},
   {.x = 0.2500000000000000, .y = 0.0000000000000000, .z = 0.0000000000000000},
   {.x = 0.2165063509461096, .y = 0.0000000000000000, .z = 0.2165063509461096},
   {.x = 0.1443375672974065, .y = 0.1443375672974065, .z = -0.2886751345948129},
   {.x = 0.2165063509461096, .y = 0.2165063509461096, .z = 0.0000000000000000},
   {.x = 0.1443375672974065, .y = 0.1443375672974064, .z = 0.2886751345948129},
   {.x = 0.2762929288806348, .y = -0.2762929288806348, .z = -0.5525858577612694},
   {.x = 0.4618065660663286, .y = -0.4618065660663286, .z = 0.0000000000000000},
   {.x = 0.2762929288806348, .y = -0.2762929288806347, .z = 0.5525858577612696},
   {.x = 0.4618065660663286, .y = 0.0000000000000000, .z = -0.4618065660663286},
   {.x = 0.6250000000000000, .y = 0.0000000000000000, .z = 0.0000000000000000},
   {.x = 0.4618065660663286, .y = 0.0000000000000000, .z = 0.4618065660663286},
   {.x = 0.2762929288806348, .y = 0.2762929288806348, .z = -0.5525858577612694},
   {.x = 0.4618065660663286, .y = 0.4618065660663286, .z = 0.0000000000000000},
   {.x = 0.2762929288806348, .y = 0.2762929288806347, .z = 0.5525858577612696},
   {.x = 0.4082482904638631, .y = -0.4082482904638630, .z = -0.8164965809277259},
   {.x = 0.7071067811865476, .y = -0.7071067811865475, .z = 0.0000000000000001},
   {.x = 0.4082482904638630, .y = -0.4082482904638630, .z = 0.8164965809277260},
   {.x = 0.7071067811865476, .y = 0.0000000000000000, .z = -0.7071067811865475},
   {.x = 1.0000000000000000, .y = 0.0000000000000000, .z = 0.0000000000000001},
   {.x = 0.7071067811865476, .y = 0.0000000000000000, .z = 0.7071067811865475},
   {.x = 0.4082482904638631, .y = 0.4082482904638630, .z = -0.8164965809277259},
   {.x = 0.7071067811865476, .y = 0.7071067811865475, .z = 0.0000000000000001},
   {.x = 0.4082482904638630, .y = 0.4082482904638630, .z = 0.8164965809277260}
  },
  // block 3
  {
   {.x = -0.4082482904638630, .y = -0.4082482904638631, .z = -0.8164965809277259},
   {.x = -0.7071067811865475, .y = -0.7071067811865476, .z = 0.0000000000000001},
   {.x = -0.4082482904638630, .y = -0.4082482904638630, .z = 0.8164965809277260},
   {.x = -0.2762929288806347, .y = -0.2762929288806348, .z = -0.5525858577612693},
   {.x = -0.4618065660663286, .y = -0.4618065660663286, .z = 0.0000000000000000},
   {.x = -0.2762929288806347, .y = -0.2762929288806347, .z = 0.5525858577612694},
   {.x = -0.1443375672974064, .y = -0.1443375672974065, .z = -0.2886751345948128},
   {.x = -0.2165063509461096, .y = -0.2165063509461096, .z = 0.0000000000000000},
   {.x = -0.1443375672974064, .y = -0.1443375672974064, .z = 0.2886751345948129},
   {.x = 0.0000000000000000, .y = -0.7071067811865476, .z = -0.7071067811865475},
   {.x = 0.0000000000000001, .y = -1.0000000000000000, .z = 0.0000000000000001},
   {.x = 0.0000000000000000, .y = -0.7071067811865476, .z = 0.7071067811865475},
   {.x = 0.0000000000000000, .y = -0.4618065660663286, .z = -0.4618065660663286},
   {.x = 0.0000000000000000, .y = -0.6250000000000000, .z = 0.0000000000000000},
   {.x = 0.0000000000000000, .y = -0.4618065660663286, .z = 0.4618065660663286},
   {.x = 0.0000000000000000, .y = -0.2165063509461096, .z = -0.2165063509461096},
   {.x = 0.0000000000000000, .y = -0.2500000000000000, .z = 0.0000000000000000},
   {.x = 0.0000000000000000, .y = -0.2165063509461096, .z = 0.2165063509461096},
   {.x = 0.4082482904638631, .y = -0.4082482904638630, .z = -0.8164965809277259},
   {.x = 0.7071067811865476, .y = -0.7071067811865475, .z = 0.0000000000000001},
   {.x = 0.4082482904638630, .y = -0.4082482904638630, .z = 0.8164965809277260},
   {.x = 0.2762929288806348, .y = -0.2762929288806347, .z = -0.5525858577612693},
   {.x = 0.4618065660663286, .y = -0.4618065660663286, .z = 0.0000000000000000},
   {.x = 0.2762929288806347, .y = -0.2762929288806347, .z = 0.5525858577612694},
   {.x = 0.1443375672974065, .y = -0.1443375672974064, .z = -0.2886751345948128},
   {.x = 0.2165063509461096, .y = -0.2165063509461096, .z = 0.0000000000000000},
   {.x = 0.1443375672974064, .y = -0.1443375672974064, .z = 0.2886751345948129}
  },
  // block 4
  {
   {.x = -0.1443375672974065, .y = 0.1443375672974065, .z = -0.2886751345948129},
   {.x = -0.2165063509461096, .y = 0.2165063509461096, .z = 0.0000000000000000},
   {.x = -0.1443375672974064, .y = 0.1443375672974065, .z = 0.2886751345948129},
   {.x = -0.2762929288806348, .y = 0.2762929288806348, .z = -0.5525858577612694},
   {.x = -0.4618065660663286, .y = 0.4618065660663286, .z = 0.0000000000000000},
   {.x = -0.2762929288806347, .y = 0.2762929288806348, .z = 0.5525858577612696},
   {.x = -0.4082482904638630, .y = 0.4082482904638631, .z = -0.8164965809277259},
   {.x = -0.7071067811865475, .y = 0.7071067811865476, .z = 0.0000000000000001},
   {.x = -0.4082482904638630, .y = 0.4082482904638630, .z = 0.8164965809277260},
   {.x = 0.0000000000000000, .y = 0.2165063509461096, .z = -0.2165063509461096},
   {.x = 0.0000000000000000, .y = 0.2500000000000000, .z = 0.0000000000000000},
   {.x = 0.0000000000000000, .y = 0.2165063509461096, .z = 0.2165063509461096},
   {.x = 0.0000000000000000, .y = 0.4618065660663286, .z = -0.4618065660663286},
   {.x = 0.0000000000000000, .y = 0.6250000000000000, .z = 0.0000000000000000},
   {.x = 0.0000000000000000, .y = 0.4618065660663286, .z = 0.4618065660663286},
   {.x = 0.0000000000000000, .y = 0.7071067811865476, .z = -0.7071067811865475},
   {.x = 0.0000000000000001, .y = 1.0000000000000000, .z = 0.0000000000000001},
   {.x = 0.0000000000000000, .y = 0.7071067811865476, .z = 0.7071067811865475},
   {.x = 0.1443375672974065, .y = 0.1443375672974065, .z = -0.2886751345948129},
   {.x = 0.2165063509461096, .y = 0.2165063509461096, .z = 0.0000000000000000},
   {.x = 0.1443375672974065, .y = 0.1443375672974064, .z = 0.2886751345948129},
   {.x = 0.2762929288806348, .y = 0.2762929288806348, .z = -0.5525858577612694},
   {.x = 0.4618065660663286, .y = 0.4618065660663286, .z = 0.0000000000000000},
   {.x = 0.2762929288806348, .y = 0.2762929288806347, .z = 0.5525858577612696},
   {.x = 0.4082482904638631, .y = 0.4082482904638630, .z = -0.8164965809277259},
   {.x = 0.7071067811865476, .y = 0.7071067811865475, .z = 0.0000000000000001},
   {.x = 0.4082482904638630, .y = 0.4082482904638630, .z = 0.8164965809277260}
  },
  // block 5
  {
   {.x = -0.2041241452319315, .y = -0.2041241452319316, .z = 0.2041241452319315},
   {.x = -0.3907372072107786, .y = -0.3907372072107787, .z = 0.3907372072107786},
   {.x = -0.5773502691896257, .y = -0.5773502691896258, .z = 0.5773502691896257},
   {.x = -0.2165063509461096, .y = 0.0000000000000000, .z = 0.2165063509461096},
   {.x = -0.4618065660663286, .y = 0.0000000000000001, .z = 0.4618065660663286},
   {.x = -0.7071067811865475, .y = 0.0000000000000001, .z = 0.7071067811865476},
   {.x = -0.2041241452319315, .y = 0.2041241452319316, .z = 0.2041241452319315},
   {.x = -0.3907372072107786, .y = 0.3907372072107787, .z = 0.3907372072107786},
   {.x = -0.5773502691896257, .y = 0.5773502691896258, .z = 0.5773502691896257},
   {.x = 0.0000000000000000, .y = 0.0000000000000000, .z = 0.3061862178478972},
   {.x = 0.0000000000000000, .y = 0.0000000000000000, .z = 0.6530931089239487},
   {.x = 0.0000000000000000, .y = 0.0000000000000000, .z = 1.0000000000000000},
   {.x = 0.0000000000000000, .y = 0.0000000000000000, .z = 0.2500000000000000},
   {.x = 0.0000000000000000, .y = 0.0000000000000000, .z = 0.6250000000000000},
   {.x = 0.0000000000000000, .y = 0.0000000000000000, .z = 1.0000000000000000},
   {.x = 0.0000000000000000, .y = 0.0000000000000000, .z = 0.3061862178478972},
   {.x = 0.0000000000000000, .y = 0.0000000000000000, .z = 0.6530931089239487},
   {.x = 0.0000000000000000, .y = 0.0000000000000000, .z = 1.0000000000000000},
   {.x = 0.2041241452319315, .y = -0.2041241452319315, .z = 0.2041241452319315},
   {.x = 0.3907372072107786, .y = -0.3907372072107786, .z = 0.3907372072107788},
   {.x = 0.5773502691896257, .y = -0.5773502691896256, .z = 0.5773502691896258},
   {.x = 0.2165063509461096, .y = 0.0000000000000000, .z = 0.2165063509461096},
   {.x = 0.4618065660663286, .y = 0.0000000000000000, .z = 0.4618065660663286},
   {.x = 0.7071067811865475, .y = 0.0000000000000000, .z = 0.7071067811865476},
   {.x = 0.2041241452319315, .y = 0.2041241452319315, .z = 0.2041241452319315},
   {.x = 0.3907372072107786, .y = 0.3907372072107786, .z = 0.3907372072107788},
   {.x = 0.5773502691896257, .y = 0.5773502691896256, .z = 0.5773502691896258}
  },
  // block 6
  {
   {.x = -0.5773502691896257, .y = -0.5773502691896258, .z = -0.5773502691896257},
   {.x = -0.3907372072107786, .y = -0.3907372072107786, .z = -0.3907372072107786},
   {.x = -0.2041241452319315, .y = -0.2041241452319315, .z = -0.2041241452319315},
   {.x = -0.7071067811865476, .y = 0.0000000000000001, .z = -0.7071067811865475},
   {.x = -0.4618065660663286, .y = 0.0000000000000001, .z = -0.4618065660663286},
   {.x = -0.2165063509461096, .y = 0.0000000000000000, .z = -0.2165063509461096},
   {.x = -0.5773502691896257, .y = 0.5773502691896258, .z = -0.5773502691896257},
   {.x = -0.3907372072107786, .y = 0.3907372072107786, .z = -0.3907372072107786},
   {.x = -0.2041241452319315, .y = 0.2041241452319315, .z = -0.2041241452319315},
   {.x = 0.0000000000000000, .y = -0.0000000000000001, .z = -1.0000000000000000},
   {.x = 0.0000000000000000, .y = -0.0000000000000001, .z = -0.6530931089239487},
   {.x = 0.0000000000000000, .y = -0.0000000000000000, .z = -0.3061862178478972},
   {.x = 0.0000000000000001, .y = 0.0000000000000000, .z = -1.0000000000000000},
   {.x = 0.0000000000000001, .y = 0.0000000000000000, .z = -0.6250000000000000},
   {.x = 0.0000000000000000, .y = 0.0000000000000000, .z = -0.2500000000000000},
   {.x = 0.0000000000000000, .y = 0.0000000000000001, .z = -1.0000000000000000},
   {.x = 0.0000000000000000, .y = 0.0000000000000001, .z = -0.6530931089239487},
   {.x = 0.0000000000000000, .y = 0.0000000000000000, .z = -0.3061862178478972},
   {.x = 0.5773502691896258, .y = -0.5773502691896257, .z = -0.5773502691896257},
   {.x = 0.3907372072107786, .y = -0.3907372072107786, .z = -0.3907372072107786},
   {.x = 0.2041241452319315, .y = -0.2041241452319315, .z = -0.2041241452319315},
   {.x = 0.7071067811865476, .y = 0.0000000000000000, .z = -0.7071067811865475},
   {.x = 0.4618065660663286, .y = 0.0000000000000000, .z = -0.4618065660663286},
   {.x = 0.2165063509461096, .y = 0.0000000000000000, .z = -0.2165063509461096},
   {.x = 0.5773502691896258, .y = 0.5773502691896257, .z = -0.5773502691896257},
   {.x = 0.3907372072107786, .y = 0.3907372072107786, .z = -0.3907372072107786},
   {.x = 0.2041241452319315, .y = 0.2041241452319315, .z = -0.2041241452319315}
  }
};

static void test_block(int block_index)
{
  point_t zero = {.x = 0.0, .y = 0.0, .z = 0.0};
  sp_func_t* mapping = ball_mapping_new(&zero, 1.0, 0.5, block_index);
  point_t test_points[num_test_points], xis[num_test_points];

  // Construct the test points in logical space.
  int offset = 0;
  for (int i = 0; i < 3; ++i)
  {
    for (int j = 0; j < 3; ++j)
    {
      for (int k = 0; k < 3; ++k, ++offset)
      {
        xis[offset].x = 0.5*i;
        xis[offset].y = 0.5*j;
        xis[offset].z = 0.5*k;
      }
    }
  }
  ASSERT(offset == num_test_points);

  for (int i = 0; i < num_test_points; ++i)
  {
    double X[3];
    sp_func_eval(mapping, &xis[i], X);
    test_points[i].x = X[0];
    test_points[i].y = X[1];
    test_points[i].z = X[2];

    double error = point_distance(&test_points[i], &test_xs[block_index][i]);
    assert_true(error < 1e-15);
  }

  //char fn[1024];
  //snprintf(fn, 1024, "points.%d", block_index);
  //plot_points(test_points, 27, fn);
  //write_source(block_index, test_points, 27);
}

void test_block0(void** state)
{
  test_block(0);
}

void test_block1(void** state)
{
  test_block(1);
}

void test_block2(void** state)
{
  test_block(2);
}

void test_block3(void** state)
{
  test_block(3);
}

void test_block4(void** state)
{
  test_block(4);
}

void test_block5(void** state)
{
  test_block(5);
}

void test_block6(void** state)
{
  test_block(6);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const UnitTest tests[] = 
  {
    unit_test(test_block0),
    unit_test(test_block1),
    unit_test(test_block2),
    unit_test(test_block3),
    unit_test(test_block4),
    unit_test(test_block5),
    unit_test(test_block6)
  };
  return run_tests(tests);
}
