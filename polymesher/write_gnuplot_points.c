// Copyright (c) 2012-2013, Jeffrey N. Johnson
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without 
// modification, are permitted provided that the following conditions are met:
// 
// 1. Redistributions of source code must retain the above copyright notice, this 
// list of conditions and the following disclaimer.
// 
// 2. Redistributions in binary form must reproduce the above copyright notice, 
// this list of conditions and the following disclaimer in the documentation 
// and/or other materials provided with the distribution.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE 
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

// This implements polymesher's capability for writing 3D scatter point data  
// readable by Gnuplot.

#include <string.h>
#include "core/polymec.h"
#include "core/interpreter.h"

// Lua stuff.
#include "lua.h"
#include "lualib.h"
#include "lauxlib.h"

// write_gnuplot_points(args) -- This function writes a given list of points
// to a file on disk.
//
int write_gnuplot_points(lua_State* lua)
{
  // Rank 0 only.
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank != 0) return 0;

  // Check the arguments.
  int num_args = lua_gettop(lua);
  if ((num_args != 2) || 
      ((num_args == 2) && (!lua_ispointlist(lua, 1) || !lua_isstring(lua, 2))))
  {
    return luaL_error(lua, "write_gnuplot_points: invalid arguments. Usage:\n"
                      "write_gnuplot_points(points, filename)");
  }

  // Get the argument(s).
  int num_points;
  point_t* points = lua_topointlist(lua, 1, &num_points);
  const char* prefix = lua_tostring(lua, 2);

  log_info("Writing GNUPlot scatter plot with prefix '%s'...", prefix);

  // Write the data.
  char filename[1024];
  snprintf(filename, 1024, "%s.gnuplot", prefix);
  FILE* fd = fopen(filename, "w");
  fprintf(fd, "# x y z\n");
  for (int i = 0; i < num_points; ++i)
    fprintf(fd, "%g %g %g\n", points[i].x, points[i].y, points[i].z);
  fclose(fd);

  return 0;
}

