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

