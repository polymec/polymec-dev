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

// This implements a command to view the nodes in a generated mesh.

#include <string.h>
#include "core/polymec.h"
#include "core/interpreter.h"
#include "core/point.h"

// Lua stuff.
#include "lua.h"
#include "lualib.h"
#include "lauxlib.h"

// mesh_nodes(mesh) -- This function returns a pointlist representing the 
// nodes of the given mesh.
int mesh_nodes(lua_State* lua)
{
  // Check the arguments.
  int num_args = lua_gettop(lua);
  if ((num_args != 1) || !lua_ismesh(lua, 1))
  {
    return luaL_error(lua, "mesh_nodes: invalid arguments. Usage:\n"
                      "nodes = mesh_nodes(mesh)");
  }

  // Get the argument(s).
  mesh_t* mesh = lua_tomesh(lua, 1);

  point_t* nodes = malloc(sizeof(point_t) * mesh->num_nodes);
  for (int n = 0; n < mesh->num_nodes; ++n)
  {
    nodes[n].x = mesh->nodes[n].x;
    nodes[n].y = mesh->nodes[n].y;
    nodes[n].z = mesh->nodes[n].z;
  }

  lua_pushpointlist(lua, nodes, mesh->num_nodes);
  return 1;
}


