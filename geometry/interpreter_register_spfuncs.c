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

#include "core/interpreter.h"
#include "geometry/sphere.h"
#include "geometry/rect_prism.h"

// Lua stuff.
#include "lua.h"
#include "lualib.h"
#include "lauxlib.h"

static int sphere(lua_State* lua)
{
  // Check the arguments.
  int num_args = lua_gettop(lua);
  if ((num_args != 2) || !lua_ispoint(lua, 1) || !lua_isnumber(lua, 2))
  {
    return luaL_error(lua, "Invalid argument(s). Usage:\n"
                      "F = sphere(x, r)");
  }

  // Get the arguments.
  point_t* x = lua_topoint(lua, 1);
  double r = (double)lua_tonumber(lua, 2);
  if (r <= 0.0)
    return luaL_error(lua, "Sphere radius must be positive.");

  sp_func_t* s = sphere_new(x, r, INWARD_NORMAL);
  lua_pushscalarfunction(lua, st_func_from_sp_func(s));
  return 1;
}

static int rect_prism(lua_State* lua)
{
  // Check the arguments.
  int num_args = lua_gettop(lua);
  if ((num_args != 1) || !lua_isboundingbox(lua, 1))
  {
    return luaL_error(lua, "Invalid argument(s). Usage:\n"
                      "F = rect_prism(bbox)");
  }

  // Get the arguments.
  bbox_t* bbox = lua_toboundingbox(lua, 1);
  sp_func_t* prism = rect_prism_new_from_bbox(bbox);
  lua_pushscalarfunction(lua, st_func_from_sp_func(prism));
  prism = NULL;
  return 1;
}

void interpreter_register_spfuncs(interpreter_t* interp)
{
  interpreter_register_function(interp, "sphere", sphere);
  interpreter_register_function(interp, "rect_prism", rect_prism);
}

