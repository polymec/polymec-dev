#include "core/interpreter.h"
#include "geometry/sphere.h"

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
    lua_pushstring(lua, "Invalid argument(s). Usage:\n"
                        "F = sphere(x, r)");
    lua_error(lua);
    return LUA_ERRRUN;
  }

  // Get the arguments.
  point_t* x = lua_topoint(lua, 1);
  double r = (double)lua_tonumber(lua, 2);
  if (r <= 0.0)
  {
    lua_pushstring(lua, "Sphere radius must be positive.");
    lua_error(lua);
  }

  sp_func_t* s = sphere_new(x, r, INWARD_NORMAL);
  lua_pushscalarfunction(lua, st_func_from_spfunc(s));
  return 1;
}

void interpreter_register_spfuncs(interpreter_t* interp)
{
  interpreter_register_function(interp, "sphere", sphere);
}

