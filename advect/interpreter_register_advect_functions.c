#include "advect/interpreter_register_advect_functions.h"
#include "core/constant_st_func.h"

#ifdef __cplusplus
extern "C" {
#endif

// Lua stuff.
#include "lua.h"
#include "lualib.h"
#include "lauxlib.h"

static int dirichlet_bc(lua_State* lua)
{
  // Check the argument.
  int num_args = lua_gettop(lua);
  if ((num_args != 1) || (!lua_isstfunc(lua, 1) && !lua_isnumber(lua, 1)))
  {
    lua_pushstring(lua, "Invalid arguments. Usage:\nbc = dirichlet_bc(F)\nwhere F is a number or function.");
    lua_error(lua);
    return LUA_ERRRUN;
  }

  // Get the argument. 
  st_func_t* F;
  if (lua_isnumber(lua, 1))
  {
    double F0 = lua_tonumber(lua, 1);
    F = constant_st_func_new(1, &F0);
  }
  else
  {
    F = lua_tostfunc(lua, 1);
  }

  // Create a boundary condition object and push it onto the stack.
  advect_bc_t* bc = advect_bc_new(1.0, 0.0, F);
  lua_pushuserdefined(lua, bc, advect_bc_free);
  return 1;
}

static int neumann_bc(lua_State* lua)
{
  // Check the argument.
  int num_args = lua_gettop(lua);
  if ((num_args != 1) || (!lua_isstfunc(lua, 1) && !lua_isnumber(lua, 1)))
  {
    lua_pushstring(lua, "Invalid arguments. Usage:\nbc = neumann_bc(F)\nwhere F is a number or function.");
    lua_error(lua);
    return LUA_ERRRUN;
  }

  // Get the argument. 
  st_func_t* F;
  if (lua_isnumber(lua, 1))
  {
    double F0 = lua_tonumber(lua, 1);
    F = constant_st_func_new(1, &F0);
  }
  else
  {
    F = lua_tostfunc(lua, 1);
  }

  // Create a boundary condition object and push it onto the stack.
  advect_bc_t* bc = advect_bc_new(0.0, 1.0, F);
  lua_pushuserdefined(lua, bc, advect_bc_free);
  return 1;
}

static int robin_bc(lua_State* lua)
{
  // Check the argument.
  int num_args = lua_gettop(lua);
  if ((num_args != 3) || 
      !lua_isnumber(lua, 1) ||
      !lua_isnumber(lua, 2) ||
      (!lua_isstfunc(lua, 3) && !lua_isnumber(lua, 3)))
  {
    lua_pushstring(lua, "Invalid arguments. Usage:\nbc = robin_bc(alpha, beta, F)\nwhere F is a number or function.");
    lua_error(lua);
    return LUA_ERRRUN;
  }

  // Get the arguments. 
  double alpha = lua_tonumber(lua, 1);
  double beta = lua_tonumber(lua, 2);
  st_func_t* F;
  if (lua_isnumber(lua, 3))
  {
    double F0 = lua_tonumber(lua, 3);
    F = constant_st_func_new(1, &F0);
  }
  else
  {
    F = lua_tostfunc(lua, 3);
  }

  // Create a boundary condition object and push it onto the stack.
  advect_bc_t* bc = advect_bc_new(alpha, beta, F);
  lua_pushuserdefined(lua, bc, advect_bc_free);
  return 1;
}

void interpreter_register_advect_functions(interpreter_t* interpreter)
{
  interpreter_register_function(interpreter, "dirichlet_bc", dirichlet_bc);
  interpreter_register_function(interpreter, "neumann_bc", neumann_bc);
  interpreter_register_function(interpreter, "robin_bc", robin_bc);
}

#ifdef __cplusplus
}
#endif

