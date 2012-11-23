#include "core/interpreter.h"
#include "core/unordered_map.h"

#ifdef __cplusplus
extern "C" {
#endif

// Lua stuff.
#include "lua.h"
#include "lualib.h"
#include "lauxlib.h"

// Interpreter data structure.
struct interpreter_t
{
  // The Lua interpreter.
  lua_State* lua;

  // The data store.
  str_ptr_unordered_map_t* store;

  // A list of valid inputs and their types.
  int num_valid_inputs;
  interpreter_validation_t* valid_inputs;
};

interpreter_t* interpreter_new(interpreter_validation_t* valid_inputs, int num_valid_inputs)
{
  ASSERT(num_valid_inputs >= 0);
  ASSERT((valid_inputs != NULL) || (num_valid_inputs == 0));

  interpreter_t* interp = malloc(sizeof(interpreter_t));

  // Initialize the Lua interpreter.
  interp->lua = luaL_newstate();
  luaL_openlibs(interp->lua);

  // Add some functions for creating data types.
  // FIXME
//  lua_register(interp->lua, "cubic_mesh", lua_cubic_mesh);

  // Initialize the data store.
  interp->store = str_ptr_unordered_map_new();

  // Copy over the valid inputs.
  interp->valid_inputs = malloc(num_valid_inputs*sizeof(interpreter_validation_t));
  interp->num_valid_inputs = num_valid_inputs;
  for (int i = 0; i < num_valid_inputs; ++i)
  {
    interp->valid_inputs[i].variable = strdup(valid_inputs[i].variable);
    interp->valid_inputs[i].type = valid_inputs[i].type;
  }

  return interp;
}

void interpreter_free(interpreter_t* interp)
{
  for (int i = 0; i < interp->num_valid_inputs; ++i)
    free(interp->valid_inputs[i].variable);
  free(interp->valid_inputs);
  str_ptr_unordered_map_free(interp->store);
  ASSERT(interp->lua != NULL);
  lua_close(interp->lua);
  free(interp);
}

void interpreter_parse(interpreter_t* interp, FILE* input)
{
}

void* interpreter_get(interpreter_t* interp, const char* name)
{
  return NULL;
}

#ifdef __cplusplus
}
#endif

