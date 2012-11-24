#include "core/interpreter.h"
#include "core/unordered_map.h"

#ifdef __cplusplus
extern "C" {
#endif

// Lua stuff.
#include "lua.h"
#include "lualib.h"
#include "lauxlib.h"

// Define a map from variable names to storage.
typedef struct 
{
  void *datum;                 // The datum associated with a variable.
  interpreter_var_type_t type; // The data type.
  void (*dtor)(void*);         // Data destructor.
} interpreter_storage_t;

DEFINE_UNORDERED_MAP(interpreter_map, char*, interpreter_storage_t*, string_hash, string_equals)

static void destroy_variable(char* key, interpreter_storage_t* value)
{
  free(key);
  if (value->dtor)
    (*value->dtor)(value->datum);
  free(value);
}

static interpreter_storage_t* store_string(const char* var)
{
  interpreter_storage_t* storage = malloc(sizeof(interpreter_storage_t));
  storage->datum = strdup(var);
  storage->type = INTERPRETER_STRING;
  storage->dtor = free;
  return storage;
}

static interpreter_storage_t* store_number(double var)
{
  interpreter_storage_t* storage = malloc(sizeof(interpreter_storage_t));
  double* dvar = malloc(sizeof(double));
  *dvar = var;
  storage->datum = dvar;
  storage->type = INTERPRETER_NUMBER;
  storage->dtor = free;
  return storage;
}

static void destroy_table(void* table)
{
  str_ptr_unordered_map_free((str_ptr_unordered_map_t*)table);
}

static interpreter_storage_t* store_table(str_ptr_unordered_map_t* table)
{
  interpreter_storage_t* storage = malloc(sizeof(interpreter_storage_t));
  storage->datum = table;
  storage->type = INTERPRETER_TABLE;
  storage->dtor = destroy_table;
  return storage;
}

// Interpreter data structure.
struct interpreter_t
{
  // The Lua interpreter.
  lua_State* lua;

  // The data store.
  interpreter_map_t* store;

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
  interp->store = interpreter_map_new();

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
  interpreter_map_free(interp->store);
  ASSERT(interp->lua != NULL);
  lua_close(interp->lua);
  free(interp);
}

static interpreter_validation_t* interpreter_validation_entry(interpreter_t* interp, const char* key)
{
  for (int i = 0; i < interp->num_valid_inputs; ++i)
  {
    if (!strcmp(key, interp->valid_inputs[i].variable))
      return &interp->valid_inputs[i];
  }
  return NULL;
}

void interpreter_parse(interpreter_t* interp, char* input)
{
  int error = luaL_dostring(interp->lua, input);
  if (error == LUA_ERRSYNTAX)
    arbi_error("Syntax error in input.");
  else if (error != LUA_OK)
    arbi_error("Internal error in interpreter.");

  // Traverse Lua's table of global variables and pull them into the
  // data store.
  lua_pushnil(interp->lua); // nil tells lua_next to start at the first key.
  while (lua_next(interp->lua, LUA_RIDX_GLOBALS)) // Traverse globals table.
  {
    // key is at index -2, value is at index -1.
    const char* key = lua_tostring(interp->lua, -2);

    // Does the key appear in our validation table?
    interpreter_validation_t* entry = interpreter_validation_entry(interp, key);
    if (entry != NULL)
    {
      // We must validate this variable against its allowed type.
      if ((entry->type == INTERPRETER_STRING) && !lua_isstring(interp->lua, -2))
        arbi_error("Type error: %s must be a string.", key);
      else if ((entry->type == INTERPRETER_NUMBER) && !lua_isnumber(interp->lua, -2))
        arbi_error("Type error: %s must be a number.", key);
      else if (entry->type == INTERPRETER_MESH)
      {
        if (!lua_islightuserdata(interp->lua, -2))
          arbi_error("Type error: %s must be a mesh.", key);
        interpreter_storage_t* var = (void*)lua_topointer(interp->lua, -2);
        if (var->type != INTERPRETER_MESH)
          arbi_error("Type error: %s must be a mesh.", key);
      }
      else if (entry->type == INTERPRETER_SP_FUNC)
      {
        if (!lua_islightuserdata(interp->lua, -2))
          arbi_error("Type error: %s must be a spatial function.", key);
        interpreter_storage_t* var = (void*)lua_topointer(interp->lua, -2);
        if (var->type != INTERPRETER_SP_FUNC)
          arbi_error("Type error: %s must be a spatial function.", key);
      }
      else if (entry->type == INTERPRETER_ST_FUNC)
      {
        if (!lua_islightuserdata(interp->lua, -2))
          arbi_error("Type error: %s must be a space-time function.", key);
        interpreter_storage_t* var = (void*)lua_topointer(interp->lua, -2);
        if (var->type != INTERPRETER_ST_FUNC)
          arbi_error("Type error: %s must be a space-time function.", key);
      }
      else if ((entry->type == INTERPRETER_TABLE) && !lua_istable(interp->lua, -2))
        arbi_error("Type error: %s must be a table mapping strings to objects.", key);
    }

    interpreter_storage_t* var = NULL;
    if (lua_isstring(interp->lua, -2))
      var = store_string(lua_tostring(interp->lua, -2));
    else if (lua_isnumber(interp->lua, -2))
      var = store_number(lua_tonumber(interp->lua, -2));
    else if (lua_istable(interp->lua, -2))
    {
      str_ptr_unordered_map_t* table = str_ptr_unordered_map_new();
      // Traverse this table.
      lua_pushnil(interp->lua);
      while (lua_next(interp->lua, -3)) // FIXME??
      {
        // Key is at index -4, value is at -5.
        if (!lua_isstring(interp->lua, -4))
          arbi_error("Type error: %s must be a table mapping strings to objects.", key);
        if (!lua_islightuserdata(interp->lua, -5))
          arbi_error("Type error: %s must be a table mapping strings to objects.", key);
        char* tkey = (char*)lua_tostring(interp->lua, -4);
        void* tval = (void*)lua_topointer(interp->lua, -5);
        // FIXME: Objects in the table must currently be garbage-collected.
        str_ptr_unordered_map_insert(table, tkey, tval);

        // Removes value from stack.
        lua_pop(interp->lua, 1);
      }
      var = store_table(table);
    }
    else 
    {
      ASSERT(lua_islightuserdata(interp->lua, -2));
      var = (void*)lua_topointer(interp->lua, -2);
    }
    interpreter_map_insert_with_dtor(interp->store, (char*)key, var, destroy_variable);

    // Removes value from stack -- key is kept for next iteration.
    lua_pop(interp->lua, 1);
  }
  lua_gettable(interp->lua, LUA_RIDX_GLOBALS);
}

void* interpreter_get(interpreter_t* interp, const char* name)
{
  return NULL;
}

#ifdef __cplusplus
}
#endif

