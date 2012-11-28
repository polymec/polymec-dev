#include "core/interpreter.h"
#include "core/constant_st_func.h"
#include "core/unordered_set.h"

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
  value->datum = NULL;
  free(value);
}

static void destroy_table_entry(char* key, void* value)
{
  free(key);
  free(value);
}

static void destroy_mesh(void* mesh)
{
  mesh_free((mesh_t*)mesh);
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
  // A registry of functions to extend the interpreter.
  int num_functions;
  char* function_names[1024];
  lua_CFunction functions[1024];

  // The data store.
  interpreter_map_t* store;

  // A list of valid inputs and their types.
  int num_valid_inputs;
  interpreter_validation_t* valid_inputs;

  // A set of pre-existing variables in Lua that will be ignored.
  str_unordered_set_t* preexisting_vars;
};

// Creates a constant (scalar-valued) function from a number.
static int constant_function(lua_State* lua)
{
  // Check the argument.
  int num_args = lua_gettop(lua);
  if ((num_args != 1) || !lua_isnumber(lua, 1))
  {
    lua_pushstring(lua, "Argument must be a number.");
    lua_error(lua);
    return LUA_ERRRUN;
  }

  // Get the argument.
  double arg = lua_tonumber(lua, 1);

  // Push a constant function onto the stack.
  st_func_t* func = constant_st_func_new(1, &arg);
  lua_pushstfunc(lua, func);
  return 1;
}

static void register_default_functions(interpreter_t* interp)
{
  interpreter_register_function(interp, "constant_function", constant_function);
}

interpreter_t* interpreter_new(interpreter_validation_t* valid_inputs)
{
  interpreter_t* interp = malloc(sizeof(interpreter_t));
  interp->num_functions = 0;

  // Add the default functions to the interpreter.
  register_default_functions(interp);

  // Initialize the data store.
  interp->store = interpreter_map_new();

  // Copy over the valid inputs.
  if (valid_inputs != NULL)
  {
    int num_valid_inputs = 0;
    while (valid_inputs[num_valid_inputs].type != INTERPRETER_TERMINUS)
      num_valid_inputs++;
    interp->valid_inputs = malloc(num_valid_inputs*sizeof(interpreter_validation_t));
    interp->num_valid_inputs = num_valid_inputs;
  }
  else
  {
    interp->valid_inputs = NULL;
    interp->num_valid_inputs = 0;
  }
  for (int i = 0; i < interp->num_valid_inputs; ++i)
  {
    interp->valid_inputs[i].variable = strdup(valid_inputs[i].variable);
    interp->valid_inputs[i].type = valid_inputs[i].type;
  }

  interp->preexisting_vars = NULL;

  return interp;
}

void interpreter_free(interpreter_t* interp)
{
  for (int i = 0; i < interp->num_functions; ++i)
    free(interp->function_names[i]);
  for (int i = 0; i < interp->num_valid_inputs; ++i)
    free(interp->valid_inputs[i].variable);
  free(interp->valid_inputs);
  interpreter_map_free(interp->store);
  free(interp);
}

void interpreter_register_function(interpreter_t* interp, const char* function_name, int (*function)(lua_State*))
{
  ASSERT(interp->num_functions < 1024);
  interp->function_names[interp->num_functions] = strdup(function_name);
  interp->functions[interp->num_functions] = function;
  interp->num_functions++;
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

static lua_State* interpreter_open_lua(interpreter_t* interp)
{
  // Initialize the Lua interpreter.
  lua_State* lua = luaL_newstate();
  ASSERT(lua != NULL);
  luaL_openlibs(lua);

  // Add whatever other functions we've had registered.
  for (int i = 0; i < interp->num_functions; ++i)
    lua_register(lua, (const char*)interp->function_names[i], interp->functions[i]);

  // Before we go adding our own variables to the interpreter, make a 
  // list of what's there so we can ignore it.
  // FIXME: This is a hack, but it's easier than learning about 
  // FIXME: lua's notions of environments, etc.
  ASSERT(interp->preexisting_vars == NULL);
  interp->preexisting_vars = str_unordered_set_new();
  lua_pushglobaltable(lua);
  ASSERT(lua_istable(lua, -1));
  lua_pushnil(lua); 
  while (lua_next(lua, -2)) // Traverse globals table.
  {
    const char* key = lua_tostring(lua, -2);
    str_unordered_set_insert(interp->preexisting_vars, (char*)key);
    lua_pop(lua, 1);
  }
  lua_pop(lua, 1);

  return lua;
}

static void interpreter_close_lua(interpreter_t* interp, lua_State* lua)
{
  str_unordered_set_free(interp->preexisting_vars);
  interp->preexisting_vars = NULL;
  lua_close(lua);
}

static void interpreter_store_chunk_contents(interpreter_t* interp, lua_State* lua)
{
  // Get the global variables table from the registry.
  lua_pushglobaltable(lua);
//  lua_rawgeti(lua, LUA_REGISTRYINDEX, LUA_RIDX_GLOBALS);
  ASSERT(lua_istable(lua, -1));

  // Traverse the table.
  lua_pushnil(lua); 
  while (lua_next(lua, -2)) // Traverse globals table.
  {
    // key is at index -2, value is at index -1.
    static const int key_index = -2;
    static const int val_index = -1;
    const char* key = lua_tostring(lua, key_index);

    // Check to see if this variable existed in Lua before we 
    // parsed our inputs.
    bool preexisting_var = str_unordered_set_contains(interp->preexisting_vars, (char*)key);
    bool skip_this_var = false;

    // Does the key appear in our validation table?
    interpreter_validation_t* entry = interpreter_validation_entry(interp, key);
    if (entry != NULL)
    {
      // We must validate this variable against its allowed type.
      if ((entry->type == INTERPRETER_STRING) && !lua_isstring(lua, val_index))
      {
        if (preexisting_var)
          skip_this_var = true;
        else
          polymec_error("Type error: %s must be a string.", key);
      }
      else if ((entry->type == INTERPRETER_NUMBER) && !lua_isnumber(lua, val_index))
      {
        if (preexisting_var)
          skip_this_var = true;
        else
          polymec_error("Type error: %s must be a number.", key);
      }
      else if (entry->type == INTERPRETER_MESH)
      {
        if (!lua_islightuserdata(lua, val_index))
        {
          if (preexisting_var)
            skip_this_var = true;
          else
            polymec_error("Type error: %s must be a mesh.", key);
        }
        interpreter_storage_t* var = (void*)lua_topointer(lua, val_index);
        if (var->type != INTERPRETER_MESH)
        {
          if (preexisting_var)
            skip_this_var = true;
          else
            polymec_error("Type error: %s must be a mesh.", key);
        }
      }
      else if (entry->type == INTERPRETER_FUNCTION)
      {
        if (!lua_islightuserdata(lua, val_index))
        {
          if (preexisting_var)
            skip_this_var = true;
          else
            polymec_error("Type error: %s must be a function.", key);
        }
        interpreter_storage_t* var = (void*)lua_topointer(lua, val_index);
        if (var->type != INTERPRETER_FUNCTION)
        {
          if (preexisting_var)
            skip_this_var = true;
          else
            polymec_error("Type error: %s must be a function.", key);
        }
      }
      else if ((entry->type == INTERPRETER_TABLE) && !lua_istable(lua, val_index))
      {
        if (preexisting_var)
          skip_this_var = true;
        else
          polymec_error("Type error: %s must be a table mapping strings to objects.", key);
      }
    }

    // Skip this variable if we need to.
    if (skip_this_var)
    {
      lua_pop(lua, 1);
      continue;
    }

    interpreter_storage_t* var = NULL;
    if (lua_isnumber(lua, val_index))
      var = store_number(lua_tonumber(lua, val_index));
    else if (lua_isstring(lua, val_index))
      var = store_string(lua_tostring(lua, val_index));
    else if (lua_istable(lua, val_index))
    {
      // Before we do anything, we validate the table.
      interpreter_var_type_t value_type = INTERPRETER_TERMINUS;
      static const char* type_names[] = {"string", "number", "mesh", "function"};
      // Traverse this table and make sure its values are all of one type.
      lua_pushnil(lua);
      while (lua_next(lua, -2))
      {
        // Key is at index -2, value is at -1.
        static const int key_index = -2;
        static const int val_index = -1;
        if (!lua_isstring(lua, key_index))
        {
          if (preexisting_var)
            skip_this_var = true;
          else
            polymec_error("Type error: %s must be a table mapping strings to objects.", key);
        }
        if (!lua_isnumber(lua, val_index) && 
            !lua_isstring(lua, val_index) && 
            !lua_islightuserdata(lua, val_index))
        {
          if (preexisting_var)
            skip_this_var = true;
          else
            polymec_error("Type error: %s must be a table mapping strings to objects.", key);
        }
        if (skip_this_var)
        {
          lua_pop(lua, 1);
          continue;
        }

        char* tkey = (char*)lua_tostring(lua, key_index);
        void* tval = (void*)lua_topointer(lua, val_index);
        if (tval == NULL)
        {
          // We don't need to do any further validation on non-pointer 
          // objects.
          lua_pop(lua, 1);
          continue;
        }

        // No tables of tables allowed!
        interpreter_storage_t* tvar = (interpreter_storage_t*)tval;
        if (tvar->type == INTERPRETER_TABLE)
        {
          if (preexisting_var)
            skip_this_var = true;
          else
            polymec_error("Value error: Key '%s' in table %s stores a table (not allowed!)", tkey, key); 
        }

        // Make sure the type of this key is the same as the others.
        if (!skip_this_var)
        {
          if (value_type == INTERPRETER_TERMINUS)
            value_type = tvar->type;
          else if (tvar->type != value_type)
          {
            polymec_error("Value error: Key '%s' in table %s stores a %s (should be %s)", 
                tkey, key, type_names[tvar->type], type_names[value_type]);
          }
        }

        // Removes value from stack.
        lua_pop(lua, 1);
      }

      // If we need to skip the table, do so.
      if (skip_this_var)
      {
        lua_pop(lua, 1);
        continue;
      }

      // Now traverse it again and update our unordered map.
      str_ptr_unordered_map_t* table = str_ptr_unordered_map_new();
      lua_pushnil(lua);
      while (lua_next(lua, -2))
      {
        // Key is at index -2, value is at -1.
        static const int key_index = -2;
        static const int val_index = -1;
        char* tkey = (char*)lua_tostring(lua, key_index);
        if (lua_isnumber(lua, val_index))
        {
          double* var = malloc(sizeof(double));
          *var = lua_tonumber(lua, val_index);
          str_ptr_unordered_map_insert_with_dtor(table, tkey, var, destroy_table_entry);
        }
        else if (lua_isstring(lua, val_index))
        {
          const char* var = lua_tostring(lua, val_index);
          str_ptr_unordered_map_insert_with_dtor(table, tkey, strdup(var), destroy_table_entry);
        }
        else if (lua_islightuserdata(lua, val_index))
        {
          void* tval = (void*)lua_topointer(lua, val_index);
          interpreter_storage_t* tvar = (interpreter_storage_t*)tval;
          str_ptr_unordered_map_insert_with_dtor(table, tkey, tvar->datum, destroy_table_entry);
        }

        // Removes value from stack.
        lua_pop(lua, 1);
      }
      var = store_table(table);
    }
    else 
    {
      // We're out of types, so we hope this one's a light user data.
      if (!lua_islightuserdata(lua, val_index))
      {
        // Skip this variable if it existed before and is the wrong type.
        if (preexisting_var)
        {
          lua_pop(lua, 1);
          continue;
        }
      }
      var = (void*)lua_topointer(lua, val_index);
    }
    interpreter_map_insert_with_dtor(interp->store, strdup(key), var, destroy_variable);

    // Removes value from stack -- key is kept for next iteration.
    lua_pop(lua, 1);
  }
  lua_pop(lua, 1); // Pops the globals table.
}

void interpreter_parse_string(interpreter_t* interp, char* input_string)
{
  // Clear the current data store.
  interpreter_map_clear(interp->store);

  // Initialize the Lua interpreter.
  lua_State* lua = interpreter_open_lua(interp);

printf("%s\n", input_string);
  // Load the input, but does not parse it.
  int error = luaL_loadstring(lua, input_string);
  if (error == LUA_ERRSYNTAX)
    polymec_error("in parsing input syntax:\n%s", lua_tostring(lua, 1));
  else if (error != LUA_OK)
    polymec_error("Internal error in interpreter.");

  // Parses the input.
  error = lua_pcall(lua, 0, LUA_MULTRET, 0);
  if (error == LUA_ERRRUN)
    polymec_error("in parsing input: %s", lua_tostring(lua, 1));
  else if (error != LUA_OK)
    polymec_error("Internal error in interpreter.");

  // Move the data from the chunk to the data store.
  interpreter_store_chunk_contents(interp, lua);

  // Put lua away.
  interpreter_close_lua(interp, lua);
}

void interpreter_parse_file(interpreter_t* interp, char* input_file)
{
  log_detail("interpreter: Looking for input in file '%s'...", input_file);
  FILE* desc = fopen(input_file, "r");
  if (desc == NULL)
    polymec_error("interpreter: Could not open input file '%s'", input_file);
  fclose(desc);

  // Clear the current data store.
  interpreter_map_clear(interp->store);

  // Initialize the Lua interpreter.
  lua_State* lua = interpreter_open_lua(interp);

  // Load the input from the file, but does not parse it.
  int error = luaL_loadfile(lua, input_file);
  if (error == LUA_ERRSYNTAX)
    polymec_error("in parsing input syntax:\n%s", lua_tostring(lua, 1));
  else if (error != LUA_OK)
    polymec_error("Internal error in interpreter.");

  // Parses the input.
  error = lua_pcall(lua, 0, LUA_MULTRET, 0);
  if (error == LUA_ERRRUN)
    polymec_error("in parsing input: %s", lua_tostring(lua, 1));
  else if (error != LUA_OK)
    polymec_error("Internal error in interpreter.");

  // Move the data from the chunk to the data store.
  interpreter_store_chunk_contents(interp, lua);

  // Put lua away.
  interpreter_close_lua(interp, lua);
}

bool interpreter_contains(interpreter_t* interp, const char* variable, interpreter_var_type_t type)
{
  ASSERT(type != INTERPRETER_TERMINUS);
  interpreter_storage_t** storage = interpreter_map_get(interp->store, (char*)variable);
  if (storage == NULL)
    return false;
  return ((*storage)->type == type);
}

char* interpreter_get_string(interpreter_t* interp, const char* name)
{
  interpreter_storage_t** storage = interpreter_map_get(interp->store, (char*)name);
  if (storage == NULL)
    return NULL;
  if ((*storage)->type != INTERPRETER_STRING)
    return NULL;
  return (char*)((*storage)->datum);
}

double interpreter_get_number(interpreter_t* interp, const char* name)
{
  interpreter_storage_t** storage = interpreter_map_get(interp->store, (char*)name);
  if (storage == NULL)
    return -FLT_MAX;
  if ((*storage)->type != INTERPRETER_NUMBER)
    return -FLT_MAX;
  return *((double*)(*storage)->datum);
}

mesh_t* interpreter_get_mesh(interpreter_t* interp, const char* name)
{
  interpreter_storage_t** storage = interpreter_map_get(interp->store, (char*)name);
  if (storage == NULL)
    return NULL;
  if ((*storage)->type != INTERPRETER_MESH)
    return NULL;
  (*storage)->dtor = NULL; // Caller assumes responsibility for mesh.
  return (mesh_t*)((*storage)->datum);
}

st_func_t* interpreter_get_function(interpreter_t* interp, const char* name)
{
  interpreter_storage_t** storage = interpreter_map_get(interp->store, (char*)name);
  if (storage == NULL)
    return NULL;
  if ((*storage)->type != INTERPRETER_FUNCTION)
    return NULL;
  return (st_func_t*)((*storage)->datum);
}

str_ptr_unordered_map_t* interpreter_get_table(interpreter_t* interp, const char* name)
{
  interpreter_storage_t** storage = interpreter_map_get(interp->store, (char*)name);
  if (storage == NULL)
    return NULL;
  if ((*storage)->type != INTERPRETER_TABLE)
    return NULL;
  (*storage)->dtor = NULL; // Caller assumes responsibility for table.
  return (str_ptr_unordered_map_t*)((*storage)->datum);
}

void* interpreter_get_user_defined(interpreter_t* interp, const char* name)
{
  interpreter_storage_t** storage = interpreter_map_get(interp->store, (char*)name);
  if (storage == NULL)
    return NULL;
  if ((*storage)->type != INTERPRETER_USER_DEFINED)
    return NULL;
  (*storage)->dtor = NULL; // Caller assumes responsibility for user-defined datum.
  return (*storage)->datum;
}

//------------------------------------------------------------------------
//                          Lua helpers.
//------------------------------------------------------------------------

bool lua_isstfunc(struct lua_State* lua, int index)
{
  if (!lua_islightuserdata(lua, index))
    return false;
  interpreter_storage_t* storage = (interpreter_storage_t*)lua_topointer(lua, index);
  return (storage->type == INTERPRETER_FUNCTION);
}

st_func_t* lua_tostfunc(struct lua_State* lua, int index)
{
  if (!lua_islightuserdata(lua, index))
    return NULL;
  interpreter_storage_t* storage = (interpreter_storage_t*)lua_topointer(lua, index);
  if (storage->type == INTERPRETER_FUNCTION)
    return (st_func_t*)storage->datum;
  else
    return NULL;
}

void lua_pushstfunc(struct lua_State* lua, st_func_t* func)
{
  // Bundle it up and store it in the given variable.
  interpreter_storage_t* storage = malloc(sizeof(interpreter_storage_t));
  storage->type = INTERPRETER_FUNCTION;
  storage->datum = (void*)func;
  storage->dtor = NULL;
  lua_pushlightuserdata(lua, (void*)storage);
}

bool lua_ismesh(struct lua_State* lua, int index)
{
  if (!lua_islightuserdata(lua, index))
    return false;
  interpreter_storage_t* storage = (interpreter_storage_t*)lua_topointer(lua, index);
  return (storage->type == INTERPRETER_MESH);
}

mesh_t* lua_tomesh(struct lua_State* lua, int index)
{
  if (!lua_islightuserdata(lua, index))
    return NULL;
  interpreter_storage_t* storage = (interpreter_storage_t*)lua_topointer(lua, index);
  if (storage->type == INTERPRETER_MESH)
    return (mesh_t*)storage->datum;
  else
    return NULL;
}

void lua_pushmesh(struct lua_State* lua, mesh_t* mesh)
{
  // Bundle it up and store it in the given variable.
  interpreter_storage_t* storage = malloc(sizeof(interpreter_storage_t));
  storage->type = INTERPRETER_MESH;
  storage->datum = (void*)mesh;
  storage->dtor = destroy_mesh;
  lua_pushlightuserdata(lua, (void*)storage);
}

bool lua_isuserdefined(struct lua_State* lua, int index)
{
  if (!lua_islightuserdata(lua, index))
    return false;
  interpreter_storage_t* storage = (interpreter_storage_t*)lua_topointer(lua, index);
  return (storage->type == INTERPRETER_USER_DEFINED);
}

void* lua_touserdefined(struct lua_State* lua, int index)
{
  if (!lua_islightuserdata(lua, index))
    return NULL;
  interpreter_storage_t* storage = (interpreter_storage_t*)lua_topointer(lua, index);
  if (storage->type == INTERPRETER_USER_DEFINED)
    return (void*)storage->datum;
  else
    return NULL;
}

void lua_pushuserdefined(struct lua_State* lua, void* userdefined, void (*dtor)(void*))
{
  // Bundle it up and store it in the given variable.
  interpreter_storage_t* storage = malloc(sizeof(interpreter_storage_t));
  storage->type = INTERPRETER_USER_DEFINED;
  storage->datum = (void*)userdefined;
  storage->dtor = dtor;
  lua_pushlightuserdata(lua, (void*)storage);
}

#ifdef __cplusplus
}
#endif

