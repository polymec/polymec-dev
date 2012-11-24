#include "core/interpreter.h"

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

static void destroy_table_entry(char* key, void* value)
{
  free(key);
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
  // The data store.
  interpreter_map_t* store;

  // A list of valid inputs and their types.
  int num_valid_inputs;
  interpreter_validation_t* valid_inputs;
};

interpreter_t* interpreter_new(interpreter_validation_t* valid_inputs)
{
  interpreter_t* interp = malloc(sizeof(interpreter_t));

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

  return interp;
}

void interpreter_free(interpreter_t* interp)
{
  for (int i = 0; i < interp->num_valid_inputs; ++i)
    free(interp->valid_inputs[i].variable);
  free(interp->valid_inputs);
  interpreter_map_free(interp->store);
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

void interpreter_parse_string(interpreter_t* interp, char* input_string)
{
  // Clear the current data store.
  interpreter_map_clear(interp->store);

  // Initialize the Lua interpreter.
  lua_State* lua = luaL_newstate();
  ASSERT(lua != NULL);
  luaL_openlibs(lua);

  // Add some functions for creating data types.
  // FIXME
//  lua_register(lua, "cubic_mesh", lua_cubic_mesh);

  int error = luaL_dostring(lua, input_string);
  if (error == LUA_ERRSYNTAX)
    arbi_error("Syntax error in input.");
  else if (error != LUA_OK)
    arbi_error("Internal error in interpreter.");

  // Traverse Lua's table of global variables and pull them into the
  // data store.
  lua_pushnil(lua); // nil tells lua_next to start at the first key.
  while (lua_next(lua, LUA_RIDX_GLOBALS)) // Traverse globals table.
  {
    // key is at index -2, value is at index -1.
    const char* key = lua_tostring(lua, -2);

    // Does the key appear in our validation table?
    interpreter_validation_t* entry = interpreter_validation_entry(interp, key);
    if (entry != NULL)
    {
      // We must validate this variable against its allowed type.
      if ((entry->type == INTERPRETER_STRING) && !lua_isstring(lua, -2))
        arbi_error("Type error: %s must be a string.", key);
      else if ((entry->type == INTERPRETER_NUMBER) && !lua_isnumber(lua, -2))
        arbi_error("Type error: %s must be a number.", key);
      else if (entry->type == INTERPRETER_MESH)
      {
        if (!lua_islightuserdata(lua, -2))
          arbi_error("Type error: %s must be a mesh.", key);
        interpreter_storage_t* var = (void*)lua_topointer(lua, -2);
        if (var->type != INTERPRETER_MESH)
          arbi_error("Type error: %s must be a mesh.", key);
      }
      else if (entry->type == INTERPRETER_FUNCTION)
      {
        if (!lua_islightuserdata(lua, -2))
          arbi_error("Type error: %s must be a function.", key);
        interpreter_storage_t* var = (void*)lua_topointer(lua, -2);
        if (var->type != INTERPRETER_FUNCTION)
          arbi_error("Type error: %s must be a function.", key);
      }
      else if ((entry->type == INTERPRETER_TABLE) && !lua_istable(lua, -2))
        arbi_error("Type error: %s must be a table mapping strings to objects.", key);
    }

    interpreter_storage_t* var = NULL;
    if (lua_isstring(lua, -2))
      var = store_string(lua_tostring(lua, -2));
    else if (lua_isnumber(lua, -2))
      var = store_number(lua_tonumber(lua, -2));
    else if (lua_istable(lua, -2))
    {
      str_ptr_unordered_map_t* table = str_ptr_unordered_map_new();
      // Traverse this table.
      lua_pushnil(lua);
      while (lua_next(lua, -3)) // FIXME??
      {
        // Key is at index -4, value is at -5.
        if (!lua_isstring(lua, -4))
          arbi_error("Type error: %s must be a table mapping strings to objects.", key);
        if (!lua_islightuserdata(lua, -5))
          arbi_error("Type error: %s must be a table mapping strings to objects.", key);
        char* tkey = (char*)lua_tostring(lua, -4);
        void* tval = (void*)lua_topointer(lua, -5);
        interpreter_storage_t* tvar = (interpreter_storage_t*)tval;
        str_ptr_unordered_map_insert_with_dtor(table, tkey, tvar->datum, destroy_table_entry);

        // Removes value from stack.
        lua_pop(lua, 1);
      }
      var = store_table(table);
    }
    else 
    {
      ASSERT(lua_islightuserdata(lua, -2));
      var = (void*)lua_topointer(lua, -2);
    }
    interpreter_map_insert_with_dtor(interp->store, (char*)key, var, destroy_variable);

    // Removes value from stack -- key is kept for next iteration.
    lua_pop(lua, 1);
  }

  // Put lua away.
  lua_close(lua);
}

void interpreter_parse_file(interpreter_t* interp, char* input_file)
{
  log_detail("interpreter: Looking for input in file '%s'...", input_file);
  FILE* desc = fopen(input_file, "r");
  if (desc == NULL)
    arbi_error("interpreter: Could not open input file '%s'", input_file);

  // Read the contents of the file into a string.
  fseek(desc, 0, SEEK_END);
  long size = ftell(desc);
  rewind(desc);
  char* input = malloc(size*sizeof(char));
  fread(input, sizeof(char), size, desc);
  fclose(desc);

  // Parse the input script.
  interpreter_parse_string(interp, input);

  // Clean up.
  free(input);
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
  (*storage)->dtor = NULL; // Caller assumes responsibility for table.
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

#ifdef __cplusplus
}
#endif

