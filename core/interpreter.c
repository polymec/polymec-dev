#include "core/interpreter.h"
#include "core/constant_st_func.h"
#include "core/unordered_set.h"
#include "core/boundary_cell_map.h"

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
  int size;                    // The size of the datum (if any).
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

static interpreter_storage_t* store_point(point_t* point)
{
  interpreter_storage_t* storage = malloc(sizeof(interpreter_storage_t));
  storage->datum = point;
  storage->type = INTERPRETER_POINT;
  storage->dtor = NULL;
  return storage;
}

static interpreter_storage_t* store_pointlist(point_t* points, int size)
{
  interpreter_storage_t* storage = malloc(sizeof(interpreter_storage_t));
  storage->datum = points;
  storage->size = size;
  storage->type = INTERPRETER_POINT_LIST;
  storage->dtor = free;
  return storage;
}

static interpreter_storage_t* store_vector(vector_t* vec)
{
  interpreter_storage_t* storage = malloc(sizeof(interpreter_storage_t));
  storage->datum = vec;
  storage->type = INTERPRETER_VECTOR;
  storage->dtor = NULL;
  return storage;
}

static interpreter_storage_t* store_vectorlist(vector_t* vectors, int size)
{
  interpreter_storage_t* storage = malloc(sizeof(interpreter_storage_t));
  storage->datum = vectors;
  storage->size = size;
  storage->type = INTERPRETER_VECTOR_LIST;
  storage->dtor = free;
  return storage;
}

static interpreter_storage_t* store_boundingbox(bbox_t* var)
{
  interpreter_storage_t* storage = malloc(sizeof(interpreter_storage_t));
  storage->type = INTERPRETER_BOUNDING_BOX;
  storage->datum = (void*)var;
  storage->dtor = NULL;
  return storage;
}

static void destroy_mesh(void* mesh)
{
  mesh_free((mesh_t*)mesh);
}

static interpreter_storage_t* store_mesh(mesh_t* var)
{
  interpreter_storage_t* storage = malloc(sizeof(interpreter_storage_t));
  storage->type = INTERPRETER_MESH;
  storage->datum = (void*)var;
  storage->dtor = destroy_mesh;
  return storage;
}

static interpreter_storage_t* store_scalar_function(st_func_t* var)
{
  ASSERT(st_func_num_comp(var) == 1);
  interpreter_storage_t* storage = malloc(sizeof(interpreter_storage_t));
  storage->type = INTERPRETER_SCALAR_FUNCTION;
  storage->datum = (void*)var;
  storage->dtor = NULL;
  return storage;
}

static interpreter_storage_t* store_vector_function(st_func_t* var)
{
  ASSERT(st_func_num_comp(var) == 3);
  interpreter_storage_t* storage = malloc(sizeof(interpreter_storage_t));
  storage->type = INTERPRETER_VECTOR_FUNCTION;
  storage->datum = (void*)var;
  storage->dtor = NULL;
  return storage;
}

static void destroy_table(void* table)
{
  str_ptr_unordered_map_free((str_ptr_unordered_map_t*)table);
}

static interpreter_storage_t* store_sequence(double* sequence, int len)
{
  interpreter_storage_t* storage = malloc(sizeof(interpreter_storage_t));
  storage->datum = sequence;
  storage->type = INTERPRETER_SEQUENCE;
  storage->dtor = free;
  storage->size = len;
  return storage;
}

static interpreter_storage_t* store_table(str_ptr_unordered_map_t* table)
{
  interpreter_storage_t* storage = malloc(sizeof(interpreter_storage_t));
  storage->datum = table;
  storage->type = INTERPRETER_TABLE;
  storage->dtor = destroy_table;
  return storage;
}

static interpreter_storage_t* store_user_defined(void* user_defined, void (*dtor)(void*))
{
  interpreter_storage_t* storage = malloc(sizeof(interpreter_storage_t));
  storage->type = INTERPRETER_USER_DEFINED;
  storage->datum = (void*)user_defined;
  storage->dtor = dtor;
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

// Creates a bounding box from a table.
static int bounding_box(lua_State* lua)
{
  // Check the arguments.
  int num_args = lua_gettop(lua);
  if (num_args != 1)
  {
    if (!lua_istable(lua, 1))
    {
      lua_pushstring(lua, "Argument must be a table containing x1, x2, y1, y2, z1, z2 values.");
      lua_error(lua);
      return LUA_ERRRUN;
    }
  }

  // Look for x1, x2, y1, y2, z1, z2 in the table.
  bbox_t* bbox = bbox_new(0.0, 1.0, 0.0, 1.0, 0.0, 1.0);
  const char* entries[] = {"x1", "x2", "y1", "y2", "z1", "z2"};
  for (int i = 1; i <= 6; ++i)
  {
    lua_pushstring(lua, entries[i]);
    lua_gettable(lua, 1); // Reads name from top, replaces with bounds[name].
    if (!lua_isnumber(lua, -1))
    {
      lua_pushstring(lua, "x1, x2, y1, y2, z1, z2, must all be numbers.");
      lua_error(lua);
      return LUA_ERRRUN;
    }
    switch(i)
    {
      case 0: bbox->x1 = lua_tonumber(lua, -1);
              break;
      case 1: bbox->x2 = lua_tonumber(lua, -1);
              break;
      case 2: bbox->y1 = lua_tonumber(lua, -1);
              break;
      case 3: bbox->y2 = lua_tonumber(lua, -1);
              break;
      case 4: bbox->z1 = lua_tonumber(lua, -1);
              break;
      case 5: bbox->z2 = lua_tonumber(lua, -1);
              break;
      default: break;
    }
    lua_pop(lua, 1); 
  }

  // Push the bounding box onto the stack.
  lua_pushboundingbox(lua, bbox);
  return 1;
}

// Creates a constant function from a number or a 3-tuple.
static int constant_function(lua_State* lua)
{
  // Check the argument.
  int num_args = lua_gettop(lua);
  if (num_args == 1) // Scalar-valued constant.
  {
    if (!lua_isnumber(lua, 1))
    {
      lua_pushstring(lua, "Argument must be a number.");
      lua_error(lua);
      return LUA_ERRRUN;
    }
  }
  else if (num_args == 3) // Vector-valued constant.
  {
    for (int i = 0; i < 3; ++i)
    {
      if (!lua_isnumber(lua, i+1))
      {
        lua_pushfstring(lua, "Argument %d must be a number.", i);
        lua_error(lua);
        return LUA_ERRRUN;
      }
    }
  }
  else
  {
    lua_pushstring(lua, "Argument must be a 1 or 3 numbers.");
    lua_error(lua);
    return LUA_ERRRUN;
  }

  // Get the arguments.
  double args[3];
  for (int i = 0; i < num_args; ++i)
    args[i] = lua_tonumber(lua, i+1);

  // Push a constant function onto the stack.
  st_func_t* func = constant_st_func_new(num_args, args);
  if (num_args == 1)
    lua_pushscalarfunction(lua, func);
  else
    lua_pushvectorfunction(lua, func);
  return 1;
}

// Creates a vector-valued function from 3 scalar functions.
static int vector_function(lua_State* lua)
{
  // Check the argument.
  int num_args = lua_gettop(lua);
  if (num_args != 3)
  {
    lua_pushstring(lua, "Arguments must be 3 scalar functions.");
    lua_error(lua);
    return LUA_ERRRUN;
  }
  for (int i = 1; i <= 3; ++i)
  {
    if (!lua_isscalarfunction(lua, i))
    {
      lua_pushfstring(lua, "Argument %d must be a scalar function.", i);
      lua_error(lua);
      return LUA_ERRRUN;
    }
  }

  st_func_t* functions[3];
  for (int i = 1; i <= 3; ++i)
    functions[i-1] = lua_toscalarfunction(lua, i);
  lua_pushvectorfunction(lua, multicomp_st_func_from_funcs("vector function", functions, 3));
  return 1;
}

// Creates a periodic boundary condition from a pair of tags.
static int periodic_bc(lua_State* lua)
{
  // Check the argument.
  int num_args = lua_gettop(lua);
  if (num_args != 2)
  {
    lua_pushstring(lua, "Arguments must be 2 boundary mesh (face) tags.");
    lua_error(lua);
    return LUA_ERRRUN;
  }
  for (int i = 1; i <= 3; ++i)
  {
    if (!lua_isstring(lua, i))
    {
      lua_pushfstring(lua, "Argument %d must be a face tag.", i);
      lua_error(lua);
      return LUA_ERRRUN;
    }
  }

  const char* tag1 = lua_tostring(lua, 1);
  const char* tag2 = lua_tostring(lua, 2);
  periodic_bc_t* bc = periodic_bc_new(tag1, tag2);
  lua_pushuserdefined(lua, bc, NULL);
  return 1;
}

static void register_default_functions(interpreter_t* interp)
{
  interpreter_register_function(interp, "bounding_box", bounding_box);
  interpreter_register_function(interp, "constant_function", constant_function);
  interpreter_register_function(interp, "vector_function", vector_function);
  interpreter_register_function(interp, "periodic_bc", periodic_bc);
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
      else if (entry->type == INTERPRETER_POINT)
      {
        if (!lua_ispoint(lua, val_index))
        {
          if (preexisting_var)
            skip_this_var = true;
          else
            polymec_error("Type error: %s must be a point.", key);
        }
      }
      else if (entry->type == INTERPRETER_POINT_LIST)
      {
        if (!lua_ispointlist(lua, val_index))
        {
          if (preexisting_var)
            skip_this_var = true;
          else
            polymec_error("Type error: %s must be a list of points.", key);
        }
      }
      else if (entry->type == INTERPRETER_VECTOR)
      {
        if (!lua_isvector(lua, val_index))
        {
          if (preexisting_var)
            skip_this_var = true;
          else
            polymec_error("Type error: %s must be a vector.", key);
        }
      }
      else if (entry->type == INTERPRETER_VECTOR_LIST)
      {
        if (!lua_isvectorlist(lua, val_index))
        {
          if (preexisting_var)
            skip_this_var = true;
          else
            polymec_error("Type error: %s must be a list of vectors.", key);
        }
      }
      else if (entry->type == INTERPRETER_BOUNDING_BOX)
      {
        if (!lua_isboundingbox(lua, val_index))
        {
          if (preexisting_var)
            skip_this_var = true;
          else
            polymec_error("Type error: %s must be a bounding box.", key);
        }
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
      else if (entry->type == INTERPRETER_SCALAR_FUNCTION)
      {
        if (!lua_islightuserdata(lua, val_index))
        {
          if (preexisting_var)
            skip_this_var = true;
          else
            polymec_error("Type error: %s must be a scalar-valued function.", key);
        }
        interpreter_storage_t* var = (void*)lua_topointer(lua, val_index);
        if (var->type != INTERPRETER_SCALAR_FUNCTION)
        {
          if (preexisting_var)
            skip_this_var = true;
          else
            polymec_error("Type error: %s must be a scalar-valued function.", key);
        }
      }
      else if (entry->type == INTERPRETER_VECTOR_FUNCTION)
      {
        if (!lua_islightuserdata(lua, val_index))
        {
          if (preexisting_var)
            skip_this_var = true;
          else
            polymec_error("Type error: %s must be a vector-valued function.", key);
        }
        interpreter_storage_t* var = (void*)lua_topointer(lua, val_index);
        if (var->type != INTERPRETER_VECTOR_FUNCTION)
        {
          if (preexisting_var)
            skip_this_var = true;
          else
            polymec_error("Type error: %s must be a vector-valued function.", key);
        }
      }
      else if ((entry->type == INTERPRETER_SEQUENCE) && !lua_issequence(lua, val_index))
      {
        if (preexisting_var)
          skip_this_var = true;
        else
          polymec_error("Type error: %s must be a sequence of numbers.", key);
      }
      else if ((entry->type == INTERPRETER_TABLE) && !lua_istable(lua, val_index))
      {
        if (preexisting_var)
          skip_this_var = true;
        else
          polymec_error("Type error: %s must be a table mapping strings to objects.", key);
      }

      // Skip empty tables, as they cause problems.
      else if (lua_istable(lua, val_index) && (lua_rawlen(lua, val_index) == 0))
      {
        skip_this_var = true;
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
    else if (lua_issequence(lua, val_index))
    {
      // Sequences can be interpreted in many ways. Are we asked to 
      // interpret this sequence as something other than just a sequence?
      int len;
      double* seq = lua_tosequence(lua, val_index, &len);
      if ((entry == NULL) || (entry->type == INTERPRETER_SEQUENCE))
      {
        var = store_sequence(seq, len);
      }
      else if (entry->type == INTERPRETER_POINT)
      {
        point_t* p = point_new(seq[0], seq[1], seq[2]);
        var = store_point(p);
      }
      else if (entry->type == INTERPRETER_VECTOR)
      {
        vector_t* v = vector_new(seq[0], seq[1], seq[2]);
        var = store_vector(v);
      }

      // Clean up if necessary.
      if ((entry != NULL) && (entry->type != INTERPRETER_SEQUENCE))
        free(seq);
    }
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
            polymec_error("Type error: %s must be a sequence or a table mapping strings to objects.", key);
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

      // If we need to skip the table/sequence, do so.
      if (skip_this_var)
      {
        lua_pop(lua, 1);
        continue;
      }

      // Now traverse the table again and create a C data structure.
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
          str_ptr_unordered_map_insert_with_kv_dtor(table, tkey, var, destroy_table_entry);
        }
        else if (lua_isstring(lua, val_index))
        {
          const char* var = lua_tostring(lua, val_index);
          str_ptr_unordered_map_insert_with_kv_dtor(table, tkey, strdup(var), destroy_table_entry);
        }
        else if (lua_islightuserdata(lua, val_index))
        {
          void* tval = (void*)lua_topointer(lua, val_index);
          interpreter_storage_t* tvar = (interpreter_storage_t*)tval;
          str_ptr_unordered_map_insert_with_kv_dtor(table, tkey, tvar->datum, destroy_table_entry);
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
    interpreter_map_insert_with_kv_dtor(interp->store, strdup(key), var, destroy_variable);

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

void interpreter_set_string(interpreter_t* interp, const char* name, const char* value)
{
  interpreter_storage_t* storage = store_string(value);
  interpreter_map_insert_with_kv_dtor(interp->store, strdup(name), storage, destroy_variable);
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

void interpreter_set_number(interpreter_t* interp, const char* name, double value)
{
  interpreter_storage_t* storage = store_number(value);
  interpreter_map_insert_with_kv_dtor(interp->store, strdup(name), storage, destroy_variable);
}

point_t* interpreter_get_point(interpreter_t* interp, const char* name)
{
  interpreter_storage_t** storage = interpreter_map_get(interp->store, (char*)name);
  if (storage == NULL)
    return NULL;
  if ((*storage)->type != INTERPRETER_POINT)
    return NULL;
  return (point_t*)((*storage)->datum);
}

void interpreter_set_point(interpreter_t* interp, const char* name, point_t* value)
{
  interpreter_storage_t* storage = store_point(value);
  interpreter_map_insert_with_kv_dtor(interp->store, strdup(name), storage, destroy_variable);
}

point_t* interpreter_get_pointlist(interpreter_t* interp, const char* name, int* num_points)
{
  interpreter_storage_t** storage = interpreter_map_get(interp->store, (char*)name);
  if (storage == NULL)
    return NULL;
  if ((*storage)->type != INTERPRETER_POINT_LIST)
    return NULL;
  *num_points = (*storage)->size;
  return (point_t*)((*storage)->datum);
}

void interpreter_set_pointlist(interpreter_t* interp, const char* name, point_t* points, int num_points)
{
  interpreter_storage_t* storage = store_pointlist(points, num_points);
  interpreter_map_insert_with_kv_dtor(interp->store, strdup(name), storage, destroy_variable);
}

vector_t* interpreter_get_vector(interpreter_t* interp, const char* name)
{
  interpreter_storage_t** storage = interpreter_map_get(interp->store, (char*)name);
  if (storage == NULL)
    return NULL;
  if ((*storage)->type != INTERPRETER_VECTOR)
    return NULL;
  return (vector_t*)((*storage)->datum);
}

void interpreter_set_vector(interpreter_t* interp, const char* name, vector_t* value)
{
  interpreter_storage_t* storage = store_vector(value);
  interpreter_map_insert_with_kv_dtor(interp->store, strdup(name), storage, destroy_variable);
}

vector_t* interpreter_get_vectorlist(interpreter_t* interp, const char* name, int* num_vectors)
{
  interpreter_storage_t** storage = interpreter_map_get(interp->store, (char*)name);
  if (storage == NULL)
    return NULL;
  if ((*storage)->type != INTERPRETER_VECTOR_LIST)
    return NULL;
  *num_vectors = (*storage)->size;
  return (vector_t*)((*storage)->datum);
}

void interpreter_set_vectorlist(interpreter_t* interp, const char* name, vector_t* vectors, int num_vectors)
{
  interpreter_storage_t* storage = store_vectorlist(vectors, num_vectors);
  interpreter_map_insert_with_kv_dtor(interp->store, strdup(name), storage, destroy_variable);
}

bbox_t* interpreter_get_boundingbox(interpreter_t* interp, const char* name)
{
  interpreter_storage_t** storage = interpreter_map_get(interp->store, (char*)name);
  if (storage == NULL)
    return NULL;
  if ((*storage)->type != INTERPRETER_BOUNDING_BOX)
    return NULL;
  return (bbox_t*)((*storage)->datum);
}

void interpreter_set_boundingbox(interpreter_t* interp, const char* name, bbox_t* value)
{
  interpreter_storage_t* storage = store_boundingbox(value);
  interpreter_map_insert_with_kv_dtor(interp->store, strdup(name), storage, destroy_variable);
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

void interpreter_set_mesh(interpreter_t* interp, const char* name, mesh_t* value)
{
  interpreter_storage_t* storage = store_mesh(value);
  interpreter_map_insert_with_kv_dtor(interp->store, strdup(name), storage, destroy_variable);
}

st_func_t* interpreter_get_scalar_function(interpreter_t* interp, const char* name)
{
  interpreter_storage_t** storage = interpreter_map_get(interp->store, (char*)name);
  if (storage == NULL)
    return NULL;
  if ((*storage)->type != INTERPRETER_SCALAR_FUNCTION)
    return NULL;
  st_func_t* func = (st_func_t*)((*storage)->datum);
  ASSERT(st_func_num_comp(func) == 1);
  return func;
}

void interpreter_set_scalar_function(interpreter_t* interp, const char* name, st_func_t* value)
{
  interpreter_storage_t* storage = store_scalar_function(value);
  interpreter_map_insert_with_kv_dtor(interp->store, strdup(name), storage, destroy_variable);
}

st_func_t* interpreter_get_vector_function(interpreter_t* interp, const char* name)
{
  interpreter_storage_t** storage = interpreter_map_get(interp->store, (char*)name);
  if (storage == NULL)
    return NULL;
  if ((*storage)->type != INTERPRETER_VECTOR_FUNCTION)
    return NULL;
  st_func_t* func = (st_func_t*)((*storage)->datum);
  ASSERT(st_func_num_comp(func) == 3);
  return func;
}

void interpreter_set_vector_function(interpreter_t* interp, const char* name, st_func_t* value)
{
  interpreter_storage_t* storage = store_vector_function(value);
  interpreter_map_insert_with_kv_dtor(interp->store, strdup(name), storage, destroy_variable);
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

void interpreter_set_sequence(interpreter_t* interp, const char* name, double* sequence, int len)
{
  interpreter_storage_t* storage = store_sequence(sequence, len);
  interpreter_map_insert_with_kv_dtor(interp->store, strdup(name), storage, destroy_variable);
}

void interpreter_set_table(interpreter_t* interp, const char* name, str_ptr_unordered_map_t* value)
{
  interpreter_storage_t* storage = store_table(value);
  interpreter_map_insert_with_kv_dtor(interp->store, strdup(name), storage, destroy_variable);
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

void interpreter_set_user_defined(interpreter_t* interp, const char* name, void* value, void (*dtor)(void*))
{
  interpreter_storage_t* storage = store_user_defined(value, dtor);
  interpreter_map_insert_with_kv_dtor(interp->store, strdup(name), storage, destroy_variable);
}

//------------------------------------------------------------------------
//                          Lua helpers.
//------------------------------------------------------------------------

bool lua_issequence(struct lua_State* lua, int index)
{
  index = lua_absindex(lua, index);
  if (lua_istable(lua, index))
  {
    int len = lua_rawlen(lua, index);
    for (int i = 1; i <= len; ++i)
    {
      lua_pushinteger(lua, (lua_Integer)i);
      lua_gettable(lua, index);
      bool is_number = lua_isnumber(lua, -1);
      lua_pop(lua, 1);
      if (!is_number)
        return false;
    }
    return true;
  }
  if (!lua_islightuserdata(lua, index))
    return false;
  interpreter_storage_t* storage = (interpreter_storage_t*)lua_topointer(lua, index);
  return (storage->type == INTERPRETER_SEQUENCE);
}

double* lua_tosequence(struct lua_State* lua, int index, int* len)
{
  if (!lua_issequence(lua, index))
    return NULL;
  index = lua_absindex(lua, index);
  if (lua_istable(lua, index))
  {
    *len = lua_rawlen(lua, index);
    double* seq = malloc(sizeof(double)*(*len));
    for (int i = 1; i <= *len; ++i)
    {
      lua_pushinteger(lua, (lua_Integer)i);
      lua_gettable(lua, index);
      seq[i-1] = lua_tonumber(lua, -1);
      lua_pop(lua, 1);
    }
    return seq;
  }
  interpreter_storage_t* storage = (interpreter_storage_t*)lua_topointer(lua, index);
  if (storage->type == INTERPRETER_SEQUENCE)
    return (double*)storage->datum;
  else
    return NULL;
}

void lua_pushsequence(struct lua_State* lua, double* sequence, int len)
{
  // Bundle it up and store it in the given variable.
  interpreter_storage_t* storage = store_sequence(sequence, len);
  lua_pushlightuserdata(lua, (void*)storage);
}

bool lua_ispoint(struct lua_State* lua, int index)
{
  // A sequence with 3 numbers in it will work.
  if (lua_issequence(lua, index) && (lua_rawlen(lua, index) == 3))
    return true;
  if (!lua_islightuserdata(lua, index))
    return false;
  interpreter_storage_t* storage = (interpreter_storage_t*)lua_topointer(lua, index);
  return (storage->type == INTERPRETER_POINT);
}

point_t* lua_topoint(struct lua_State* lua, int index)
{
  if (!lua_ispoint(lua, index))
    return NULL;
  if (lua_istable(lua, index))
  {
    int len;
    double* seq = lua_tosequence(lua, index, &len);
    point_t* p = point_new(seq[0], seq[1], seq[2]);
    free(seq);
    return p;
  }
  interpreter_storage_t* storage = (interpreter_storage_t*)lua_topointer(lua, index);
  if (storage->type == INTERPRETER_POINT)
    return (point_t*)storage->datum;
  else if (storage->type == INTERPRETER_SEQUENCE)
  {
    double* seq = storage->datum;
    return point_new(seq[0], seq[1], seq[2]);
  }
  else
    return NULL;
}

void lua_pushpoint(struct lua_State* lua, point_t* point)
{
  // Bundle it up and store it in the given variable.
  interpreter_storage_t* storage = store_point(point);
  lua_pushlightuserdata(lua, (void*)storage);
}

bool lua_ispointlist(struct lua_State* lua, int index)
{
  if (lua_istable(lua, index)) // A table with points in it will work.
  {
    size_t len = lua_rawlen(lua, index);
    for (size_t i = 1; i <= len; ++i)
    {
      lua_pushinteger(lua, (lua_Integer)i);
      lua_gettable(lua, -1);
      bool is_point = lua_ispoint(lua, -1);
      lua_pop(lua, 1);
      if (!is_point) 
        return false;
    }
    return true;
  }
  if (!lua_islightuserdata(lua, index))
    return false;
  interpreter_storage_t* storage = (interpreter_storage_t*)lua_topointer(lua, index);
  return (storage->type == INTERPRETER_POINT_LIST);
}

point_t* lua_topointlist(struct lua_State* lua, int index, int* size)
{
  if (!lua_ispointlist(lua, index))
    return NULL;
  if (lua_istable(lua, index))
  {
    *size = (int)lua_rawlen(lua, index);
    point_t* points = malloc(sizeof(point_t) * (*size));
    for (int i = 1; i <= *size; ++i)
    {
      lua_pushinteger(lua, (lua_Integer)i);
      lua_gettable(lua, -1);
      point_t* p = lua_topoint(lua, -1);
      lua_pop(lua, 1);
      point_copy(&points[i], p);
      free(p);
    }
    return points;
  }
  interpreter_storage_t* storage = (interpreter_storage_t*)lua_topointer(lua, index);
  if (storage->type == INTERPRETER_POINT_LIST)
  {
    *size = storage->size;
    return (point_t*)storage->datum;
  }
  else
    return NULL;
}

void lua_pushpointlist(struct lua_State* lua, point_t* points, int size)
{
  // Bundle it up and store it in the given variable.
  interpreter_storage_t* storage = store_pointlist(points, size);
  lua_pushlightuserdata(lua, (void*)storage);
}

bool lua_isvector(struct lua_State* lua, int index)
{
  // A sequence with 3 numbers in it will work.
  if (lua_issequence(lua, index) && (lua_rawlen(lua, index) == 3))
    return true;
  if (!lua_islightuserdata(lua, index))
    return false;
  interpreter_storage_t* storage = (interpreter_storage_t*)lua_topointer(lua, index);
  return (storage->type == INTERPRETER_VECTOR);
}

vector_t* lua_tovector(struct lua_State* lua, int index)
{
  if (!lua_isvector(lua, index))
    return NULL;
  if (lua_istable(lua, index))
  {
    int len;
    double* seq = lua_tosequence(lua, index, &len);
    vector_t* v = vector_new(seq[0], seq[1], seq[2]);
    free(seq);
    return v;
  }
  interpreter_storage_t* storage = (interpreter_storage_t*)lua_topointer(lua, index);
  if (storage->type == INTERPRETER_VECTOR)
    return (vector_t*)storage->datum;
  else if (storage->type == INTERPRETER_SEQUENCE)
  {
    double* seq = storage->datum;
    return vector_new(seq[0], seq[1], seq[2]);
  }
  else
    return NULL;
  if (!lua_isvector(lua, index))
    return NULL;
}

void lua_pushvector(struct lua_State* lua, vector_t* vec)
{
  // Bundle it up and store it in the given variable.
  interpreter_storage_t* storage = store_vector(vec);
  lua_pushlightuserdata(lua, (void*)storage);
}

bool lua_isvectorlist(struct lua_State* lua, int index)
{
  if (lua_istable(lua, index)) // A table with vectors in it will work.
  {
    size_t len = lua_rawlen(lua, index);
    for (size_t i = 1; i <= len; ++i)
    {
      lua_pushinteger(lua, (lua_Integer)i);
      lua_gettable(lua, -1);
      bool is_vector = lua_isvector(lua, -1);
      lua_pop(lua, 1);
      if (!is_vector) 
        return false;
    }
    return true;
  }
  if (!lua_islightuserdata(lua, index))
    return false;
  interpreter_storage_t* storage = (interpreter_storage_t*)lua_topointer(lua, index);
  return (storage->type == INTERPRETER_VECTOR_LIST);
}

vector_t* lua_tovectorlist(struct lua_State* lua, int index, int* size)
{
  if (!lua_isvectorlist(lua, index))
    return NULL;
  if (lua_istable(lua, index))
  {
    *size = (int)lua_rawlen(lua, index);
    vector_t* vectors = malloc(sizeof(vector_t) * (*size));
    for (int i = 1; i <= *size; ++i)
    {
      lua_pushinteger(lua, (lua_Integer)i);
      lua_gettable(lua, -1);
      vector_t* v = lua_tovector(lua, -1);
      lua_pop(lua, 1);
      vector_copy(&vectors[i], v);
      free(v);
    }
    return vectors;
  }
  interpreter_storage_t* storage = (interpreter_storage_t*)lua_topointer(lua, index);
  if (storage->type == INTERPRETER_VECTOR_LIST)
  {
    *size = storage->size;
    return (vector_t*)storage->datum;
  }
  else
    return NULL;
}

void lua_pushvectorlist(struct lua_State* lua, vector_t* vectors, int size)
{
  // Bundle it up and store it in the given variable.
  interpreter_storage_t* storage = store_vectorlist(vectors, size);
  lua_pushlightuserdata(lua, (void*)storage);
}

bool lua_isboundingbox(struct lua_State* lua, int index)
{
  if (!lua_islightuserdata(lua, index))
    return false;
  interpreter_storage_t* storage = (interpreter_storage_t*)lua_topointer(lua, index);
  return (storage->type == INTERPRETER_BOUNDING_BOX);
}

bbox_t* lua_toboundingbox(struct lua_State* lua, int index)
{
  if (!lua_islightuserdata(lua, index))
    return NULL;
  interpreter_storage_t* storage = (interpreter_storage_t*)lua_topointer(lua, index);
  if (storage->type == INTERPRETER_BOUNDING_BOX)
    return (bbox_t*)storage->datum;
  else
    return NULL;
}

void lua_pushboundingbox(struct lua_State* lua, bbox_t* bbox)
{
  // Bundle it up and store it in the given variable.
  interpreter_storage_t* storage = store_boundingbox(bbox);
  lua_pushlightuserdata(lua, (void*)storage);
}

bool lua_isscalarfunction(struct lua_State* lua, int index)
{
  if (!lua_islightuserdata(lua, index))
    return false;
  interpreter_storage_t* storage = (interpreter_storage_t*)lua_topointer(lua, index);
  return (storage->type == INTERPRETER_SCALAR_FUNCTION);
}

st_func_t* lua_toscalarfunction(struct lua_State* lua, int index)
{
  if (!lua_islightuserdata(lua, index))
    return NULL;
  interpreter_storage_t* storage = (interpreter_storage_t*)lua_topointer(lua, index);
  if (storage->type == INTERPRETER_SCALAR_FUNCTION)
    return (st_func_t*)storage->datum;
  else
    return NULL;
}

void lua_pushscalarfunction(struct lua_State* lua, st_func_t* func)
{
  // Only single-component functions are allowed.
  ASSERT(st_func_num_comp(func) == 1); 
  // Bundle it up and store it in the given variable.
  interpreter_storage_t* storage = store_scalar_function(func);
  lua_pushlightuserdata(lua, (void*)storage);
}

bool lua_isvectorfunction(struct lua_State* lua, int index)
{
  if (!lua_islightuserdata(lua, index))
    return false;
  interpreter_storage_t* storage = (interpreter_storage_t*)lua_topointer(lua, index);
  return (storage->type == INTERPRETER_VECTOR_FUNCTION);
}

st_func_t* lua_tovectorfunction(struct lua_State* lua, int index)
{
  if (!lua_islightuserdata(lua, index))
    return NULL;
  interpreter_storage_t* storage = (interpreter_storage_t*)lua_topointer(lua, index);
  if (storage->type == INTERPRETER_VECTOR_FUNCTION)
    return (st_func_t*)storage->datum;
  else
    return NULL;
}

void lua_pushvectorfunction(struct lua_State* lua, st_func_t* func)
{
  // Only 3-component functions are allowed.
  ASSERT(st_func_num_comp(func) == 3); 
  // Bundle it up and store it in the given variable.
  interpreter_storage_t* storage = malloc(sizeof(interpreter_storage_t));
  storage->type = INTERPRETER_VECTOR_FUNCTION;
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
  interpreter_storage_t* storage = store_mesh(mesh);
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
  interpreter_storage_t* storage = store_user_defined(userdefined, dtor);
  lua_pushlightuserdata(lua, (void*)storage);
}

#ifdef __cplusplus
}
#endif

