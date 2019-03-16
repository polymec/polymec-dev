// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_LUA_TYPE_H
#define POLYMEC_LUA_TYPE_H

#include "core/polymec.h"

// This is a forward declaration of the Lua interpreter state, needed to
// implement functions in Lua.
typedef struct lua_State lua_State;

/// \addtogroup core core
///@{

/// \addtogroup lua lua
/// Functions for exposing types to the Lua interpreter.
///@{

/// \struct lua_module_field
/// This type represents a field in a Lua module, with a name,
/// a getter, and a setter (if any). A field must have a getter
/// but doesn't need a setter if it is read only.
typedef struct
{
  /// The name of the field in the lua module.
  const char* name;
  /// The function for retrieving this lua module field.
  int (*getter)(lua_State* L);
  /// The function for setting the value of the lua module field.
  int (*setter)(lua_State* L);
} lua_module_field;

/// \struct lua_module_function
/// This type represents a function for a Lua module.
typedef struct
{
  /// The name of the function in the lua module.
  const char* name;
  /// The C function implementing the lua module function.
  int (*func)(lua_State* L);
  /// A string documenting the lua module function.
  const char* doc;
} lua_module_function;

/// Registers a set of fields and functions with the interpreter L in the
/// module with the given name.
/// \param [in] L The Lua interpreter in which the module is registered.
/// \param [in] module_name The name of the module being registered.
/// \param [in] module_doc A string containing documentation for the module.
/// \param [in] fields An array of fields (get/set functions) that expose data within
///                    the module.
/// \param [in] functions An array of module functions.
void lua_register_module(lua_State* L,
                         const char* module_name,
                         const char* module_doc,
                         lua_module_field fields[],
                         lua_module_function functions[]);

/// Registers a set of functions in a named table within a module. Useful for
/// factory functions.
/// \param [in] L The Lua interpreter in which the module functions are registered.
/// \param [in] module_name The module in which the functions are registered.
/// \param [in] table_name The table within which the functions are registered.
/// \param [in] table_doc A string containing documentation for the module function table.
/// \param [in] funcs An array of module functions to be registered.
void lua_register_module_function_table(lua_State* L,
                                        const char* module_name,
                                        const char* table_name,
                                        const char* table_doc,
                                        lua_module_function funcs[]);

/// \enum lua_ownership_t
/// This type specifies whether a Lua object is owned by the Lua
/// environment or by some C component.
typedef enum
{
  LUA_OWNED_BY_LUA, // owned by Lua environment
  LUA_OWNED_BY_C,   // owned by C environment
} lua_ownership_t;

/// \struct lua_class_field
/// This type represents a field in a Lua class, with a name,
/// a getter, and a setter (if any). A field must have a getter
/// but doesn't need a setter if it is read only.
typedef struct
{
  /// The name of the field belonging to a lua class.
  const char* name;
  /// The function that retrieves the lua class field.
  int (*getter)(lua_State* L);
  /// The function that sets the value of the lua class field.
  int (*setter)(lua_State* L);
} lua_class_field;

/// \struct lua_class_method
/// This type represents a method for a Lua class.
typedef struct
{
  /// The name of a method in a lua class.
  const char* name;
  /// The C function that implements the lua class method.
  int (*method)(lua_State* L);
  /// A string documenting the lua class method.
  const char* doc;
} lua_class_method;

/// Registers a new Lua class with the interpreter L, giving it a name,
/// a set of static functions and a set of methods. The functions live in
/// a module named after the type.
/// \param [in] L The Lua interpreter in which the Lua class is registered.
/// \param [in] class_name The name of the Lua class. Must not correspond to an already-
///                        registered Lua class.
/// \param [in] class_doc A string containing documentation for the class.
/// \param [in] functions An array of functions exposed within the module named after
///                       this class. Can be NULL if no module functions are desired.
/// \param [in] fields An array of accessible fields with get/set functions, available
///                    on each object of this class's type. Can be NULL if these objects
///                    have no fields.
/// \param [in] methods An array of methods callable on each object of this class's type.
///                     Can be NULL if these objects have no methods.
/// \param [in] c_dtor A C destructor function that frees the data in the context pointer
///                    provided in the call to \ref lua_push_object that creates an object
///                    of this class's type.
void lua_register_class(lua_State* L,
                        const char* class_name,
                        const char* class_doc,
                        lua_module_function functions[],
                        lua_class_field fields[],
                        lua_class_method methods[],
                        void (*c_dtor)(void* context));

/// Pushes a new (polymec) Lua object of the given class to the top of the
/// stack in the interpreter L, associating it with a context pointer and a
/// destructor to be called when the object is garbage-collected.
/// \param [in] L The Lua interpreter in which an object of the given type is to be created.
/// \param [in] class_name The name of the class identifying the type of object to create.
/// \param [in] context A pointer holding the initial state of the created object.
void lua_push_object(lua_State* L,
                     const char* class_name,
                     void* context);

/// Returns true if the object at the given index in the interpreter is of
/// the type identified by the given class name, false if not.
/// \param [in] L The Lua interpreter in which the object is queried.
/// \param [in] index The index of the object being queried within the interpreter L.
/// \param [in] class_name The class identifying the type in question.
bool lua_is_object(lua_State* L,
                   int index,
                   const char* class_name);

/// Returns the context pointer associated with the (polymec) lua object of the
/// given type at the given index in the interpreter L, or NULL if the value at
/// that index is not of that type.
/// \param [in] L The Lua interpreter from which the object is retrieved.
/// \param [in] index The index of the object being retrieved within the interpreter L.
/// \param [in] class_name The class identifying the type of the object being retrieved.
void* lua_to_object(lua_State* L,
                    int index,
                    const char* class_name);

/// This is a wrapper around lua_to_object that throws an error if the value at
/// the given index is not an object of the correct class. Analogous to
/// luaL_checkudata.
/// \param [in] L The Lua interpreter in which the object is checked.
/// \param [in] index The index of the object being checked within the interpreter L.
/// \param [in] class_name The class identifying the type of the object being checked.
void* lua_check_object(lua_State* L,
                       int index,
                       const char* class_name);

/// This function transfers the ownership of the object with the given type
/// to C or to Lua. By default, Lua owns all objects pushed to the stack.
/// \param [in] L The Lua interpreter to or from which the object's is transferred.
/// \param [in] index The index of the object being transferred within the interpreter L.
/// \param [in] class_name The class identifying the type of the object being transferred.
/// \param [in] ownership The direction of transfer of ownership for the object.
void lua_transfer_object(lua_State* L,
                         int index,
                         const char* class_name,
                         lua_ownership_t ownership);

///@}

///@}

#endif

