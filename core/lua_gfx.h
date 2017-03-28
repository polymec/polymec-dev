// Copyright (c) 2012-2017, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_LUA_GFX_H
#define POLYMEC_LUA_GFX_H

#include "core/polymec.h"
#include "core/lua_types.h"
#include "core/gfx.h"

// Pushes a gfx figure onto L's stack.
void lua_push_gfx_figure(lua_State* L, gfx_figure_t* fig);

// Returns true if the item at the given index on L's stack is a gfx figure, 
// false if not.
bool lua_is_gfx_figure(lua_State* L, int index);

// Returns the gfx figure at the given index on L's stack, or NULL if the item 
// there is not a gfx figure.
gfx_figure_t* lua_to_gfx_figure(lua_State* L, int index);

// Pushes a gfx page onto L's stack.
void lua_push_gfx_page(lua_State* L, gfx_page_t* page);

// Returns true if the item at the given index on L's stack is a gfx page, 
// false if not.
bool lua_is_gfx_page(lua_State* L, int index);

// Returns the gfx page at the given index on L's stack, or NULL if the item 
// there is not a gfx page.
gfx_page_t* lua_to_gfx_page(lua_State* L, int index);

#endif

