// Copyright (c) 2012-2014, Jeffrey N. Johnson
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without 
// modification, are permitted provided that the following conditions are met:
// 
// 1. Redistributions of source code must retain the above copyright notice, this 
// list of conditions and the following disclaimer.
// 
// 2. Redistributions in binary form must reproduce the above copyright notice, 
// this list of conditions and the following disclaimer in the documentation 
// and/or other materials provided with the distribution.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE 
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#ifndef POLYMEC_TAGGER_H
#define POLYMEC_TAGGER_H

#include <stdlib.h>
#include <stdbool.h>

// A tagger is an object that stores tag and property information. A tag is 
// a set of integer indices associated with a name. It can be used to store 
// material or geometry information in the context of a spatial discretization.
typedef struct tagger_t tagger_t;

// Creates a new tagger that stores information.
tagger_t* tagger_new();

// Destroys the given tagger.
void tagger_free(tagger_t* tags);

// Copies all tags in src to dest. Does not copy property information, since 
// the sizes of properties are not stored.
void tagger_copy(tagger_t* dest, tagger_t* src);

// Creates a new tag with the given name and size within the tagger, returning 
// storage for the tag indices. If the tag already exists, this returns NULL.
int* tagger_create_tag(tagger_t* tagger, const char* tag, int num_indices);

// Returns the tag with the given name in the tagger, and stores its size 
// in num_indices. Returns NULL if the tag was not found.
int* tagger_tag(tagger_t* tagger, const char* tag, int* num_indices);

// Returns true if the tagger contains a tag with the given name, false if not.
bool tagger_has_tag(tagger_t* tagger, const char* tag);

// Associates a "property" (of unspecified size and type) with the given tag 
// within the tagger. NOTE that since we do not store size/type information, 
// properties are not copied automatically between tags.
bool tagger_set_property(tagger_t* tagger, const char* tag, const char* property, void* data, void (*destructor)(void*));

// Returns a pointer to the data for the property of the given name 
// associated with the given tag. Returns NULL if the property or tag are 
// not found.
void* tagger_property(tagger_t* tagger, const char* tag, const char* property);

// Deletes the property from the tag, calling its destructor. Has no effect 
// if the property and/or tag are not found.
void tagger_delete_property(tagger_t* tagger, const char* tag, const char* property);

// Renames a tag within the tagger, preserving all information and associations.
void tagger_rename_tag(tagger_t* tagger, const char* old_tag, const char* new_tag);

// Deletes a tag from the tagger, destroying all of its information (including properties).
void tagger_delete_tag(tagger_t* tagger, const char* tag);

// Allows iteration over the tags in the tagger, returning the name, indices, 
// and size of the tag. Returns true if a tag was found, or false if the 
// iteration has ended. Set pos to 0 to reset the iteration.
bool tagger_next_tag(tagger_t* tagger, int* pos, char** tag_name, int** tag_indices, int* tag_size);

#endif
