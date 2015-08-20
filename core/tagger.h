// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_TAGGER_H
#define POLYMEC_TAGGER_H

#include <stdlib.h>
#include <stdbool.h>
#include "core/serializer.h"

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

// Resizes the tag with the given name if it exists in the tagger, preserving 
// existing data. If the given tag does not exist, this function does nothing.
void tagger_resize_tag(tagger_t* tagger, const char* tag, int new_num_indices);

// Associates a "property" (of unspecified size and type) with the given tag 
// within the tagger. A serializer should be given so that properties can be
// copied automatically between tags.
bool tagger_set_property(tagger_t* tagger, 
                         const char* tag, 
                         const char* property, 
                         void* data, 
                         serializer_t* serializer);

// Returns a pointer to the data for the property of the given name 
// associated with the given tag. Returns NULL if the property or tag are 
// not found.
void* tagger_property(tagger_t* tagger, const char* tag, const char* property);

// Deletes the property from the tag, calling its destructor. Has no effect 
// if the property and/or tag are not found.
void tagger_delete_property(tagger_t* tagger, const char* tag, const char* property);

// Allows the traversal of properties for a tag. Set *pos to 0 to reset the 
// iteration.
bool tagger_next_property(tagger_t* tagger, const char* tag, int* pos, 
                          char** prop_name, void** prop_data, 
                          serializer_t** prop_serializer);

// Renames a tag within the tagger, preserving all information and associations.
void tagger_rename_tag(tagger_t* tagger, const char* old_tag, const char* new_tag);

// Deletes a tag from the tagger, destroying all of its information (including properties).
void tagger_delete_tag(tagger_t* tagger, const char* tag);

// Copies the indices of the given "source" tag to a "destination" tag. If 
// the "source" tag does not exist, this function has no effect. Otherwise, 
// the destination is overwritten. NOTE: No properties are copied.
void tagger_copy_tag(tagger_t* tagger, const char* source_tag, const char* destination_tag);

// Replaces the given tag with its union with the other tag. If the other 
// tag does not exist, this function has no effect. Otherwise, the given 
// tag will be created if it does not yet exist. 
void tagger_unite_tag(tagger_t* tagger, const char* tag, const char* other);

// Replaces the given tag with its intersection with the other tag. If either 
// of the given tags does not exist in this tagger, this function has no effect.
void tagger_intersect_tag(tagger_t* tagger, const char* tag, const char* other);

// Replaces the given tag with its difference with the other tag (i.e. tag 
// will contain no indices that appear in other after this operation). If 
// either of the given tags does not exist in this tagger, this function has no 
// effect.
void tagger_difference_tag(tagger_t* tagger, const char* tag, const char* other);

// Allows iteration over the tags in the tagger, returning the name, indices, 
// and size of the tag. Returns true if a tag was found, or false if the 
// iteration has ended. Set pos to 0 to reset the iteration.
bool tagger_next_tag(tagger_t* tagger, int* pos, char** tag_name, int** tag_indices, int* tag_size);

// Returns a serializer object that can read/write taggers from/to byte arrays.
serializer_t* tagger_serializer();

#endif
