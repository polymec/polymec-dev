// Copyright (c) 2012-2019, Jeffrey N. Johnson
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

/// \addtogroup geometry geometry
///@{

/// \class tagger
/// A tagger is an object that stores tag information. A tag is
/// a set of integer indices associated with a name. It can be used to store
/// material or geometry information in the context of a spatial discretization.
typedef struct tagger_t tagger_t;

/// Creates a new tagger that stores information.
/// \memberof tagger
tagger_t* tagger_new(void);

/// Destroys the given tagger.
/// \memberof tagger
void tagger_free(tagger_t* tags);

/// Copies all tags in src to dest.
/// \memberof tagger
void tagger_copy(tagger_t* dest, tagger_t* src);

/// Creates a new tag with the given name and size within the tagger, returning
/// storage for the tag indices. If the tag already exists, this returns NULL.
/// \memberof tagger
int* tagger_create_tag(tagger_t* tagger, const char* tag, size_t num_indices);

/// Returns the tag with the given name in the tagger, and stores its size
/// in num_indices. Returns NULL if the tag was not found.
/// \memberof tagger
int* tagger_tag(tagger_t* tagger, const char* tag, size_t* num_indices);

/// Returns true if the tagger contains a tag with the given name, false if not.
/// \memberof tagger
bool tagger_has_tag(tagger_t* tagger, const char* tag);

/// Resizes the tag with the given name if it exists in the tagger, preserving
/// existing data. If the given tag does not exist, this function does nothing.
/// \memberof tagger
void tagger_resize_tag(tagger_t* tagger, const char* tag, size_t new_num_indices);

/// Renames a tag within the tagger, preserving all information and associations.
/// \memberof tagger
void tagger_rename_tag(tagger_t* tagger, const char* old_tag, const char* new_tag);

/// Deletes a tag from the tagger, destroying all of its information (including properties).
/// \memberof tagger
void tagger_delete_tag(tagger_t* tagger, const char* tag);

/// Copies the indices of the given "source" tag to a "destination" tag. If
/// the "source" tag does not exist, this function has no effect. Otherwise,
/// the destination is overwritten. Any serializeable properties are copied.
/// \memberof tagger
void tagger_copy_tag(tagger_t* tagger, const char* source_tag, const char* destination_tag);

/// Replaces the given tag with its union with the other tag. If the other
/// tag does not exist, this function has no effect. Otherwise, the given
/// tag will be created if it does not yet exist.
/// \memberof tagger
void tagger_unite_tag(tagger_t* tagger, const char* tag, const char* other);

/// Replaces the given tag with its intersection with the other tag. If either
/// of the given tags does not exist in this tagger, this function has no effect.
/// \memberof tagger
void tagger_intersect_tag(tagger_t* tagger, const char* tag, const char* other);

/// Replaces the given tag with its difference with the other tag (i.e. tag
/// will contain no indices that appear in other after this operation). If
/// either of the given tags does not exist in this tagger, this function has no
/// effect.
/// \memberof tagger
void tagger_difference_tag(tagger_t* tagger, const char* tag, const char* other);

/// Allows iteration over the tags in the tagger, returning the name, indices,
/// and size of the tag. Returns true if a tag was found, or false if the
/// iteration has ended. Set pos to 0 to reset the iteration.
/// \memberof tagger
bool tagger_next_tag(tagger_t* tagger, int* pos, char** tag_name, int** tag_indices, size_t* tag_size);

/// Returns a serializer object that can read/write taggers from/to byte arrays.
/// \memberof tagger
serializer_t* tagger_serializer(void);

///@}

#endif
