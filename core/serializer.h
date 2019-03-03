// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_SERIALIZER_H
#define POLYMEC_SERIALIZER_H

#include "core/point.h"
#include "core/point2.h"
#include "core/array.h"

/// \addtogroup core core
///@{

/// This prototype defines a function that measures the size of an object
/// to be serialized.
typedef size_t (*serializer_size_func)(void* context);

/// This prototype defines a function that reads an object from a byte array
/// at the given offset.
typedef void* (*serializer_read_func)(byte_array_t* stream, size_t* offset);

/// This prototype defines a function that writes an object to a byte array
/// at the given offset.
typedef void (*serializer_write_func)(void* context, byte_array_t* stream, size_t* offset);

/// This type represents a generic destructor associated with the cargo of
/// a serializer.
typedef void (*serializer_dtor_func)(void*);

/// \class serializer
/// This type represents an object that manipulates byte streams, reading
/// and/or writing objects from them.
typedef struct serializer_t serializer_t;

/// Constructs a new serializer with behavior defined by the given functions.
/// Every serializer has a unique name that identifies it within a global
/// registry of serializers.
/// \param [in] The name of the serializer (usually the type it serializes).
/// \param [in] size_func A function that returns the size of a given object
///                       to be serialized by the serializer.
/// \param [in] read_func A function that reads the contents of a byte_array
///                       and produces a newly-created object.
/// \param [in] write_func A function that writes the contents of an object to
///                        a byte_array.
/// \param [in] destructor_func A function that allows a serializer to destroy
///                        an object.
/// \memberof serializer
serializer_t* serializer_new(const char* name,
                             serializer_size_func size_func,
                             serializer_read_func read_func,
                             serializer_write_func write_func,
                             serializer_dtor_func destructor_func);

/// Returns the registered name of the serializer.
/// \memberof serializer
const char* serializer_name(serializer_t* s);

/// Returns the globally registered serializer with the given name, or
/// NULL if no such serializer has been registered.
/// \memberof serializer
serializer_t* serializer_from_name(const char* name);

/// Returns the size of the given object in bytes.
/// \memberof serializer
size_t serializer_size(serializer_t* s, void* object);

/// Writes the given object to the byte stream at the given offset, changing
/// the offset in place.
/// \param [in] object The object to be written to the byte stream.
/// \param [in,out] byte_stream The byte stream that stores the written object.
/// \param [in,out] offset Stores the offset in the byte stream to which the
///                        object is written.
/// \memberof serializer
void serializer_write(serializer_t* s, void* object, byte_array_t* byte_stream, size_t* offset);

/// Returns an object read from the byte stream at the given offset, changing
/// the offset in place.
/// \param [in,out] byte_stream The byte stream that stores the object to be read.
/// \param [in,out] offset Stores the offset in the byte stream from which the
///                        object is read.
/// \memberof serializer
void* serializer_read(serializer_t* s, byte_array_t* byte_stream, size_t* offset);

/// Creates and returns a copy of the given object via serialization.
/// \memberof serializer
void* serializer_clone_object(serializer_t* s, void* object);

/// Destroys the given object using the destructor registered with this
/// serializer at construction time.
/// \memberof serializer
void serializer_destroy_object(serializer_t* s, void* object);

//------------------------------------------------------------------------
// The following functions allow us to read/write primitives from/to
// byte arrays.
//------------------------------------------------------------------------

/// Reads n characters from the byte array at the given offset,
/// placing it into the data array and updating the offset.
/// \memberof byte_array
void byte_array_read_chars(byte_array_t* byte_stream, size_t n, char* data, size_t* offset);

/// Writes n characters to the byte array at the given offset, updating the offset.
/// \memberof byte_array
void byte_array_write_chars(byte_array_t* byte_stream, size_t n, char* data, size_t* offset);

/// Reads n ints from the byte array at the given offset,
/// placing it into the data array and updating the offset.
/// \memberof byte_array
void byte_array_read_ints(byte_array_t* byte_stream, size_t n, int* data, size_t* offset);

/// Writes n ints to the byte array at the given offset, updating the offset.
/// \memberof byte_array
void byte_array_write_ints(byte_array_t* byte_stream, size_t n, int* data, size_t* offset);

/// Reads n longs from the byte array at the given offset,
/// placing it into the data array and updating the offset.
/// \memberof byte_array
void byte_array_read_longs(byte_array_t* byte_stream, size_t n, long* data, size_t* offset);

/// Writes n longs to the byte array at the given offset, updating the offset.
/// \memberof byte_array
void byte_array_write_longs(byte_array_t* byte_stream, size_t n, long* data, size_t* offset);

/// Reads n long longs from the byte array at the given offset,
/// placing it into the data array and updating the offset.
/// \memberof byte_array
void byte_array_read_long_longs(byte_array_t* byte_stream, size_t n, long long* data, size_t* offset);

/// Writes n long longs to the byte array at the given offset, updating the offset.
/// \memberof byte_array
void byte_array_write_long_longs(byte_array_t* byte_stream, size_t n, long long* data, size_t* offset);

/// Reads n unsigned long longs from the byte array at the given offset,
/// placing it into the data array and updating the offset.
/// \memberof byte_array
void byte_array_read_unsigned_long_longs(byte_array_t* byte_stream, size_t n, unsigned long long* data, size_t* offset);

/// Writes n unsigned long longs to the byte array at the given offset, updating the offset.
/// \memberof byte_array
void byte_array_write_unsigned_long_longs(byte_array_t* byte_stream, size_t n, unsigned long long* data, size_t* offset);

/// Reads n index_ts from the byte array at the given offset,
/// placing it into the data array and updating the offset.
/// \memberof byte_array
void byte_array_read_index_ts(byte_array_t* byte_stream, size_t n, index_t* data, size_t* offset);

/// Writes n index_ts to the byte array at the given offset, updating the offset.
/// \memberof byte_array
void byte_array_write_index_ts(byte_array_t* byte_stream, size_t n, index_t* data, size_t* offset);

/// Reads n uint64_ts from the byte array at the given offset,
/// placing it into the data array and updating the offset.
/// \memberof byte_array
void byte_array_read_uint64_ts(byte_array_t* byte_stream, size_t n, uint64_t* data, size_t* offset);

/// Writes n uint64_ts to the byte array at the given offset, updating the offset.
/// \memberof byte_array
void byte_array_write_uint64_ts(byte_array_t* byte_stream, size_t n, uint64_t* data, size_t* offset);

/// Reads n size_ts from the byte array at the given offset,
/// placing it into the data array and updating the offset.
/// \memberof byte_array
void byte_array_read_size_ts(byte_array_t* byte_stream, size_t n, size_t* data, size_t* offset);

/// Writes n uint64s to the byte array at the given offset, updating the offset.
/// \memberof byte_array
void byte_array_write_size_ts(byte_array_t* byte_stream, size_t n, size_t* data, size_t* offset);

/// Reads n reals from the byte array at the given offset,
/// placing it into the data array and updating the offset.
/// \memberof byte_array
void byte_array_read_real_ts(byte_array_t* byte_stream, size_t n, real_t* data, size_t* offset);

/// Writes n reals to the byte array at the given offset, updating the offset.
/// \memberof byte_array
void byte_array_write_real_ts(byte_array_t* byte_stream, size_t n, real_t* data, size_t* offset);

#ifndef __cplusplus
/// Reads n complex numbers from the byte array at the given offset,
/// placing it into the data array and updating the offset.
/// \memberof byte_array
void byte_array_read_complex_ts(byte_array_t* byte_stream, size_t n, complex_t* data, size_t* offset);

/// Writes n complex numbers to the byte array at the given offset, updating the offset.
/// \memberof byte_array
void byte_array_write_complex_ts(byte_array_t* byte_stream, size_t n, complex_t* data, size_t* offset);
#endif

/// Reads n floats from the byte array at the given offset,
/// placing it into the data array and updating the offset.
/// \memberof byte_array
void byte_array_read_floats(byte_array_t* byte_stream, size_t n, float* data, size_t* offset);

/// Writes n floats to the byte array at the given offset, updating the offset.
/// \memberof byte_array
void byte_array_write_floats(byte_array_t* byte_stream, size_t n, float* data, size_t* offset);

/// Reads n doubles from the byte array at the given offset,
/// placing it into the data array and updating the offset.
/// \memberof byte_array
void byte_array_read_doubles(byte_array_t* byte_stream, size_t n, double* data, size_t* offset);

/// Writes n doubles to the byte array at the given offset, updating the offset.
/// \memberof byte_array
void byte_array_write_doubles(byte_array_t* byte_stream, size_t n, double* data, size_t* offset);

/// Reads n points from the byte array at the given offset,
/// placing it into the data array and updating the offset.
/// \memberof byte_array
void byte_array_read_points(byte_array_t* byte_stream, size_t n, point_t* data, size_t* offset);

/// Writes n points to the byte array at the given offset, updating the offset.
/// \memberof byte_array
void byte_array_write_points(byte_array_t* byte_stream, size_t n, point_t* data, size_t* offset);

/// Reads n 2D points from the byte array at the given offset,
/// placing it into the data array and updating the offset.
/// \memberof byte_array
void byte_array_read_point2s(byte_array_t* byte_stream, size_t n, point2_t* data, size_t* offset);

/// Writes n 2D points to the byte array at the given offset, updating the offset.
/// \memberof byte_array
void byte_array_write_point2s(byte_array_t* byte_stream, size_t n, point2_t* data, size_t* offset);

/// Reads n vectors from the byte array at the given offset,
/// placing it into the data array and updating the offset.
/// \memberof byte_array
void byte_array_read_vectors(byte_array_t* byte_stream, size_t n, vector_t* data, size_t* offset);

/// Writes n vectors to the byte array at the given offset, updating the offset.
/// \memberof byte_array
void byte_array_write_vectors(byte_array_t* byte_stream, size_t n, vector_t* data, size_t* offset);

///@{
/// Speaking of primitives, we bundle serializers for some primitive types.
serializer_t* string_serializer(void);
serializer_t* bbox_serializer(void);
serializer_t* int_array_serializer(void);
serializer_t* index_array_serializer(void);
serializer_t* real_array_serializer(void);
serializer_t* string_array_serializer(void);
// FIXME: More to come as needed.
///@}

///@}

#endif
