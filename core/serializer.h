// Copyright (c) 2012-2015, Jeffrey N. Johnson
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

#ifndef POLYMEC_SERIALIZER_H
#define POLYMEC_SERIALIZER_H

#include "core/point.h"
#include "core/array.h"

// This prototype defines a function that measures the size of an object 
// to be serialized.
typedef size_t (*serializer_size_func)(void* context);

// This prototype defines a function that reads an object from a byte array
// at the given offset.
typedef void* (*serializer_read_func)(byte_array_t* stream, size_t* offset);

// This prototype defines a function that writes an object to a byte array 
// at the given offset.
typedef void (*serializer_write_func)(void* context, byte_array_t* stream, size_t* offset);

// This type represents a generic destructor associated with the cargo of 
// a serializer.
typedef void (*serializer_dtor_func)(void*);

// This type represents an object that manipulates byte streams, reading 
// and/or writing objects from them. Objects of this type are garbage-collected.
typedef struct serializer_t serializer_t;

// Constructs a new serializer with behavior defined by the given functions.
// Every serializer should have a unique name so that it can be added to a 
// global registry of serializers.
serializer_t* serializer_new(const char* name,
                             serializer_size_func size_func,
                             serializer_read_func read_func,
                             serializer_write_func write_func,
                             serializer_dtor_func destructor_func);

// Returns the registered name of the serializer.
const char* serializer_name(serializer_t* s);

// Returns the globally registered serializer with the given name, or 
// NULL if no such serializer has been registered.
serializer_t* serializer_from_name(const char* name);

// Returns the size of the given object in bytes.
size_t serializer_size(serializer_t* s, void* object);

// Writes the given object to the byte stream at the given offset, changing 
// the offset in place.
void serializer_write(serializer_t* s, void* object, byte_array_t* byte_stream, size_t* offset);

// Returns an object read from the byte stream at the given offset, changing 
// the offset in place.
void* serializer_read(serializer_t* s, byte_array_t* byte_stream, size_t* offset);

// Destroys the given object using the destructor registered with this 
// serializer at construction time.
void serializer_destroy(serializer_t* s, void* object);

//------------------------------------------------------------------------
// The following functions allow us to read/write primitives from/to 
// byte arrays.
//------------------------------------------------------------------------

// Reads n characters from the byte array at the given offset, 
// placing it into the data array and updating the offset.
void byte_array_read_chars(byte_array_t* byte_stream, size_t n, char* data, size_t* offset);

// Writes n characters to the byte array at the given offset, updating the offset.
void byte_array_write_chars(byte_array_t* byte_stream, size_t n, char* data, size_t* offset);

// Reads n ints from the byte array at the given offset, 
// placing it into the data array and updating the offset.
void byte_array_read_ints(byte_array_t* byte_stream, size_t n, int* data, size_t* offset);

// Writes n ints to the byte array at the given offset, updating the offset.
void byte_array_write_ints(byte_array_t* byte_stream, size_t n, int* data, size_t* offset);

// Reads n longs from the byte array at the given offset, 
// placing it into the data array and updating the offset.
void byte_array_read_longs(byte_array_t* byte_stream, size_t n, long* data, size_t* offset);

// Writes n longs to the byte array at the given offset, updating the offset.
void byte_array_write_longs(byte_array_t* byte_stream, size_t n, long* data, size_t* offset);

// Reads n long longs from the byte array at the given offset, 
// placing it into the data array and updating the offset.
void byte_array_read_long_longs(byte_array_t* byte_stream, size_t n, long long* data, size_t* offset);

// Writes n long longs to the byte array at the given offset, updating the offset.
void byte_array_write_long_longs(byte_array_t* byte_stream, size_t n, long long* data, size_t* offset);

// Reads n uint64s from the byte array at the given offset, 
// placing it into the data array and updating the offset.
void byte_array_read_uint64s(byte_array_t* byte_stream, size_t n, uint64_t* data, size_t* offset);

// Writes n uint64s to the byte array at the given offset, updating the offset.
void byte_array_write_uint64s(byte_array_t* byte_stream, size_t n, uint64_t* data, size_t* offset);

// Reads n reals from the byte array at the given offset, 
// placing it into the data array and updating the offset.
void byte_array_read_reals(byte_array_t* byte_stream, size_t n, real_t* data, size_t* offset);

// Writes n reals to the byte array at the given offset, updating the offset.
void byte_array_write_reals(byte_array_t* byte_stream, size_t n, real_t* data, size_t* offset);

// Reads n points from the byte array at the given offset, 
// placing it into the data array and updating the offset.
void byte_array_read_points(byte_array_t* byte_stream, size_t n, point_t* data, size_t* offset);

// Writes n points to the byte array at the given offset, updating the offset.
void byte_array_write_points(byte_array_t* byte_stream, size_t n, point_t* data, size_t* offset);

// Reads n vectors from the byte array at the given offset, 
// placing it into the data array and updating the offset.
void byte_array_read_vectors(byte_array_t* byte_stream, size_t n, vector_t* data, size_t* offset);

// Writes n vectors to the byte array at the given offset, updating the offset.
void byte_array_write_vectors(byte_array_t* byte_stream, size_t n, vector_t* data, size_t* offset);

// Speaking of primitives, we bundle serializers for some primitive types.
serializer_t* string_serializer();
serializer_t* bbox_serializer();

// FIXME: More to come as needed.

#endif
