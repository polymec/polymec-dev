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

// This type represents an object that manipulates byte streams, reading 
// and/or writing objects from them. Objects of this type are garbage-collected.
typedef struct serializer_t serializer_t;

// Constructs a new serializer with behavior defined by the given functions.
serializer_t* serializer_new(serializer_size_func size_func,
                             serializer_read_func read_func,
                             serializer_write_func write_func);

// Writes the given object to the byte stream at the given offset, changing 
// the offset in place.
void serializer_write(serializer_t* s, void* object, byte_array_t* byte_stream, size_t* offset);

// Returns an object read from the byte stream at the given offset, changing 
// the offset in place.
void* serializer_read(serializer_t* s, byte_array_t* byte_stream, size_t* offset);

// The following functions allow us to read/write primitives from/to byte arrays.

// Reads n characters from the byte array at the given offset, 
// placing it into the data array and updating the offset.
void byte_array_read_chars(byte_array_t* byte_stream, size_t n, size_t* offset, char* data);

// Writes n characters from the byte array at the given offset, updating the offset.
void byte_array_write_chars(byte_array_t* byte_stream, size_t n, char* data, size_t* offset);

// Reads n ints from the byte array at the given offset, 
// placing it into the data array and updating the offset.
void byte_array_read_ints(byte_array_t* byte_stream, size_t n, size_t* offset, int* data);

// Writes n ints from the byte array at the given offset, updating the offset.
void byte_array_write_ints(byte_array_t* byte_stream, size_t n, int* data, size_t* offset);

// Reads n longs from the byte array at the given offset, 
// placing it into the data array and updating the offset.
void byte_array_read_longs(byte_array_t* byte_stream, size_t n, size_t* offset, long* data);

// Writes n longs from the byte array at the given offset, updating the offset.
void byte_array_write_longs(byte_array_t* byte_stream, size_t n, long* data, size_t* offset);

// Reads n long longs from the byte array at the given offset, 
// placing it into the data array and updating the offset.
void byte_array_read_long_longs(byte_array_t* byte_stream, size_t n, size_t* offset, long long* data);

// Writes n long longs from the byte array at the given offset, updating the offset.
void byte_array_write_long_longs(byte_array_t* byte_stream, size_t n, long long* data, size_t* offset);

// Reads n reals from the byte array at the given offset, 
// placing it into the data array and updating the offset.
void byte_array_read_reals(byte_array_t* byte_stream, size_t n, size_t* offset, real_t* data);

// Writes n reals from the byte array at the given offset, updating the offset.
void byte_array_write_reals(byte_array_t* byte_stream, size_t n, real_t* data, size_t* offset);

// Reads n points from the byte array at the given offset, 
// placing it into the data array and updating the offset.
void byte_array_read_points(byte_array_t* byte_stream, size_t n, size_t* offset, point_t* data);

// Writes n points from the byte array at the given offset, updating the offset.
void byte_array_write_points(byte_array_t* byte_stream, size_t n, point_t* data, size_t* offset);

// Reads n vectors from the byte array at the given offset, 
// placing it into the data array and updating the offset.
void byte_array_read_vectors(byte_array_t* byte_stream, size_t n, size_t* offset, vector_t* data);

// Writes n vectors from the byte array at the given offset, updating the offset.
void byte_array_write_vectors(byte_array_t* byte_stream, size_t n, vector_t* data, size_t* offset);

#endif
