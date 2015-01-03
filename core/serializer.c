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

#include <gc/gc.h>
#include "core/serializer.h"
#include "core/unordered_map.h"

// A serializer is just a set of read/write functions.
struct serializer_t 
{
  char* name;
  serializer_size_func  size;
  serializer_read_func  read;
  serializer_write_func write;
  serializer_dtor_func  dtor;
};

// We keep a global registry of serializers so that objects can be sent 
// across the network portably.
DEFINE_UNORDERED_MAP(serializer_registry, char*, serializer_t*, string_hash, string_equals)
static serializer_registry_t* registry = NULL;

static void destroy_serializer_registry()
{
  ASSERT(registry != NULL);
  serializer_registry_free(registry);
}

static void serializer_free(void* ctx, void* dummy)
{
  serializer_t* s = ctx;
  string_free(s->name);
}

serializer_t* serializer_new(const char* name,
                             serializer_size_func size_func,
                             serializer_read_func read_func,
                             serializer_write_func write_func,
                             serializer_dtor_func destructor_func)
{
  ASSERT(size_func != NULL);
  ASSERT(read_func != NULL);
  ASSERT(write_func != NULL);

  // If we haven't created our registry yet, do so now.
  if (registry == NULL)
  {
    registry = serializer_registry_new();
    polymec_atexit(destroy_serializer_registry);
  }

  // If our global registry does not contain this serializer, make one 
  // and insert it.
  serializer_t* s = NULL;
  serializer_t** s_ptr = serializer_registry_get(registry, (char*)name);
  if (s_ptr == NULL)
  {
    s = GC_MALLOC(sizeof(serializer_t));
    s->name = string_dup(name);
    s->size = size_func;
    s->read = read_func;
    s->write = write_func;
    s->dtor = destructor_func;
    GC_register_finalizer(s, serializer_free, s, NULL, NULL);
    serializer_registry_insert_with_k_dtor(registry, string_dup(name), s, string_free);
  }
  else
  {
    // Just fetch the entry in the table.
    s = *s_ptr;
  }

  return s;
}

serializer_t* serializer_from_name(const char* name)
{
  if (registry == NULL)
    return NULL;
  serializer_t** s_ptr = serializer_registry_get(registry, (char*)name);
  if (s_ptr != NULL)
    return *s_ptr;
  else
    return NULL;
}

const char* serializer_name(serializer_t* s)
{
  return (const char*)s->name;
}

size_t serializer_size(serializer_t* s, void* object)
{
  return s->size(object);
}

void serializer_write(serializer_t* s, void* object, byte_array_t* byte_stream, size_t* offset)
{
  s->write(object, byte_stream, offset);
  ASSERT(byte_stream->size >= s->size(object));
}

void* serializer_read(serializer_t* s, byte_array_t* byte_stream, size_t* offset)
{
  void* object = s->read(byte_stream, offset);
  return object;
}

void serializer_destroy(serializer_t* s, void* object)
{
  if ((s->dtor != NULL) && (object != NULL))
    s->dtor(object);
}

void byte_array_read_chars(byte_array_t* byte_stream, size_t n, char* data, size_t* offset)
{
  ASSERT(*offset <= byte_stream->size - sizeof(char) * n);
  memcpy(data, &byte_stream->data[*offset], sizeof(char) * n);
  *offset += n * sizeof(char);
}

void byte_array_write_chars(byte_array_t* byte_stream, size_t n, char* data, size_t* offset)
{
  byte_array_resize(byte_stream, *offset + n * sizeof(char));
  memcpy(&byte_stream->data[*offset], data, sizeof(char) * n);
  *offset += sizeof(char) * n;
}

void byte_array_read_ints(byte_array_t* byte_stream, size_t n, int* data, size_t* offset)
{
  ASSERT(*offset <= byte_stream->size - sizeof(int) * n);
  memcpy(data, &byte_stream->data[*offset], sizeof(int) * n);
  *offset += n * sizeof(int);
}

void byte_array_write_ints(byte_array_t* byte_stream, size_t n, int* data, size_t* offset)
{
  byte_array_resize(byte_stream, *offset + n * sizeof(int));
  memcpy(&byte_stream->data[*offset], data, sizeof(int) * n);
  *offset += sizeof(int) * n;
}

void byte_array_read_longs(byte_array_t* byte_stream, size_t n, long* data, size_t* offset)
{
  ASSERT(*offset <= byte_stream->size - sizeof(long) * n);
  memcpy(data, &byte_stream->data[*offset], sizeof(long) * n);
  *offset += n * sizeof(long);
}

void byte_array_write_longs(byte_array_t* byte_stream, size_t n, long* data, size_t* offset)
{
  byte_array_resize(byte_stream, *offset + n * sizeof(long));
  memcpy(&byte_stream->data[*offset], data, sizeof(long) * n);
  *offset += sizeof(long) * n;
}

void byte_array_read_long_longs(byte_array_t* byte_stream, size_t n, long long* data, size_t* offset)
{
  ASSERT(*offset <= byte_stream->size - sizeof(long long) * n);
  memcpy(data, &byte_stream->data[*offset], sizeof(long long) * n);
  *offset += n * sizeof(long long);
}

void byte_array_write_long_longs(byte_array_t* byte_stream, size_t n, long long* data, size_t* offset)
{
  byte_array_resize(byte_stream, *offset + n * sizeof(long long));
  memcpy(&byte_stream->data[*offset], data, sizeof(long long) * n);
  *offset += sizeof(long long) * n;
}

void byte_array_read_uint64s(byte_array_t* byte_stream, size_t n, uint64_t* data, size_t* offset)
{
  ASSERT(*offset <= byte_stream->size - sizeof(uint64_t) * n);
  memcpy(data, &byte_stream->data[*offset], sizeof(uint64_t) * n);
  *offset += n * sizeof(uint64_t);
}

void byte_array_write_uint64s(byte_array_t* byte_stream, size_t n, uint64_t* data, size_t* offset)
{
  byte_array_resize(byte_stream, *offset + n * sizeof(uint64_t));
  memcpy(&byte_stream->data[*offset], data, sizeof(uint64_t) * n);
  *offset += sizeof(uint64_t) * n;
}

void byte_array_read_reals(byte_array_t* byte_stream, size_t n, real_t* data, size_t* offset)
{
  ASSERT(*offset <= byte_stream->size - sizeof(real_t) * n);
  memcpy(data, &byte_stream->data[*offset], sizeof(real_t) * n);
  *offset += n * sizeof(real_t);
}

void byte_array_write_reals(byte_array_t* byte_stream, size_t n, real_t* data, size_t* offset)
{
  byte_array_resize(byte_stream, *offset + n * sizeof(real_t));
  memcpy(&byte_stream->data[*offset], data, sizeof(real_t) * n);
  *offset += sizeof(real_t) * n;
}


void byte_array_read_points(byte_array_t* byte_stream, size_t n, point_t* data, size_t* offset)
{
  ASSERT(*offset <= byte_stream->size - sizeof(point_t) * n);
  memcpy(data, &byte_stream->data[*offset], sizeof(point_t) * n);
  *offset += n * sizeof(point_t);
}

void byte_array_write_points(byte_array_t* byte_stream, size_t n, point_t* data, size_t* offset)
{
  byte_array_resize(byte_stream, *offset + n * sizeof(point_t));
  memcpy(&byte_stream->data[*offset], data, sizeof(point_t) * n);
  *offset += sizeof(point_t) * n;
}

void byte_array_read_vectors(byte_array_t* byte_stream, size_t n, vector_t* data, size_t* offset)
{
  ASSERT(*offset <= byte_stream->size - sizeof(vector_t) * n);
  memcpy(data, &byte_stream->data[*offset], sizeof(vector_t) * n);
  *offset += n * sizeof(vector_t);
}

void byte_array_write_vectors(byte_array_t* byte_stream, size_t n, vector_t* data, size_t* offset)
{
  byte_array_resize(byte_stream, *offset + n * sizeof(vector_t));
  memcpy(&byte_stream->data[*offset], data, sizeof(vector_t) * n);
  *offset += sizeof(vector_t) * n;
}

static size_t string_byte_size(void* obj)
{
  char* s = obj;
  return (size_t)(sizeof(int) + strlen(s)*sizeof(char));
}

static void* string_byte_read(byte_array_t* bytes, size_t* offset)
{
  int len;
  byte_array_read_ints(bytes, 1, &len, offset);
  char* s = polymec_malloc(sizeof(char) * (len+1));
  byte_array_read_chars(bytes, len, s, offset);
  s[len] = '\0';
  return s;
}

static void string_byte_write(void* obj, byte_array_t* bytes, size_t* offset)
{
  char* s = obj;
  int len = strlen(s);
  byte_array_write_ints(bytes, 1, &len, offset);
  byte_array_write_chars(bytes, len, s, offset);
}

serializer_t* string_serializer()
{
  return serializer_new("string", string_byte_size, string_byte_read, string_byte_write, DTOR(string_free));
}

static size_t bb_byte_size(void* obj)
{
  return 6 * sizeof(real_t);
}

static void* bb_byte_read(byte_array_t* bytes, size_t* offset)
{
  real_t data[6];
  byte_array_read_reals(bytes, 6, data, offset);
  return bbox_new(data[0], data[1], data[2], data[3], data[4], data[5]);
}

static void bb_byte_write(void* obj, byte_array_t* bytes, size_t* offset)
{
  real_t* b = obj;
  byte_array_write_reals(bytes, 6, b, offset);
}

serializer_t* bbox_serializer()
{
  return serializer_new("bbox", bb_byte_size,
                        bb_byte_read, bb_byte_write, NULL);
}

