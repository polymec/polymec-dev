// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

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

static void serializer_free(serializer_t* s)
{
  string_free(s->name);
  polymec_free(s);
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
    s = polymec_malloc(sizeof(serializer_t));
    s->name = string_dup(name);
    s->size = size_func;
    s->read = read_func;
    s->write = write_func;
    s->dtor = destructor_func;
    serializer_registry_insert_with_kv_dtors(registry, string_dup(name), s, string_free, serializer_free);
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

void* serializer_clone_object(serializer_t* s, void* object)
{
  byte_array_t* bytes = byte_array_new();
  size_t offset = 0;
  serializer_write(s, object, bytes, &offset);
  ASSERT(offset == serializer_size(s, object));

  offset = 0;
  void* clone = serializer_read(s, bytes, &offset);
  byte_array_free(bytes);
  return clone;
}

void serializer_destroy_object(serializer_t* s, void* object)
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

void byte_array_read_unsigned_long_longs(byte_array_t* byte_stream, size_t n, unsigned long long* data, size_t* offset)
{
  ASSERT(*offset <= byte_stream->size - sizeof(unsigned long long) * n);
  memcpy(data, &byte_stream->data[*offset], sizeof(unsigned long long) * n);
  *offset += n * sizeof(unsigned long long);
}

void byte_array_write_unsigned_long_longs(byte_array_t* byte_stream, size_t n, unsigned long long* data, size_t* offset)
{
  byte_array_resize(byte_stream, *offset + n * sizeof(unsigned long long));
  memcpy(&byte_stream->data[*offset], data, sizeof(unsigned long long) * n);
  *offset += sizeof(unsigned long long) * n;
}

void byte_array_read_index_ts(byte_array_t* byte_stream, size_t n, index_t* data, size_t* offset)
{
  ASSERT(*offset <= byte_stream->size - sizeof(index_t) * n);
  memcpy(data, &byte_stream->data[*offset], sizeof(index_t) * n);
  *offset += n * sizeof(index_t);
}

void byte_array_write_index_ts(byte_array_t* byte_stream, size_t n, index_t* data, size_t* offset)
{
  byte_array_resize(byte_stream, *offset + n * sizeof(index_t));
  memcpy(&byte_stream->data[*offset], data, sizeof(index_t) * n);
  *offset += sizeof(index_t) * n;
}

void byte_array_read_uint64_ts(byte_array_t* byte_stream, size_t n, uint64_t* data, size_t* offset)
{
  ASSERT(*offset <= byte_stream->size - sizeof(uint64_t) * n);
  memcpy(data, &byte_stream->data[*offset], sizeof(uint64_t) * n);
  *offset += n * sizeof(uint64_t);
}

void byte_array_write_uint64_ts(byte_array_t* byte_stream, size_t n, uint64_t* data, size_t* offset)
{
  byte_array_resize(byte_stream, *offset + n * sizeof(uint64_t));
  memcpy(&byte_stream->data[*offset], data, sizeof(uint64_t) * n);
  *offset += sizeof(uint64_t) * n;
}

void byte_array_read_size_ts(byte_array_t* byte_stream, size_t n, size_t* data, size_t* offset)
{
  ASSERT(*offset <= byte_stream->size - sizeof(size_t) * n);
  memcpy(data, &byte_stream->data[*offset], sizeof(size_t) * n);
  *offset += n * sizeof(size_t);
}

void byte_array_write_size_ts(byte_array_t* byte_stream, size_t n, size_t* data, size_t* offset)
{
  byte_array_resize(byte_stream, *offset + n * sizeof(size_t));
  memcpy(&byte_stream->data[*offset], data, sizeof(size_t) * n);
  *offset += sizeof(size_t) * n;
}

void byte_array_read_real_ts(byte_array_t* byte_stream, size_t n, real_t* data, size_t* offset)
{
  ASSERT(*offset <= byte_stream->size - sizeof(real_t) * n);
  memcpy(data, &byte_stream->data[*offset], sizeof(real_t) * n);
  *offset += n * sizeof(real_t);
}

void byte_array_write_real_ts(byte_array_t* byte_stream, size_t n, real_t* data, size_t* offset)
{
  byte_array_resize(byte_stream, *offset + n * sizeof(real_t));
  memcpy(&byte_stream->data[*offset], data, sizeof(real_t) * n);
  *offset += sizeof(real_t) * n;
}

void byte_array_read_complex_ts(byte_array_t* byte_stream, size_t n, complex_t* data, size_t* offset)
{
  ASSERT(*offset <= byte_stream->size - sizeof(complex_t) * n);
  memcpy(data, &byte_stream->data[*offset], sizeof(complex_t) * n);
  *offset += n * sizeof(complex_t);
}

void byte_array_write_complex_ts(byte_array_t* byte_stream, size_t n, complex_t* data, size_t* offset)
{
  byte_array_resize(byte_stream, *offset + n * sizeof(complex_t));
  memcpy(&byte_stream->data[*offset], data, sizeof(complex_t) * n);
  *offset += sizeof(complex_t) * n;
}

void byte_array_read_floats(byte_array_t* byte_stream, size_t n, float* data, size_t* offset)
{
  ASSERT(*offset <= byte_stream->size - sizeof(float) * n);
  memcpy(data, &byte_stream->data[*offset], sizeof(float) * n);
  *offset += n * sizeof(float);
}

void byte_array_write_floats(byte_array_t* byte_stream, size_t n, float* data, size_t* offset)
{
  byte_array_resize(byte_stream, *offset + n * sizeof(float));
  memcpy(&byte_stream->data[*offset], data, sizeof(float) * n);
  *offset += sizeof(float) * n;
}

void byte_array_read_doubles(byte_array_t* byte_stream, size_t n, double* data, size_t* offset)
{
  ASSERT(*offset <= byte_stream->size - sizeof(double) * n);
  memcpy(data, &byte_stream->data[*offset], sizeof(double) * n);
  *offset += n * sizeof(double);
}

void byte_array_write_doubles(byte_array_t* byte_stream, size_t n, double* data, size_t* offset)
{
  byte_array_resize(byte_stream, *offset + n * sizeof(double));
  memcpy(&byte_stream->data[*offset], data, sizeof(double) * n);
  *offset += sizeof(double) * n;
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

void byte_array_read_point2s(byte_array_t* byte_stream, size_t n, point2_t* data, size_t* offset)
{
  ASSERT(*offset <= byte_stream->size - sizeof(point2_t) * n);
  memcpy(data, &byte_stream->data[*offset], sizeof(point2_t) * n);
  *offset += n * sizeof(point2_t);
}

void byte_array_write_point2s(byte_array_t* byte_stream, size_t n, point2_t* data, size_t* offset)
{
  byte_array_resize(byte_stream, *offset + n * sizeof(point2_t));
  memcpy(&byte_stream->data[*offset], data, sizeof(point2_t) * n);
  *offset += sizeof(point2_t) * n;
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
  int len = (int)strlen(s);
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
  byte_array_read_real_ts(bytes, 6, data, offset);
  return bbox_new(data[0], data[1], data[2], data[3], data[4], data[5]);
}

static void bb_byte_write(void* obj, byte_array_t* bytes, size_t* offset)
{
  real_t* b = obj;
  byte_array_write_real_ts(bytes, 6, b, offset);
}

serializer_t* bbox_serializer()
{
  return serializer_new("bbox", bb_byte_size,
                        bb_byte_read, bb_byte_write, NULL);
}

#define DEFINE_ARRAY_SERIALIZER(array_name, element) \
static size_t array_name##_byte_size(void* obj) \
{ \
  array_name##_t* a = obj; \
  return sizeof(size_t) + sizeof(size_t) + sizeof(element) * a->size; \
} \
\
static void* array_name##_byte_read(byte_array_t* bytes, size_t* offset) \
{ \
  size_t size, capacity; \
  byte_array_read_size_ts(bytes, 1, &size, offset); \
  byte_array_read_size_ts(bytes, 1, &capacity, offset); \
  ASSERT(capacity >= size); \
  array_name##_t* a = array_name##_new_with_capacity(capacity); \
  byte_array_read_##element##s(bytes, size, a->data, offset); \
  return a; \
} \
\
static void array_name##_byte_write(void* obj, byte_array_t* bytes, size_t* offset) \
{ \
  array_name##_t* a = obj; \
  byte_array_write_size_ts(bytes, 1, &a->size, offset); \
  byte_array_write_size_ts(bytes, 1, &a->capacity, offset); \
  byte_array_write_##element##s(bytes, a->size, a->data, offset); \
} \
\
serializer_t* array_name##_serializer() \
{ \
  return serializer_new(#array_name, array_name##_byte_size, \
                        array_name##_byte_read, array_name##_byte_write, \
                        DTOR(array_name##_free)); \
} \

DEFINE_ARRAY_SERIALIZER(int_array, int)
DEFINE_ARRAY_SERIALIZER(index_array, index_t)
DEFINE_ARRAY_SERIALIZER(real_array, real_t)

// String array stuff -- a one-off.
static size_t string_array_byte_size(void* obj)
{
  string_array_t* a = obj;
  return sizeof(size_t) + sizeof(size_t) + sizeof(char*) * a->size;
}

static void* string_array_byte_read(byte_array_t* bytes, size_t* offset)
{
  size_t size, capacity;
  byte_array_read_size_ts(bytes, 1, &size, offset);
  byte_array_read_size_ts(bytes, 1, &capacity, offset);
  ASSERT(capacity >= size);
  string_array_t* a = string_array_new_with_capacity(capacity);
  for (int i = 0; i < size; ++i)
    string_array_append_with_dtor(a, string_byte_read(bytes, offset), string_free);
  return a;
}

static void string_array_byte_write(void* obj, byte_array_t* bytes, size_t* offset)
{
  string_array_t* a = obj;
  byte_array_write_size_ts(bytes, 1, &a->size, offset);
  byte_array_write_size_ts(bytes, 1, &a->capacity, offset);
  for (int i = 0; i < a->size; ++i)
    string_byte_write(a->data[i], bytes, offset);
}

serializer_t* string_array_serializer()
{
  return serializer_new("string_array", string_array_byte_size,
                        string_array_byte_read, string_array_byte_write,
                        DTOR(string_array_free));
}

