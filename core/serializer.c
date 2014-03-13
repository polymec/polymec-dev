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

#include <gc/gc.h>
#include "core/serializer.h"

struct serializer_t 
{
  serializer_size_func  size;
  serializer_read_func  read;
  serializer_write_func write;
};

serializer_t* serializer_new(serializer_size_func size_func,
                             serializer_read_func read_func,
                             serializer_write_func write_func)
{
  ASSERT(size_func != NULL);
  ASSERT(read_func != NULL);
  ASSERT(write_func != NULL);
  serializer_t* s = GC_MALLOC(sizeof(serializer_t));
  s->size = size_func;
  s->read = read_func;
  s->write = write_func;
  return s;
}

void serializer_write(serializer_t* s, void* object, byte_array_t* byte_stream, size_t* offset)
{
  size_t size = s->size(object);
  s->write(object, byte_stream, *offset);
  *offset += size;
}

void* serializer_read(serializer_t* s, byte_array_t* byte_stream, size_t* offset)
{
  void* object = s->read(byte_stream, *offset);
  size_t size = s->size(object);
  *offset += size;
  return object;
}

void byte_array_read_chars(byte_array_t* byte_stream, size_t n, size_t* offset, char* data)
{
  ASSERT(*offset < byte_stream->size + sizeof(char) * n);
  for (size_t i = 0; i < n; ++i)
  {
    char* datum = (char*)&byte_stream->data[*offset];
    data[i] = *datum;
    *offset += sizeof(char);
  }
}

void byte_array_read_ints(byte_array_t* byte_stream, size_t n, size_t* offset, int* data)
{
  ASSERT(*offset < byte_stream->size + sizeof(int) * n);
  for (size_t i = 0; i < n; ++i)
  {
    int* datum = (int*)&byte_stream->data[*offset];
    data[i] = *datum;
    *offset += sizeof(int);
  }
}

void byte_array_read_longs(byte_array_t* byte_stream, size_t n, size_t* offset, long* data)
{
  ASSERT(*offset < byte_stream->size + sizeof(long) * n);
  for (size_t i = 0; i < n; ++i)
  {
    long* datum = (long*)&byte_stream->data[*offset];
    data[i] = *datum;
    *offset += sizeof(long);
  }
}

void byte_array_read_long_longs(byte_array_t* byte_stream, size_t n, size_t* offset, long long* data)
{
  ASSERT(*offset < byte_stream->size + sizeof(long long) * n);
  for (size_t i = 0; i < n; ++i)
  {
    long long* datum = (long long*)&byte_stream->data[*offset];
    data[i] = *datum;
    *offset += sizeof(long long);
  }
}

void byte_array_read_reals(byte_array_t* byte_stream, size_t n, size_t* offset, real_t* data)
{
  ASSERT(*offset < byte_stream->size + sizeof(real_t) * n);
  for (size_t i = 0; i < n; ++i)
  {
    double* datum = (double*)&byte_stream->data[*offset];
    data[i] = *datum;
    *offset += sizeof(real_t);
  }
}

