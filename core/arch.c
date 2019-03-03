// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <sys/mman.h>
#include <errno.h>
#include "core/arch.h"
#include "core/polymec.h"

#ifdef APPLE

typedef struct
{
  size_t pos;
  size_t size;
  char* buffer;
  bool own_buff;
} fmem_t;

static int fmem_read(void* context, char* buf, int size)
{
  fmem_t* mem = context;
  size_t available = mem->size - mem->pos;

  if (size > available)
    size = (int)available;

  memcpy(buf, mem->buffer + mem->pos, sizeof(char) * size);
  mem->pos += size;

  return size;
}

static int fmem_write(void* context, const char* buf, int size)
{
  fmem_t* mem = context;
  size_t available = mem->size - mem->pos;

  if (size > available)
    return -1;

  memcpy(mem->buffer + mem->pos, buf, sizeof(char)*size);
  mem->pos += size;

  return size;
}

static fpos_t fmem_seek(void* context, fpos_t offset, int whence)
{
  size_t pos;
  fmem_t* mem = context;

  switch (whence)
  {
    case SEEK_SET: pos = offset; break;
    case SEEK_CUR: pos = mem->pos + offset; break;
    case SEEK_END: pos = mem->size + offset; break;
    default: errno = EINVAL; return -1;
  }

  if (pos > mem->size)
    return -1;

  mem->pos = pos;
  return (fpos_t)pos;
}

static int fmem_close(void *context)
{
  fmem_t* fmem = context;
  if (fmem->own_buff)
    polymec_free(fmem->buffer);
  polymec_free(fmem);
  return 0;
}

FILE *fmemopen(void *buf, size_t size, const char *mode)
{
  fmem_t* fmem = polymec_malloc(sizeof(fmem_t)); // Context pointer.
  fmem->pos = 0;
  fmem->size = size;
  if (buf != NULL)
  {
    fmem->buffer = buf;
    fmem->own_buff = false;
  }
  else
  {
    fmem->buffer = polymec_malloc(size);
    fmem->own_buff = true;
  }
  return funopen(fmem, fmem_read, fmem_write, fmem_seek, fmem_close);
}

typedef struct
{
  char **buf;
  size_t *len;
  size_t pos;
  size_t eof;
  size_t allocated;
  char c;
} memstream_t;

static int memstream_write(void* context, const char *buf, int n)
{
  memstream_t* stream = context;
  char *cbuf = *stream->buf;

  // Prevent overflow.
  if ((ssize_t) (stream->pos + n) < 0)
  {
    errno = EFBIG;
    return EOF;
  }

  if (stream->allocated <= stream->pos + n)
  {
    // Grow the buffer if needed.
    size_t newsize = stream->allocated * 3 / 2;
    if (newsize < stream->pos + n + 1)
      newsize = stream->pos + n + 1;
    cbuf = realloc(cbuf, newsize);
    if (!cbuf)
      return EOF;
    *stream->buf = cbuf;
    stream->allocated = newsize;
  }

  // Ensure that everything beyond EOF is null.
  if (stream->eof < stream->pos)
    memset(cbuf + stream->eof, '\0', stream->pos - stream->eof);

  memcpy(cbuf + stream->pos, buf, n);
  stream->pos += n;
  if (stream->eof < stream->pos)
    stream->eof = stream->pos;
  else
    stream->c = cbuf[stream->pos];
  cbuf[stream->pos] = '\0';
  *stream->len = stream->pos;
  return n;
}

static fpos_t memstream_seek(void *context, fpos_t pos, int whence)
{
  memstream_t *stream = context;
  off_t offset = pos;

  if (whence == SEEK_CUR)
    offset += stream->pos;
  else if (whence == SEEK_END)
    offset += stream->eof;
  if (offset < 0)
  {
    errno = EINVAL;
    offset = -1;
  }
  else if ((size_t)offset != offset)
  {
    errno = ENOSPC;
    offset = -1;
  }
  else
  {
    if (stream->pos < stream->eof)
    {
      (*stream->buf)[stream->pos] = stream->c;
      stream->c = '\0';
    }
    stream->pos = offset;
    if (stream->pos < stream->eof)
    {
      stream->c = (*stream->buf)[stream->pos];
      (*stream->buf)[stream->pos] = '\0';
      *stream->len = stream->pos;
    }
    else
      *stream->len = stream->eof;
  }
  return offset;
}

static int memstream_close(void *c)
{
  memstream_t *stream = c;
  char *buf;

  buf = realloc(*stream->buf, *stream->len + 1);
  if (buf != NULL)
    *stream->buf = buf;
  polymec_free(stream);
  return 0;
}

FILE* open_memstream(char **buf, size_t *len)
{
  FILE *f;
  memstream_t *stream;

  if (!buf || !len)
  {
    errno = EINVAL;
    return NULL;
  }
  stream = polymec_malloc(sizeof(memstream_t));
  *buf = malloc(32 * sizeof(char));
  **buf = '\0';
  *len = 0;

  f = funopen(stream, NULL, memstream_write, memstream_seek, memstream_close);
  if (f == NULL)
  {
    int saved_errno = errno;
    polymec_free(stream);
    errno = saved_errno;
  }
  else
  {
    stream->buf = buf;
    stream->len = len;
    stream->pos = 0;
    stream->eof = 0;
    stream->c = '\0';
    stream->allocated = 32;
  }
  return f;
}

#endif

