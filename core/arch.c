//
// Copyright 2012 Jeff Verkoeyen
// Originally ported from https://github.com/ingenuitas/python-tesseract/blob/master/fmemopen.c
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>
#include <assert.h>
#include <errno.h>

struct fmem {
  size_t pos;
  size_t size;
  char *buffer;
};
typedef struct fmem fmem_t;

static int readfn(void *handler, char *buf, int size) {
  fmem_t *mem = handler;
  size_t available = mem->size - mem->pos;
  
  if (size > available) {
    size = available;
  }
  memcpy(buf, mem->buffer + mem->pos, sizeof(char) * size);
  mem->pos += size;
  
  return size;
}

static int writefn(void *handler, const char *buf, int size) {
  fmem_t *mem = handler;
  size_t available = mem->size - mem->pos;

  if (size > available) {
    size = available;
  }
  memcpy(mem->buffer + mem->pos, buf, sizeof(char) * size);
  mem->pos += size;

  return size;
}

static fpos_t seekfn(void *handler, fpos_t offset, int whence) {
  size_t pos;
  fmem_t *mem = handler;

  switch (whence) {
    case SEEK_SET: pos = offset; break;
    case SEEK_CUR: pos = mem->pos + offset; break;
    case SEEK_END: pos = mem->size + offset; break;
    default: return -1;
  }

  if (pos > mem->size) {
    return -1;
  }

  mem->pos = pos;
  return (fpos_t)pos;
}

static int closefn(void *handler) {
  free(handler);
  return 0;
}

FILE *fmemopen(void *buf, size_t size, const char *mode) {
  // This data is released on fclose.
  fmem_t* mem = (fmem_t *) malloc(sizeof(fmem_t));

  // Zero-out the structure.
  memset(mem, 0, sizeof(fmem_t));

  mem->size = size;
  mem->buffer = buf;

  // funopen's man page: https://developer.apple.com/library/mac/#documentation/Darwin/Reference/ManPages/man3/funopen.3.html
  return funopen(mem, readfn, writefn, seekfn, closefn);
}

/* Open a write stream around a malloc'd string.
   Copyright (C) 2010 Free Software Foundation, Inc.

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.  */

/* Written by Eric Blake <address@hidden>, 2010.  */

/* Specification.  */

# define INITIAL_ALLOC 64

struct data
{
  char **buf; /* User's argument.  */
  size_t *len; /* User's argument.  Smaller of pos or eof.  */
  size_t pos; /* Current position.  */
  size_t eof; /* End-of-file position.  */
  size_t allocated; /* Allocated size of *buf, always > eof.  */
  char c; /* Temporary storage for byte overwritten by NUL, if pos < eof.  */
};
typedef struct data data;

/* Stupid BSD interface uses int/int instead of ssize_t/size_t.  */
//verify (sizeof (int) <= sizeof (size_t));
//verify (sizeof (int) <= sizeof (ssize_t));

static int mem_write(void *c, const char *buf, int n)
{
  data *context = c;
  char *cbuf = *context->buf;

  /* Be sure we don't overflow.  */
  if ((ssize_t) (context->pos + n) < 0)
  {
    errno = EFBIG;
    return EOF;
  }
  /* Grow the buffer, if necessary.  Use geometric growth to avoid
     quadratic realloc behavior.  Overallocate, to accomodate the
     requirement to always place a trailing NUL not counted by length.
     Thus, we want max(prev_size*1.5, context->pos+n+1).  */
  if (context->allocated <= context->pos + n)
  {
    size_t newsize = context->allocated * 3 / 2;
    if (newsize < context->pos + n + 1)
      newsize = context->pos + n + 1;
    cbuf = realloc(cbuf, newsize);
    if (!cbuf)
      return EOF;
    *context->buf = cbuf;
    context->allocated = newsize;
  }
  /* If we have previously done a seek beyond eof, ensure all
     intermediate bytges are NUL.  */
  if (context->eof < context->pos)
    memset (cbuf + context->eof, '\0', context->pos - context->eof);
  memcpy (cbuf + context->pos, buf, n);
  context->pos += n;
  /* If the user has previously written beyond the current position,
     remember what the trailing NUL is overwriting.  Otherwise,
     extend the stream.  */
  if (context->eof < context->pos)
    context->eof = context->pos;
  else
    context->c = cbuf[context->pos];
  cbuf[context->pos] = '\0';
  *context->len = context->pos;
  return n;
}

static fpos_t mem_seek(void *c, fpos_t pos, int whence) 
{
  data *context = c;
  off_t offset = pos;

  if (whence == SEEK_CUR)
    offset += context->pos;
  else if (whence == SEEK_END)
    offset += context->eof;
  if (offset < 0)
  {
    errno = EINVAL;
    offset = -1;
  }
  else if ((size_t) offset != offset)
  {
    errno = ENOSPC;
    offset = -1;
  }
  else
  {
    if (context->pos < context->eof)
    {
      (*context->buf)[context->pos] = context->c;
      context->c = '\0';
    }
    context->pos = offset;
    if (context->pos < context->eof)
    {
      context->c = (*context->buf)[context->pos];
      (*context->buf)[context->pos] = '\0';
      *context->len = context->pos;
    }
    else
      *context->len = context->eof;
  }
  return offset;
}

static int mem_close(void *c)
{
  data *context = c;
  char *buf;

  /* Be nice and try to reduce excess memory.  */
  buf = realloc(*context->buf, *context->len + 1);
  if (buf)
    *context->buf = buf;
  free (context);
  return 0;
}

FILE* open_memstream(char **buf, size_t *len)
{
  FILE *f;
  data *context;

  if (!buf || !len)
  {
    errno = EINVAL;
    return NULL;
  }
  context = malloc(sizeof(data));
  *buf = malloc(INITIAL_ALLOC);
  **buf = '\0';
  *len = 0;

  f = funopen(context, NULL, mem_write, mem_seek, mem_close);
  if (!f)
  {
    int saved_errno = errno;
    free (context);
    errno = saved_errno;
  }
  else
  {
    context->buf = buf;
    context->len = len;
    context->pos = 0;
    context->eof = 0;
    context->c = '\0';
    context->allocated = INITIAL_ALLOC;
  }
  return f;
}

