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

#ifndef POLYMEC_LOGGERS_H
#define POLYMEC_LOGGERS_H

#include <stdio.h>

// Types of log messages / log levels. The levels are increasing in amount 
// of output: LOG_WARNING logs only warnings, LOG_INFO logs warnings and 
// informational messages, and so on and so forth.
typedef enum
{
  LOG_NONE,                   // Log no messages
  LOG_URGENT,                 // Log only urgent messages
  LOG_INFO,                   // Log informational message
  LOG_DETAIL,                 // Log lower-level detail messages
  LOG_DEBUG                   // Log debugging messages
} log_level_t;
  
// Sets the logging level.
void set_log_level(log_level_t level);

// Retrieves the current logging level.
log_level_t log_level();

// Sets the buffering parameters for logging, including the size limit on log messages 
// and the flush frequency.
void set_log_buffering(log_level_t log_type, int size_limit, int num_messages_between_flush);

// Sets the output stream for the given type of log message.
void set_log_stream(log_level_t log_type, FILE* stream);

// Sets the output MPI rank for parallel logging output. Output will only 
// be reported on this rank.
void set_log_mpi_rank(log_level_t log_type, int mpi_rank);

// Returns the output stream for the given type of log message if the 
// type is enabled, NULL if it is not.
FILE* log_stream(log_level_t log_type);

// Issues a debug message.
void log_debug(const char* message, ...);

// Issues a detail message (more probing than informational).
void log_detail(const char* message, ...);

// Issues an informational message.
void log_info(const char* message, ...);

// Issues an urgent.
void log_urgent(const char* message, ...);

#endif
