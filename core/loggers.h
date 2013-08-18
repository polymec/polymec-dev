// Copyright 2012-2013 Jeffrey Johnson.
// 
// This file is part of Polymec, and is licensed under the Apache License, 
// Version 2.0 (the "License"); you may not use this file except in 
// compliance with the License. You may may find the text of the license in 
// the LICENSE file at the top-level source directory, or obtain a copy of 
// it at
// 
// http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

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

// Issues a debug message.
void log_debug(const char* message, ...);

// Issues a detail message (more probing than informational).
void log_detail(const char* message, ...);

// Issues an informational message.
void log_info(const char* message, ...);

// Issues an urgent.
void log_urgent(const char* message, ...);

#endif
