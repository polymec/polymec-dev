// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_LOGGERS_H
#define POLYMEC_LOGGERS_H

#include <stdio.h>

// Types of log messages / log levels. The levels are increasing in amount 
// of output: LOG_URGENT logs only urgent messages, LOG_INFO logs those and 
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

// Logging output modes to control degrees of output.
typedef enum
{
  LOG_TO_SINGLE_RANK,         // Logs all messages only to a single MPI rank
  LOG_TO_ALL_RANKS            // Logs all messages to all MPI ranks
} log_mode_t;

// Sets the logging model. By default, LOG_TO_SINGLE_RANK (rank 0) is used.
void set_log_mode(log_mode_t mode);

// Retrieves the current logging mode.
log_mode_t log_mode();

// Sets the output MPI rank for parallel logging output. Output will only 
// be reported on this rank if the log mode is set to LOG_TO_SINGLE_RANK.
void set_log_mpi_rank(log_level_t log_type, int rank);

// Retrieves the current buffering parameters for logging, including the size limit 
// on log messages and the flush frequency. If a logger is not found at the 
// given level, -1 is returned for each of these.
void get_log_buffering(log_level_t level, int* message_size_limit, int* num_messages_between_flush);

// Sets the buffering parameters for logging, including the size limit on log messages 
// and the flush frequency.
void set_log_buffering(log_level_t log_type, int size_limit, int num_messages_between_flush);

// Sets the output stream for the given type of log message.
void set_log_stream(log_level_t log_type, FILE* stream);

// Returns the output stream for the given type of log message if the 
// type is enabled, NULL if it is not.
FILE* log_stream(log_level_t log_type);

// Issues a debug message.
void log_debug(const char* message, ...);

// Issues a debug message without formatting. Use this for extremely large
// log messages, or for messages that contain text with formatting characters.
void log_debug_literal(const char* message);

// Issues a detail message (more probing than informational).
void log_detail(const char* message, ...);

// Issues a detail message without formatting. Use this for extremely large
// log messages, or for messages that contain text with formatting characters.
void log_detail_literal(const char* message);

// Issues an informational message.
void log_info(const char* message, ...);

// Issues an informational message without formatting. Use this for extremely 
// large log messages, or for messages that contain text with formatting characters.
void log_info_literal(const char* message);

// Issues an urgent message that cannot be filtered.
void log_urgent(const char* message, ...);

// Issues an urgent message without formatting. Use this for extremely large
// log messages, or for messages that contain text with formatting characters.
void log_urgent_literal(const char* message);

#endif
