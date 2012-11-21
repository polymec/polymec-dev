#ifndef LOGGERS_H
#define LOGGERS_H

#include <stdio.h>

// Types of log messages / log levels. The levels are increasing in amount 
// of output: LOG_WARNING logs only warnings, LOG_INFO logs warnings and 
// informational messages, and so on and so forth.
typedef enum
{
  LOG_NONE,
  LOG_WARNING,
  LOG_INFO,
  LOG_DETAIL,
  LOG_DEBUG
} log_level_t;
  
// Sets the logging level.
void set_log_level(log_level_t level);

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

// Issues an warning.
void log_warning(const char* message, ...);

#ifdef __cplusplus
}
#endif

#endif
