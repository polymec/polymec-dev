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
  LOG_DEBUG
} log_level_t;
  
// Sets the logging level.
void set_log_level(log_level_t level);

// Sets the size limit on log messages for the logger (in bytes).
void set_log_message_size_limit(log_level_t level, int size_limit);

// Sets the flushing period (in number of messages) for the given log level.
void set_log_flush_period(log_level_t level, int num_messages_between_flush);

// Sets the output stream for the given type of log message.
void set_log_stream(log_level_t log_type, FILE* stream);

// Issues a debug message.
void log_debug(const char* message, ...);

// Issues an informational message.
void log_info(const char* message, ...);

// Issues an warning.
void log_warning(const char* message, ...);

#ifdef __cplusplus
}
#endif

#endif
