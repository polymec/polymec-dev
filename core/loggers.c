#include <stdlib.h>
#include <stdarg.h>
#include <mpi.h>
#include "core/arbi.h"
#include "core/loggers.h"

typedef struct 
{
  int message_size_limit;
  int flush_every;
  FILE* stream;
  int mpi_rank;
} logger_t;

static bool first_time = true;
static log_level_t logging_level = LOG_INFO;
static logger_t* debug_logger = NULL;
static logger_t* info_logger = NULL;
static logger_t* warning_logger = NULL;

static void delete_loggers()
{
  if (debug_logger != NULL)
  {
    free(debug_logger);
    debug_logger = NULL;
  }
  if (info_logger != NULL)
  {
    free(info_logger);
    info_logger = NULL;
  }
  if (warning_logger != NULL)
  {
    free(warning_logger);
    warning_logger = NULL;
  }
}

static logger_t* create_logger()
{
  logger_t* logger = malloc(sizeof(logger_t));
  logger->message_size_limit = 1024;
  logger->flush_every = 32;
  logger->stream = stdout;
  logger->mpi_rank = 0;
  if (first_time)
  {
    arbi_atexit(delete_loggers);
    first_time = false;
  }
  return logger;
}

void set_log_level(log_level_t level)
{
  logging_level = level;
}

static logger_t* get_logger(log_level_t level)
{
  logger_t* logger = NULL;
  switch (level) 
  {
    case LOG_DEBUG:
      if (debug_logger == NULL)
        debug_logger = create_logger();
      logger = debug_logger;
      break;
    case LOG_INFO:
      if (info_logger == NULL)
        info_logger = create_logger();
      logger = info_logger;
      break;
    case LOG_WARNING:
      if (info_logger == NULL)
        info_logger = create_logger();
      logger = info_logger;
      break;
    case LOG_NONE:
      break;
  }
  return logger;
}

void set_log_message_size_limit(log_level_t level, int size_limit)
{
  ASSERT(size_limit > 0);
  logger_t* logger = get_logger(level);
  if (logger != NULL)
    logger->message_size_limit = size_limit;
}

void set_log_flush_period(log_level_t level, int num_messages_between_flush)
{
  ASSERT(num_messages_between_flush > 0);
  logger_t* logger = get_logger(level);
  if (logger != NULL)
    logger->flush_every = num_messages_between_flush;
}

void set_log_stream(log_level_t log_type, FILE* stream)
{
  logger_t* logger = get_logger(log_type);
  if (logger != NULL)
    logger->stream = stream;
}

void set_log_mpi_rank(log_level_t log_type, int mpi_rank)
{
  ASSERT(mpi_rank >= 0);
  int nprocs = MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  ASSERT(mpi_rank < nprocs);
  logger_t* logger = get_logger(log_type);
  if (logger != NULL)
    logger->mpi_rank = mpi_rank;
}

void log_debug(const char* message, ...)
{
  logger_t* logger = get_logger(LOG_DEBUG);
  if (logging_level < LOG_DEBUG) return;
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank != logger->mpi_rank) return;
  if (logger->stream != NULL)
  {
    // Extract the variadic arguments and splat them into a string.
    char m[logger->message_size_limit];
    va_list argp;
    va_start(argp, message);
    vsnprintf(m, logger->message_size_limit, message, argp);
    va_end(argp);

    // Log it.
    fprintf(logger->stream, "%s\n", m);
  }
}

void log_info(const char* message, ...)
{
  logger_t* logger = get_logger(LOG_INFO);
  if (logging_level < LOG_INFO) return;
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank != logger->mpi_rank) return;
  if (logger->stream != NULL)
  {
    // Extract the variadic arguments and splat them into a string.
    char m[logger->message_size_limit];
    va_list argp;
    va_start(argp, message);
    vsnprintf(m, logger->message_size_limit, message, argp);
    va_end(argp);

    // Log it.
    fprintf(logger->stream, "%s\n", m);
  }
}

void log_warning(const char* message, ...)
{
  logger_t* logger = get_logger(LOG_WARNING);
  if (logging_level < LOG_WARNING) return;
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank != logger->mpi_rank) return;
  if (logger->stream != NULL)
  {
    // Extract the variadic arguments and splat them into a string.
    char m[logger->message_size_limit];
    va_list argp;
    va_start(argp, message);
    vsnprintf(m, logger->message_size_limit, message, argp);
    va_end(argp);

    // Log it.
    fprintf(logger->stream, "%s\n", m);
  }
}

#ifdef __cplusplus
}
#endif

