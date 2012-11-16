#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <mpi.h>
#include "core/arbi.h"
#include "core/loggers.h"

typedef struct 
{
  int message_size_limit;
  int flush_every;
  char** buffers;
  int message_counter;
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

static void logger_set_buffering(logger_t* logger, int size_limit, int flush_every)
{
  ASSERT(logger != NULL);
  ASSERT(size_limit > 0);
  ASSERT(flush_every > 0);

  // Toss out any old messages.
  if (logger->buffers != NULL)
  {
    for (int i = 0; i < logger->flush_every; ++i)
      free(logger->buffers[i]);
    free(logger->buffers);
  }

  // In with the new.
  logger->flush_every = flush_every;
  logger->message_size_limit = size_limit;
  logger->buffers = malloc(sizeof(char*)*logger->flush_every);
  for (int i = 0; i < logger->flush_every; ++i)
    logger->buffers[i] = malloc(sizeof(char)*size_limit+1);
  logger->message_counter = 0;
}

static void logger_flush(logger_t* logger)
{
  if (logger != NULL)
  {
    if (logger->stream != NULL)
    {
      for (int i = 0; i < logger->message_counter; ++i)
        fprintf(logger->stream, "%s\n", logger->buffers[i]);
    }
    logger->message_counter = 0;
  }
}

static void logger_log(logger_t* logger, char* message)
{
  strncpy(logger->buffers[logger->message_counter], message, logger->message_size_limit);
  logger->message_counter++;
  if (logger->message_counter == logger->flush_every)
    logger_flush(logger);
}

static logger_t* create_logger()
{
  logger_t* logger = malloc(sizeof(logger_t));
  logger->buffers = NULL;
  logger->stream = stdout;
  logger->mpi_rank = 0;
  logger_set_buffering(logger, 1024, 1);
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

void set_log_buffering(log_level_t level, int message_size_limit, int num_messages_between_flush)
{
  ASSERT(num_messages_between_flush > 0);
  logger_t* logger = get_logger(level);
  if (logger != NULL)
    logger_set_buffering(logger, message_size_limit, num_messages_between_flush);
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
    logger_log(logger, m);
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
    logger_log(logger, m);
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
    logger_log(logger, m);
  }
}

#ifdef __cplusplus
}
#endif

