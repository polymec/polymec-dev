#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <mpi.h>
#include "core/polymec.h"
#include "core/loggers.h"

typedef struct 
{
  int message_size_limit;
  int flush_every;
  char* buffer;
  int message_counter;
  FILE* stream;
  int mpi_rank;
} logger_t;

static bool first_time = true;
static log_level_t logging_level = LOG_INFO;
static logger_t* loggers[] = {NULL, NULL, NULL, NULL, NULL};

static void delete_loggers()
{
  for (int i = 0; i < 5; ++i)
  {
    if (loggers[i] != NULL)
    {
      free(loggers[i]);
      loggers[i] = NULL;
    }
  }
}

static void logger_set_buffering(logger_t* logger, int size_limit, int flush_every)
{
  ASSERT(logger != NULL);
  ASSERT(size_limit > 0);
  ASSERT(flush_every > 0);

  // Flush what we've got.

  // In with the new.
  logger->flush_every = flush_every;
  logger->message_size_limit = size_limit;
  logger->buffer = realloc(logger->buffer, sizeof(char*)*(logger->flush_every*(logger->message_size_limit+1)));
  logger->message_counter = 0;
}

static void logger_flush(logger_t* logger)
{
  if (logger != NULL)
  {
    if (logger->stream != NULL)
    {
      for (int i = 0; i < logger->message_counter; ++i)
        fprintf(logger->stream, "%s\n", &logger->buffer[i*(logger->message_size_limit+1)]);
    }
    logger->message_counter = 0;
  }
}

static void logger_log(logger_t* logger, char* message)
{
  strncpy(&logger->buffer[logger->message_counter*(logger->message_size_limit+1)], message, logger->message_size_limit);
  logger->message_counter++;
  if (logger->message_counter == logger->flush_every)
    logger_flush(logger);
}

static logger_t* create_logger()
{
  logger_t* logger = malloc(sizeof(logger_t));
  logger->buffer = NULL;
  logger->stream = stdout;
  logger->mpi_rank = 0;
  logger_set_buffering(logger, 1024, 1);
  if (first_time)
  {
    polymec_atexit(delete_loggers);
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
  if (loggers[level] == NULL)
    loggers[level] = create_logger();
  return loggers[level];
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

void log_detail(const char* message, ...)
{
  logger_t* logger = get_logger(LOG_DETAIL);
  if (logging_level < LOG_DETAIL) return;
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

void log_urgent(const char* message, ...)
{
  logger_t* logger = get_logger(LOG_URGENT);
  if (logging_level < LOG_URGENT) return;
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

