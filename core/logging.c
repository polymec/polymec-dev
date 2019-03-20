// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include "core/polymec.h"

typedef struct
{
  // Message buffering.
  int message_size_limit;
  int flush_every;
  char* buffer;
  int message_counter;

  // File stream and MPI info.
  FILE* stream;
  int mpi_rank;

  // Indentation.
  char* indent_prefix;
  int indent_level;
} logger_t;

static bool first_time = true;
static log_level_t logging_level = LOG_INFO;
static log_mode_t logging_mode = LOG_TO_SINGLE_RANK;
static logger_t* loggers[] = {NULL, NULL, NULL, NULL, NULL};

// MPI stuff.
static int mpi_nproc = -1;
static int mpi_rank = -1;
static int max_mpi_rank_digits = -1;

static void logger_free(logger_t* logger)
{
  if (logger->stream != NULL)
    fclose(logger->stream);
  polymec_free(logger->buffer);
  if (logger->indent_prefix != NULL)
    string_free(logger->indent_prefix);
  polymec_free(logger);
}

static void delete_loggers()
{
  for (int i = 0; i < 5; ++i)
  {
    if (loggers[i] != NULL)
    {
      logger_free(loggers[i]);
      loggers[i] = NULL;
    }
  }
}

static void logger_set_buffering(logger_t* logger, int size_limit, int flush_every)
{
  ASSERT(logger != NULL);
  ASSERT(size_limit > 0);
  ASSERT(flush_every > 0);

  // In with the new.
  logger->flush_every = flush_every;
  logger->message_size_limit = size_limit;
  logger->buffer = polymec_realloc(logger->buffer, sizeof(char)*(logger->flush_every*(logger->message_size_limit+1)));
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

static void logger_log(logger_t* logger, const char* message)
{
  // Apply indenting to the message.
  size_t orig_message_len = strlen(message);
  size_t prefix_len = strlen(logger->indent_prefix);
  size_t indentation_len = prefix_len * logger->indent_level;
  char indented_message[indentation_len + orig_message_len + 1];
  if (indentation_len > 0)
  {
    int k = 0;
    for (int i = 0; i < logger->indent_level; ++i)
      for (size_t j = 0; j < prefix_len; ++j, ++k)
        indented_message[k] = logger->indent_prefix[j];
    strcpy(&indented_message[k], message);
  }
  else
    strcpy(indented_message, message);

  if (logging_mode == LOG_TO_ALL_RANKS)
  {
    size_t message_len = (size_t)max_mpi_rank_digits + 5 + strlen(indented_message);
    char message_with_rank[message_len+1];
    switch(max_mpi_rank_digits)
    {
      case 1: snprintf(message_with_rank, message_len, "%1d: %s", mpi_rank, indented_message); break;
      case 2: snprintf(message_with_rank, message_len, "%02d: %s", mpi_rank, indented_message); break;
      case 3: snprintf(message_with_rank, message_len, "%03d: %s", mpi_rank, indented_message); break;
      case 4: snprintf(message_with_rank, message_len, "%04d: %s", mpi_rank, indented_message); break;
      case 5: snprintf(message_with_rank, message_len, "%05d: %s", mpi_rank, indented_message); break;
      case 6: snprintf(message_with_rank, message_len, "%06d: %s", mpi_rank, indented_message); break;
      default: snprintf(message_with_rank, message_len, "%d: %s", mpi_rank, indented_message); break;
    }
    strncpy(&logger->buffer[logger->message_counter*(logger->message_size_limit+1)], message_with_rank, logger->message_size_limit);
  }
  else
    strncpy(&logger->buffer[logger->message_counter*(logger->message_size_limit+1)], indented_message, logger->message_size_limit);
  logger->message_counter++;
  if (logger->message_counter == logger->flush_every)
    logger_flush(logger);
}

static logger_t* logger_new()
{
  logger_t* logger = polymec_malloc(sizeof(logger_t));
  logger->buffer = NULL;
  logger->stream = stdout;
  logger->mpi_rank = 0;
  logger_set_buffering(logger, 1024, 1);
  if (first_time)
  {
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    max_mpi_rank_digits = (int)(log10(mpi_nproc)) + 1;
    polymec_atexit(delete_loggers);
    first_time = false;
  }

  logger->indent_level = 0;
  logger->indent_prefix = string_dup(" ");

  return logger;
}

void set_log_level(log_level_t level)
{
  logging_level = level;
}

log_level_t log_level()
{
  return logging_level;
}

static logger_t* get_logger(log_level_t level)
{
  if (loggers[level] == NULL)
    loggers[level] = logger_new();
  return loggers[level];
}

void set_log_mode(log_mode_t mode)
{
  logging_mode = mode;
}

log_mode_t log_mode()
{
  return logging_mode;
}

void get_log_buffering(log_level_t level, int* message_size_limit, int* num_messages_between_flush)
{
  logger_t* logger = get_logger(level);
  if (logger != NULL)
  {
    *message_size_limit = logger->message_size_limit;
    *num_messages_between_flush = logger->flush_every;
  }
  else
  {
    *message_size_limit = -1;
    *num_messages_between_flush = -1;
  }
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
  if (((mpi_rank == logger->mpi_rank) ||
       (logging_mode == LOG_TO_ALL_RANKS)) &&
      (logger != NULL))
  {
    if ((logger->stream != NULL) &&
        (logger->stream != stdout) && // don't close any system streams(!)
        (logger->stream != stderr))
      fclose(logger->stream);
    logger->stream = stream;
  }
}

void set_log_mpi_rank(log_level_t log_type, int rank)
{
  ASSERT(rank >= 0);
  ASSERT(rank < mpi_nproc);
  logger_t* logger = get_logger(log_type);
  if (logger != NULL)
    logger->mpi_rank = rank;
}

int log_mpi_rank(log_level_t log_type)
{
  if (log_mode() == LOG_TO_SINGLE_RANK)
  {
    logger_t* logger = get_logger(log_type);
    if (logger != NULL)
      return logger->mpi_rank;
    else
      return -1;
  }
  else
    return -1;
}

FILE* log_stream(log_level_t log_type)
{
  if (logging_level < log_type)
    return NULL;
  else
  {
    logger_t* logger = get_logger(log_type);
    if (logger != NULL)
      return logger->stream;
    else
      return NULL;
  }
}

void log_flush(log_level_t log_type)
{
  logger_t* logger = get_logger(log_type);
  if ((logger != NULL) && (logger->stream != NULL))
  {
    logger_flush(logger);
    fflush(logger->stream);
  }
}

void set_log_indentation_prefix(log_level_t log_type, const char* prefix)
{
  logger_t* logger = get_logger(log_type);
  if (logger != NULL)
  {
    if (logger->indent_prefix != NULL)
      string_free(logger->indent_prefix);
    logger->indent_prefix = string_dup(prefix);
  }
}

void log_indent(log_level_t log_type)
{
  logger_t* logger = get_logger(log_type);
  if (logger != NULL)
    logger->indent_level++;
}

void log_unindent(log_level_t log_type)
{
  logger_t* logger = get_logger(log_type);
  if ((logger != NULL) && (logger->indent_level > 0))
    logger->indent_level--;
}

void log_debug(const char* message, ...)
{
  logger_t* logger = get_logger(LOG_DEBUG);
  if (logging_level < LOG_DEBUG) return;
  if ((logging_mode == LOG_TO_SINGLE_RANK) &&
      (mpi_rank != logger->mpi_rank)) return;
  if (logger->stream != NULL)
  {
    // Extract the variadic arguments and splat them into a string.
    char m[logger->message_size_limit+1];
    va_list argp;
    va_start(argp, message);
    vsnprintf(m, logger->message_size_limit, message, argp);
    va_end(argp);

    // Log it.
    logger_log(logger, m);
  }
}

void log_debug_literal(const char* message)
{
  logger_t* logger = get_logger(LOG_DEBUG);
  if (logging_level < LOG_DEBUG) return;
  if ((logging_mode == LOG_TO_SINGLE_RANK) &&
      (mpi_rank != logger->mpi_rank)) return;
  if (logger->stream != NULL)
  {
    // Log it directly.
    logger_log(logger, message);
  }
}

void log_detail(const char* message, ...)
{
  logger_t* logger = get_logger(LOG_DETAIL);
  if (logging_level < LOG_DETAIL) return;
  if ((logging_mode == LOG_TO_SINGLE_RANK) &&
      (mpi_rank != logger->mpi_rank)) return;
  if (logger->stream != NULL)
  {
    // Extract the variadic arguments and splat them into a string.
    char m[logger->message_size_limit+1];
    va_list argp;
    va_start(argp, message);
    vsnprintf(m, logger->message_size_limit, message, argp);
    va_end(argp);

    // Log it.
    logger_log(logger, m);
  }
}

void log_detail_literal(const char* message)
{
  logger_t* logger = get_logger(LOG_DETAIL);
  if (logging_level < LOG_DETAIL) return;
  if ((logging_mode == LOG_TO_SINGLE_RANK) &&
      (mpi_rank != logger->mpi_rank)) return;
  if (logger->stream != NULL)
  {
    // Log it directly.
    logger_log(logger, message);
  }
}

void log_info(const char* message, ...)
{
  logger_t* logger = get_logger(LOG_INFO);
  if (logging_level < LOG_INFO) return;
  if ((logging_mode == LOG_TO_SINGLE_RANK) &&
      (mpi_rank != logger->mpi_rank)) return;
  if (logger->stream != NULL)
  {
    // Extract the variadic arguments and splat them into a string.
    char m[logger->message_size_limit+1];
    va_list argp;
    va_start(argp, message);
    vsnprintf(m, logger->message_size_limit, message, argp);
    va_end(argp);

    // Log it.
    logger_log(logger, m);
  }
}

void log_info_literal(const char* message)
{
  logger_t* logger = get_logger(LOG_INFO);
  if (logging_level < LOG_INFO) return;
  if ((logging_mode == LOG_TO_SINGLE_RANK) &&
      (mpi_rank != logger->mpi_rank)) return;
  if (logger->stream != NULL)
  {
    // Log it directly.
    logger_log(logger, message);
  }
}

void log_urgent(const char* message, ...)
{
  logger_t* logger = get_logger(LOG_URGENT);
  if (logging_level < LOG_URGENT) return;
  if ((logging_mode == LOG_TO_SINGLE_RANK) &&
      (mpi_rank != logger->mpi_rank)) return;
  if (logger->stream != NULL)
  {
    // Extract the variadic arguments and splat them into a string.
    char m[logger->message_size_limit+1];
    va_list argp;
    va_start(argp, message);
    vsnprintf(m, logger->message_size_limit, message, argp);
    va_end(argp);

    // Log it.
    logger_log(logger, m);
  }
}

void log_urgent_literal(const char* message)
{
  logger_t* logger = get_logger(LOG_URGENT);
  if (logging_level < LOG_URGENT) return;
  if ((logging_mode == LOG_TO_SINGLE_RANK) &&
      (mpi_rank != logger->mpi_rank)) return;
  if (logger->stream != NULL)
  {
    // Log it directly.
    logger_log(logger, message);
  }
}

