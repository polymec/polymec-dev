// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// Copyright (c) 2012-2017, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <string.h>
#include "cmocka.h"
#include "core/polymec.h"
#include "core/memory_info.h"

static void test_set_log_level(void** state)
{
  log_level_t levels[] = {LOG_NONE, LOG_URGENT, LOG_INFO, LOG_DETAIL, LOG_DEBUG};
  for (int i = 0; i < 5; ++i)
  {
    set_log_level(levels[i]);
    assert_true(log_level() == levels[i]);
    log_urgent("Log message");
    log_info("Log message");
    log_detail("Log message");
    log_debug("Log message");
  }
}

static void test_set_log_mode(void** state)
{
  log_mode_t modes[] = {LOG_TO_SINGLE_RANK, LOG_TO_ALL_RANKS};
  log_level_t levels[] = {LOG_NONE, LOG_URGENT, LOG_INFO, LOG_DETAIL, LOG_DEBUG};
  for (int i = 0; i < 2; ++i)
  {
    for (int l = 0; l < 5; ++l)
    {
      if (modes[i] == LOG_TO_SINGLE_RANK)
        set_log_mpi_rank(levels[l], 0);
      set_log_mode(modes[i]);
      assert_true(log_mode() == modes[i]);
      log_urgent("Log message");
      log_info("Log message");
      log_detail("Log message");
      log_debug("Log message");
    }
  }
}

static void test_set_log_buffering(void** state)
{
  log_level_t levels[] = {LOG_NONE, LOG_URGENT, LOG_INFO, LOG_DETAIL, LOG_DEBUG};
  for (int i = 0; i < 5; ++i)
  {
    int orig_size_limit, orig_buffer_len;
    get_log_buffering(levels[i], &orig_size_limit, &orig_buffer_len);
    set_log_buffering(levels[i], 1024, 10);
    int new_size_limit, new_buffer_len;
    get_log_buffering(levels[i], &new_size_limit, &new_buffer_len);
    assert_int_equal(1024, new_size_limit);
    assert_int_equal(10, new_buffer_len);
    log_urgent("Log message");
    log_info("Log message");
    log_detail("Log message");
    log_debug("Log message");
    set_log_buffering(levels[i], orig_size_limit, orig_buffer_len);
  }
}

static void test_set_log_stream(void** state)
{
  log_level_t levels[] = {LOG_NONE, LOG_URGENT, LOG_INFO, LOG_DETAIL, LOG_DEBUG};
  for (int i = 0; i < 5; ++i)
  {
    FILE* orig_stream = log_stream(levels[i]);
    assert_true(orig_stream == stdout);
    set_log_stream(levels[i], stderr);
    assert_true(log_stream(levels[i]) == stderr);
    log_urgent("Log message");
    log_info("Log message");
    log_detail("Log message");
    log_debug("Log message");
    set_log_stream(levels[i], stdout);
  }
}

static void test_set_log_literals(void** state)
{
  log_level_t levels[] = {LOG_NONE, LOG_URGENT, LOG_INFO, LOG_DETAIL, LOG_DEBUG};
  for (int i = 0; i < 5; ++i)
  {
    set_log_level(levels[i]);
    log_urgent_literal("Log message");
    log_info_literal("Log message");
    log_detail_literal("Log message");
    log_debug_literal("Log message");
  }
  set_log_level(LOG_URGENT);
}

static void test_set_log_indentation(void** state)
{
  log_level_t levels[] = {LOG_NONE, LOG_URGENT, LOG_INFO, LOG_DETAIL, LOG_DEBUG};
  for (int i = 0; i < 5; ++i)
  {
    set_log_indentation_prefix(levels[i], "  ");
    log_indent(levels[i]);
    log_urgent("Log message");
    log_info("Log message");
    log_detail("Log message");
    log_debug("Log message");
    log_unindent(levels[i]);
  }
}

static void test_debug_during_exit(void** state)
{
  set_log_level(LOG_DEBUG);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const struct CMUnitTest tests[] = 
  {
    cmocka_unit_test(test_set_log_level),
    cmocka_unit_test(test_set_log_mode),
    cmocka_unit_test(test_set_log_buffering),
    cmocka_unit_test(test_set_log_stream),
    cmocka_unit_test(test_set_log_literals),
    cmocka_unit_test(test_set_log_indentation),
    cmocka_unit_test(test_debug_during_exit)
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
