// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_TIMER_H
#define POLYMEC_TIMER_H

// This file contains macros that manipulate timers for performance profiling.
// These timers are intended to help profile functions and tasks with
// measurable amounts of work. If you find yourself fretting over "whether the
// timer code is inlined," you're working at too low a level to use these
// timers and should be timing things manually.
//
// These timers are automatically created and managed in a hierarchy,
// recreating the portion of the call graph that is instrumented with them.
//
// The design of these timers was inspired by Brian Van Straalen's timers
// in Chombo, with helpful comments from Noel Keen.

/// \addtogroup core core
///@{

//------------------------------------------------------------------------
//                              Named timers
//------------------------------------------------------------------------
///@{
/// These macros allow the creation and use of a timer with a given name, that
/// occupies the given symbol/variable (unique to its scope). It can be started
/// and stopped any number of times, though starting a started timer or
/// stopping a stopped timer produces a run-time error.
#define CREATE_TIMER(name, symbol) \
  polymec_timer_t* symbol = polymec_timer_get(name)
#define START_TIMER(symbol) \
  polymec_timer_start(symbol)
#define STOP_TIMER(symbol) \
  polymec_timer_stop(symbol)
///@}

//------------------------------------------------------------------------
//                          Function-level timers
//------------------------------------------------------------------------
///@{
/// These macros start and stop a timer intended to measure the time that
/// elapses during a call to a function. Place START_FUNCTION_TIMER() at the
/// beginning of a function body and STOP_FUNCTION_TIMER() at the end
/// (just before a return value, say).
#define START_FUNCTION_TIMER() \
  CREATE_TIMER(__func__, polymec_function_timer); \
  START_TIMER(polymec_function_timer)
#define STOP_FUNCTION_TIMER() STOP_TIMER(polymec_function_timer)
///@}

//------------------------------------------------------------------------
//                          Internal machinery.
//------------------------------------------------------------------------
// Please do not use any of this stuff directly.
void polymec_enable_timers(void);
const char* polymec_timer_file(void);
void polymec_set_timer_file(const char* timer_file);
typedef struct polymec_timer_t polymec_timer_t;
polymec_timer_t* polymec_timer_get(const char* name);
void polymec_timer_start(polymec_timer_t* timer);
void polymec_timer_stop(polymec_timer_t* timer);
void polymec_timer_stop_all(void);
void polymec_timer_report(void);

///@}

#endif
