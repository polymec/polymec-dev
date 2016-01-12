// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <stdlib.h>
#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <string.h>
#include "polymec.h"  // High-level polymec header!

// This simply tests whether a program can be built and linked with polymec 
// using the high-level headers.
int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  return 0;
}
