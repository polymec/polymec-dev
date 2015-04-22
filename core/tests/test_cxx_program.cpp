// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <string.h>
#include "cmockery.h"
#include "polymec.h"

// This just tests the ability of the Polymec headers to successfully compile 
// using a C++ compiler.
int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  return 0;
}
