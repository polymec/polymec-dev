[![License: MPL 2.0](https://img.shields.io/badge/License-MPL%202.0-brightgreen.svg)](https://opensource.org/licenses/MPL-2.0)
[![Build Status](https://travis-ci.org/polymec/polymec-dev.svg?branch=master)](https://travis-ci.org/polymec/polymec-dev)
[![Coverage Status](https://codecov.io/gh/polymec/polymec-dev/branch/master/graph/badge.svg)](https://codecov.io/gh/polymec/polymec-dev)

# Polymec

Copyright (c) 2012-2019, Jeffrey N. Johnson
All rights reserved.

Polymec is a set of libraries that can be used to construct models for
physical systems. These libraries are designed to allow a computational
scientist to construct low-maintenance science applications with robust
capabilities.

## Features

* A comprehensive set of C tools appropriate for high performance computing,
  including containers, spatial data structures, and interfaces to important
  libraries and facilities.
* A layered structure allowing different levels of buy-in: using only
  basic data structures, adopting a framework for simulator apps, and
  everything between.
* Stiffly-accurate, implicit time integrators that use high-quality
  preconditioners without a lot of application-specific classes.
* A strategy that favors a small code footprint and good scalability, likely
  at the expense of implementations that run optimally on serial computers.

## License

Polymec is licensed under the Mozilla Public License (MPL) version 2.0, which
is defined in the LICENSE file, and also available at http://mozilla.org/MPL/2.0/.
This license allows you to use Polymec in a commercial code as long as
certain requirements are satisfied. Please read or read about the license.

The 3rd-party libraries used by Polymec have their own licenses which are
described within the various source trees in the 3rdparty directory. Source
files have been minimally modified within these trees; typically they contain
files extracted directly from distribution tarballs. You'll find instructions
for obtaining each library in its directory.

## Installing Polymec

You can use polymec on Linux or macOS systems.

### Software Requirements

The polymec libraries require the following software:

* A C compiler that adheres to (or aspires to) the C11 standard.
  GCC 4.9+ and Clang 3.4+ work well.
* A reasonable C++ compiler, for the C++ build/link test
* MPI (OpenMPI or MPICH), for parallelism
* A working set of LAPACK/BLAS libraries, for dense linear algebra
* The Bourne Again SHell (bash), for bootstrapping
* CMake 3.10+, for configuring and generating build files
* GNU Make or Ninja, for performing the actual build
* A recent version of Perl, for HDF5's installation process

The following are helpful, but not required:

* A recent version of Git, for code development
* Valgrind, for finding memory-related issues
* A Python 3.x interpreter, for running some of the targets
* Doxygen, for generating reference documentation)

You can easily install these using your favorite package manager.

### Building

To build polymec on a UNIX-like system, change to your `polymec-dev` directory
and type the following commands:

```
./bootstrap build_dir
```

where `build_dir` is the directory in which you want to build. Then just
follow the onscreen directions: you change to that build directory, edit
`config.sh` to define your build, and then start the build using your
generator's build process. For the default generator (UNIX makefiles), this
is just `make`. For Ninja (recommended if you have it), it's `ninja`.

### Installing

To install polymec, use the install command for the generator you've selected.
For example, if you're using a generator that writes UNIX makefiles, run

```
make install [-j #threads]
```

from your build directory.

### Other Targets

These targets all work with Make and Ninja.

* `test` - Runs all unit tests for the library. Use `ctest -j #threads` instead, though, to run the tests in parallel.
* `memcheck` - Runs all unit tests for the library using Valgrind if you're on Linux and it's available. Tests run in parallel.
* `clean` - Removes all build assets but retains configuration options.
* `stats` - Prints some interesting code statistics to the screen. Requires Python 3.x.
* `coverage` - Generates a code coverage report. This only works if you enable the `COVERAGE` build option in `config.sh`, and if gcov and lcov are available.
* `doc` - Generates HTML and man reference documentation pages using Doxygen. Only available when CMake finds Doxygen.

