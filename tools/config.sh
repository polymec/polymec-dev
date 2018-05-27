# ^^^^^^ location of polymec source code.

# config.sh -- A CMake configuration script.
# Edit this file to change the parameters in your build. Uncomment exactly one 
# value for each parameter.

#-----------------------------------------------------------------------------
#                             Installation prefix
#-----------------------------------------------------------------------------
PREFIX=$HOME/opt

#-----------------------------------------------------------------------------
#                                   MPI
#-----------------------------------------------------------------------------

# Build with MPI for parallel simulations.
#MPI=ON 

#-----------------------------------------------------------------------------
#                                  OPENMP                 
#-----------------------------------------------------------------------------

# Enable OpenMP threading.
#OPENMP=ON

#-----------------------------------------------------------------------------
#                             Floating Point Precision
#-----------------------------------------------------------------------------
# Choose one of the following.

# Double precision.
PRECISION=double

# Single precision.
#PRECISION=single

#-----------------------------------------------------------------------------
#                                Build type
#-----------------------------------------------------------------------------
# Choose one of the following.

# Debug executable (debugging symbols, no optimization).
BUILD_TYPE=Debug

# Release executable (No symbols, optimization).
#BUILD_TYPE=Release

#-----------------------------------------------------------------------------
#                              Shared libraries
#-----------------------------------------------------------------------------

# Uncomment to build polymec libraries as shared libraries.
#SHARED_LIBS=ON

#-----------------------------------------------------------------------------
#                              Build generator
#-----------------------------------------------------------------------------
# Choose one of the following.

# Good old-fashioned UNIX makefiles, like God intended.
GENERATOR="Unix Makefiles"

# Code::Blocks (with UNIX makefiles underneath).
#GENERATOR="CodeBlocks - Unix Makefiles"

# XCode (shudder).
#GENERATOR="XCode"

# CodeLite (with UNIX makefiles underneath).
#GENERATOR="CodeLite - Unix Makefiles"

# Eclipse CDT4 (with UNIX makefiles underneath).
#GENERATOR="Eclipse CDT4 - Unix Makefiles"

# Kate
#GENERATOR="Kate"

# Sublime Text 2 (with UNIX makefiles underneath).
#GENERATOR="Sublime Text 2 - Unix Makefiles"

#-----------------------------------------------------------------------------
#                               Verbose builds
#-----------------------------------------------------------------------------

# Uncomment this if you want really verbose builds.
#VERBOSE=ON

#-----------------------------------------------------------------------------
#                           Building on Special Machines
#-----------------------------------------------------------------------------

# Uncomment this to indicate that we're building on a special, named machine.
#MACHINE=mymachine

#-----------------------------------------------------------------------------
#                           Code Coverage Analysis
#-----------------------------------------------------------------------------

#COVERAGE=ON

#-----------------------------------------------------------------------------
#                           Continuous Integration
#-----------------------------------------------------------------------------

# Use this if you want to pretend that you're building on Travis CI.
#TRAVIS=ON

#-----------------------------------------------------------------------------
#                                   Compilers
#-----------------------------------------------------------------------------

if [ "$MPI" = "ON" ]; then
  CC=mpicc
  CXX=mpic++
  FC=mpif90
else
  CC=cc
  CXX=c++
  FC=gfortran
fi

# Override compilers here (ONLY if you know what you're doing!).

# C compiler.
#CC=cc

# C++ compiler.
#CXX=c++

# Fortran compiler.
#FC=gfortran

#-----------------------------------------------------------------------------
#                   Don't change anything below here.
#-----------------------------------------------------------------------------

OPTIONS=""
if [ "$MPI" = "ON" ]; then
  OPTIONS="-DHAVE_MPI=ON"
fi
if [ "$OPENMP" = "ON" ]; then
  OPTIONS="$OPTIONS -DUSE_OPENMP=ON"
fi
if [ "$SHARED_LIBS" = "ON" ]; then
  OPTIONS="$OPTIONS -DBUILD_SHARED_LIBS=ON"
fi
if [ "$VERBOSE" = "ON" ]; then
  OPTIONS="$OPTIONS -DCMAKE_VERBOSE_MAKEFILE=ON"
fi
if [ "$TRAVIS" = "ON" ]; then
  OPTIONS="$OPTIONS -DTRAVIS=ON"
fi
if [ "$COVERAGE" = "ON" ]; then
  OPTIONS="$OPTIONS -DCOVERAGE=ON"
fi
if [ ! "$MACHINE" = "" ]; then
  OPTIONS="$OPTIONS -DPOLYMEC_MACHINE=$MACHINE"
fi

cmake \
 -DCMAKE_INSTALL_PREFIX:PATH=$PREFIX \
 -DCMAKE_BUILD_TYPE=$BUILD_TYPE \
 -DCMAKE_C_COMPILER=$CC \
 -DCMAKE_CXX_COMPILER=$CXX \
 -DCMAKE_Fortran_COMPILER=$FC \
 -DPOLYMEC_PRECISION=$PRECISION \
 $OPTIONS \
 -G "$GENERATOR" \
 $SOURCE_DIR

