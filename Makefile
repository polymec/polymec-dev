# Makefile -- Use this to build on *NIX systems.

# Options set on command line.
prefix     = not-set
verbose    = not-set
precision  = not-set
mpi        = not-set
openmp     = not-set
debug      = not-set
shared     = not-set
machine    = not-set
travis     = not-set
coverage   = not-set
ide        = not-set
build      = not-set
CC         = not-set
CXX        = not-set
FC         = not-set

# This dispatches everything to the builddir cmake.

cputype = $(shell uname -m | sed "s/\\ /_/g")
systype = $(shell uname -s)

BUILDDIR := build/$(systype)-$(cputype)
CONFIG_FLAGS = -DUNIX=1 -Wno-dev

# Process configuration options.

# Verbose builds?
ifeq ($(verbose), 1)
  CONFIG_FLAGS += -DCMAKE_VERBOSE_MAKEFILE=1
endif

# MPI
ifeq ($(mpi), 1)
  BUILDDIR := ${BUILDDIR}-mpi
  CC = mpicc
  CXX = mpicxx
  FC = mpif90
  CONFIG_FLAGS += -DHAVE_MPI=1
else
  BUILDDIR := ${BUILDDIR}-serial
  ifeq ($(CC), not-set)
    CC  = cc
  endif
  ifeq ($(CXX), not-set)
    CXX = c++
  endif
  ifeq ($(FC), not-set)
    FC = gfortran
  endif
  CONFIG_FLAGS += -DHAVE_MPI=0
endif

# OpenMP -- try to enable by default.
ifneq ($(openmp), 0)
  CONFIG_FLAGS += -DUSE_OPENMP=1
endif

# Shared libs?
ifeq ($(shared), 1)
  BUILDDIR := ${BUILDDIR}-shared
  CONFIG_FLAGS += -DBUILD_SHARED_LIBS=ON
else
  BUILDDIR := ${BUILDDIR}-static
  CONFIG_FLAGS += -DBUILD_SHARED_LIBS=OFF
endif

# Machine configuration.
ifneq ($(machine), not-set)
  BUILDDIR := build/$(machine)
  CONFIG_FLAGS += -DPOLYMEC_MACHINE=$(machine)
endif

# Precision.
ifneq ($(precision), not-set)
  BUILDDIR := ${BUILDDIR}-$(precision)
  CONFIG_FLAGS += -DPOLYMEC_PRECISION=$(precision)
else
  BUILDDIR := ${BUILDDIR}-double
  CONFIG_FLAGS += -DPOLYMEC_PRECISION=double
endif

# Code coverage.
ifeq ($(coverage), 1)
  BUILDDIR := ${BUILDDIR}-gcov
  CONFIG_FLAGS += -DCOVERAGE=1
endif

BUILDDIR := ${BUILDDIR}-`basename ${CC}`
CONFIG_FLAGS += -DCC=${CC} -DCXX=${CXX}
ifneq ($(FC), )
  CONFIG_FLAGS += -DFC=${FC}
endif

# Debugging symbols
ifeq ($(debug), not-set)
  BUILDDIR := ${BUILDDIR}-Release
  CONFIG_FLAGS += -DCMAKE_BUILD_TYPE=Release
else
  ifeq ($(debug), 0)
    BUILDDIR := ${BUILDDIR}-Release
    CONFIG_FLAGS += -DCMAKE_BUILD_TYPE=Release
  else 
    BUILDDIR := ${BUILDDIR}-Debug
    CONFIG_FLAGS += -DCMAKE_BUILD_TYPE=Debug
  endif
endif

# Override the build name if we're given one.
ifneq ($(build), not-set)
  BUILDDIR := build/$(build)
endif

# Installation prefix.
ifneq ($(prefix), not-set)
  CONFIG_FLAGS += -DCMAKE_INSTALL_PREFIX:PATH=$(prefix)
else
  # No default installation dir.
  CONFIG_FLAGS += -DCMAKE_INSTALL_PREFIX:PATH=INSTALL_DISABLED
endif

# Special considerations for specific systems.
ifeq ($(systype), Darwin)
  CONFIG_FLAGS += -DAPPLE=1
else 
  ifeq ($(systype), Linux)
    CONFIG_FLAGS += -DLINUX=1
  endif
endif

# Travis-CI-specific settings.
ifeq ($(travis), 1)
  CONFIG_FLAGS += -DTRAVIS_CI=1
endif

# Integrated Development Environment (IDE)?
ifeq ($(ide), codeblocks)
  CONFIG_FLAGS += -G "CodeBlocks - Unix Makefiles"
else
  ifeq ($(ide), xcode)
    CONFIG_FLAGS += -GXcode
  else
    ifeq ($(ide), codelite)
      CONFIG_FLAGS += -G "CodeLite - Unix Makefiles"
    else
      ifeq ($(ide), eclipse)
        CONFIG_FLAGS += -G "Eclipse CDT4 - Unix Makefiles"
      else
        ifeq ($(ide), kate)
          CONFIG_FLAGS += -GKate 
        else
          ifeq ($(ide), sublime)
            CONFIG_FLAGS += -G "Sublime Text 2 - Unix Makefiles" 
          endif
        endif
      endif
    endif
  endif
endif

define run-config
@mkdir -p $(BUILDDIR)
@cd $(BUILDDIR) && cmake $(CURDIR) $(CONFIG_FLAGS)
endef

all clean install:
	@if [ ! -f $(BUILDDIR)/Makefile ]; then \
		more INSTALL; \
	else \
		$(MAKE) -C $(BUILDDIR) $@ --no-print-directory $(MAKEFLAGS); \
	fi

test:
	@if [ ! -f $(BUILDDIR)/Makefile ]; then \
		more INSTALL; \
	else \
		$(MAKE) -C $(BUILDDIR) $@ --no-print-directory $(MAKEFLAGS); \
	  if [ "$(coverage)" = "1" ]; then \
      echo "Writing code coverage report to lcov.info..."; \
	    lcov --base-directory $(BUILDDIR) --directory $(BUILDDIR) -c -q -o lcov.info; \
	    lcov --remove lcov.info '*/build/*' '*valgrind.h' '*/tests/*' -o lcov.info; \
      echo "Generating HTML coverage report -> ./coverage..."; \
      genhtml -o coverage lcov.info; \
    fi \
  fi

memcheck:
	@cd $(BUILDDIR) && ctest -T memcheck

config: distclean
	$(run-config)

distclean:
	@rm -rf $(BUILDDIR)

stats: 
	@python tools/gather_stats.py

prepend-license: 
	@python tools/prepend_license.py

ctags-emacs:
	@ctags -e -f ETAGS -R --exclude=.git --exclude=build 

#dist:
#	utils/mkdist.sh $(PKGNAME)

.PHONY: config distclean all clean install uninstall 
