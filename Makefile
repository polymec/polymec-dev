# Makefile -- Use this to build on *NIX systems.

# Options set on command line.
debug      = not-set
mpi        = not-set
precision  = not-set
verbose    = not-set
prefix     = not-set
sanitize   = not-set

# This proxies everything to the builddir cmake.

cputype = $(shell uname -m | sed "s/\\ /_/g")
systype = $(shell uname -s)

BUILDDIR := build/$(systype)-$(cputype)
CONFIG_FLAGS = -DUNIX=1

# Process configuration options.

# Verbose builds?
ifeq ($(verbose), 1)
  CONFIG_FLAGS += -DCMAKE_VERBOSE_MAKEFILE=1
endif

# MPI
ifeq ($(mpi), 1)
  BUILDDIR := ${BUILDDIR}-mpi
  CC = mpicc
  CONFIG_FLAGS += -DHAVE_MPI=1
else
  ifeq ($(CC), )
    CC  = cc
    CXX = c++
  endif
  CONFIG_FLAGS += -DHAVE_MPI=0
endif

# Precision.
ifneq ($(precision), not-set)
  BUILDDIR := ${BUILDDIR}-$(precision)
  CONFIG_FLAGS += -DPOLYMEC_PRECISION=$(precision)
else
  BUILDDIR := ${BUILDDIR}-double
  CONFIG_FLAGS += -DPOLYMEC_PRECISION=double
endif

BUILDDIR := ${BUILDDIR}-`basename ${CC}`
CONFIG_FLAGS += -DCC=${CC} -DCXX=${CXX}

# Debugging symbols
ifneq ($(debug), not-set)
  BUILDDIR := ${BUILDDIR}-Debug
  CONFIG_FLAGS += -DCMAKE_BUILD_TYPE=Debug
else
  BUILDDIR := ${BUILDDIR}-Release
  CONFIG_FLAGS += -DCMAKE_BUILD_TYPE=Release
endif

# Installation prefix.
ifneq ($(prefix), not-set)
  CONFIG_FLAGS += -DCMAKE_INSTALL_PREFIX:PATH=$(prefix)
else
  CONFIG_FLAGS += -DCMAKE_INSTALL_PREFIX:PATH=/usr/local
endif

# Special considerations for specific systems.
ifeq ($(systype), Darwin)
  CONFIG_FLAGS += -DAPPLE=1
else 
  ifeq ($(systype), Linux)
    CONFIG_FLAGS += -DLINUX=1
  endif
endif

# Address sanitizer.
ifeq ($(sanitize), 1)
  BUILDDIR := ${BUILDDIR}-AddressSanitizer
  CONFIG_FLAGS += -DADDRESS_SANITIZER=1
endif

define run-config
@mkdir -p $(BUILDDIR)
@cd $(BUILDDIR) && cmake $(CURDIR) $(CONFIG_FLAGS)
endef

all test clean install:
	@if [ ! -f $(BUILDDIR)/Makefile ]; then \
		more INSTALL; \
	else \
		make -C $(BUILDDIR) $@ --no-print-directory $(MAKEFLAGS); \
	fi

config: distclean
	$(run-config)

distclean:
	@rm -rf $(BUILDDIR)

# For rebuilding the polytope library, which is a common thing to do 
# at the moment.
rebuild-polytope:
	@rm -f $(BUILDDIR)/lib/libpolytope*
	@touch 3rdparty/CMakeLists.txt
	@make -C $(BUILDDIR) $(MAKEFLAGS)

stats: 
	@python tools/gather_stats.py

prepend-license: 
	@python tools/prepend_license.py

#dist:
#	utils/mkdist.sh $(PKGNAME)

.PHONY: config distclean all clean install uninstall 
