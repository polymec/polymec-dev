#
# Get CURDIR (GNU Make) from .CURDIR (pmake)
#
CURDIR	?= $(.CURDIR)


#
# Setup our targets
#
all:

.PHONY: clean uninstall


#
# Setup our build environment.
#
include $(CURDIR)/mk/Variables.mk

Prologue	?= $(CURDIR)/mk/Prologue.mk
Epilogue	?= $(CURDIR)/mk/Epilogue.mk


#
# Include local Rules.mk
#
dir	:= $(CURDIR)
include $(dir)/Rules.mk


#
# Debian package management.
#
DPKG_BUILDPACKAGE ?= dpkg-buildpackage
FAKEROOT          ?= fakeroot

.PHONY: debian debian-clean

debian:
	$(DPKG_BUILDPACKAGE) -rfakeroot -I.svn $(DPKG_BUILDPACKAGE_OPTIONS)

debian-clean:
	-$(FAKEROOT) ./debian/rules clean

