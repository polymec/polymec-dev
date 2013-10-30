
#
# PREFIX, EPREFIX, prefix and exec_prefix. Lower cases variants have higher
# priority.
#
PREFIX		?= /usr/local
prefix		?= $(PREFIX)

EPREFIX		?= $(prefix)
exec_prefix	?= $(EPREFIX)

#
# Various other installation paths.
#
bindir		?= $(exec_prefix)/bin
sbindir		?= $(exec_prefix)/sbin
libexecdir	?= $(exec_prefix)/libexec
datadir		?= $(prefix)/share
sysconfdir	?= $(prefix)/etc
sharedstatedir	?= $(prefix)/com
localstatedir	?= $(prefix)/var
includedir	?= $(prefix)/include
libdir		?= $(prefix)/lib
localedir	?= $(datadir)/locale
mandir		?= $(datadir)/man
infodir		?= $(datadir)/info


#
# CC
#
CC		?= gcc
CFLAGS		?= -g -pipe
CPPFLAGS	?= -I$(DESTDIR)$(includedir)


#
# AR
#
AR		?= ar
ARFLAGS		?= cru


#
# RANLIB
#
RANLIB		?= ranlib
RANLIBFLAGS	?= 


#
# LD
#
LD		?= ld
LDFLAGS		?= -L$(DESTDIR)$(libdir)


#
# M4
#
M4		?= m4

M4FLAGS		?=


#
# install
#
INSTALL		?= install
INSTALL_DATA	?= $(INSTALL) -m 644
INSTALL_PROGRAM	?= $(INSTALL) -m 755


#
# mkdir
#
MKDIR		?= mkdir


#
# libtool
#
LIBTOOL		?= libtool

