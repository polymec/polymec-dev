# Non-recursive prologue.
include $(Prologue)


VERSION_$(d)	:= 0.3	# libtool doesn't like the third [patch] component.

OBJS_$(d)	:= $(d)/arena.o $(d)/pool.o $(d)/proto.o $(d)/util.o
SRCS_$(d)	:= $(d)/arena.c $(d)/pool.c $(d)/proto.c $(d)/util.c
INCS_$(d)	:= $(d)/arena.h $(d)/pool.h $(d)/proto.h $(d)/rbits.h $(d)/util.h $(d)/align.h

CFLAGS_$(d)	:= $(CFLAGS)
CPPFLAGS_$(d)	:= -DLIBARENA_SOURCE $(CPPFLAGS)


#
# Compilation Targets
#
Build_$(d)	= $(CC) -c $(CPPFLAGS_$(@D)) $(CFLAGS_$(@D)) -o $@ $(@:.o=.c)

$(d)/arena.o: $(d)/arena.c
	$(Build_$(@D))

$(d)/pool.o: $(d)/pool.c
	$(Build_$(@D))

$(d)/proto.o: $(d)/util.c
	$(Build_$(@D))

$(d)/util.o: $(d)/util.c
	$(Build_$(@D))

$(d)/libarena.a: $(OBJS_$(d))
	$(AR) $(ARFLAGS) $@ $^ $>
	$(RANLIB) $(RANLIBFLAGS) $@

all: $(d)/libarena.a


#
# PIC Compilation Targets
#
BuildPIC_$(d)	= $(LIBTOOL) --mode=compile $(CC) -c $(CPPFLAGS_$(@D)) $(CFLAGS_$(@D)) -o $@ $(@:.lo=.c)

$(d)/arena.lo: $(d)/arena.c
	$(BuildPIC_$(@D))

$(d)/pool.lo: $(d)/pool.c
	$(BuildPIC_$(@D))

$(d)/proto.lo: $(d)/util.c
	$(BuildPIC_$(@D))

$(d)/util.lo: $(d)/util.c
	$(BuildPIC_$(@D))

$(d)/libarena.la: $(OBJS_$(d):.o=.lo)
	$(LIBTOOL) --mode=link $(CC) -rpath $(libdir) -version-info $$(echo $(VERSION_$(@D)) | tr '.' ':') -o $@ $^ $>

pic: $(d)/libarena.la


#
# Installation Targets
#
$(DESTDIR)$(includedir)/arena/proto.h: $(d)/proto.h
	$(MKDIR) -p $(@D)
	$(INSTALL_DATA) $^ $> $@

$(DESTDIR)$(includedir)/arena/arena.h: $(d)/arena.h
	$(MKDIR) -p $(@D)
	$(INSTALL_DATA) $^ $> $@

$(DESTDIR)$(includedir)/arena/align.h: $(d)/align.h
	$(MKDIR) -p $(@D)
	$(INSTALL_DATA) $^ $> $@

$(DESTDIR)$(includedir)/arena/pool.h: $(d)/pool.h
	$(MKDIR) -p $(@D)
	$(INSTALL_DATA) $^ $> $@

$(DESTDIR)$(includedir)/arena/util.h: $(d)/util.h
	$(MKDIR) -p $(@D)
	$(INSTALL_DATA) $^ $> $@

$(DESTDIR)$(includedir)/arena/rbits.h: $(d)/rbits.h
	$(MKDIR) -p $(@D)
	$(INSTALL_DATA) $^ $> $@

$(DESTDIR)$(libdir)/libarena.a: $(d)/libarena.a
	$(MKDIR) -p $(@D)
	$(INSTALL_DATA) $^ $> $@
	@if [ -f  $(^:.a=.la) $(>:.a=.la) ]; then \
		$(LIBTOOL) --mode=install $(INSTALL_DATA) $(^:.a=.la) $(>:.a=.la) $(@:.a=.la); \
	fi

.INTERMEDIATE: -larena

-larena: $(DESTDIR)$(libdir)/libarena.a $(DESTDIR)$(includedir)/arena/arena.h \
         $(DESTDIR)$(includedir)/arena/pool.h $(DESTDIR)$(includedir)/arena/proto.h \
         $(DESTDIR)$(includedir)/arena/rbits.h $(DESTDIR)$(includedir)/arena/util.h \
         $(DESTDIR)$(includedir)/arena/align.h

install: -larena


#
# Uninstall targets
#
.PHONY: $(d)/uninstall

$(d)/uninstall:
	rm -fr $(DESTDIR)$(includedir)/arena
	rm -f $(DESTDIR)$(libdir)/libarena.a
	@$(LIBTOOL) --mode=uninstall rm -f $(DESTDIR)$(libdir)/libarena.la 2>/dev/null || true

uninstall: $(d)/uninstall


#
# Clean targets
#
.PHONY: $(d)/clean

$(d)/clean:
	echo $(OBJS_$(@D))
	rm -f $(@D)/arena $(@D)/pool $(@D)/*.o $(@D)/*.a
	@$(LIBTOOL) --mode=clean rm -f $(OBJS_$(@D):.o=.lo) $(@D)/libarena.la 2>/dev/null || true

clean: $(d)/clean

# Non-recursive epilogue.
include $(Epilogue)

