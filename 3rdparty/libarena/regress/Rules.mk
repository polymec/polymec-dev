# Non-recursive prologue.
include $(Prologue)


CFLAGS_$(d)	:= $(CFLAGS)
CPPFLAGS_$(d)	:= -I$($(d)/..)/src -DLIBARENA_SOURCE $(CPPFLAGS)
LDFLAGS_$(d)	:= -L$($(d)/..)/src $(LDFLAGS)

$(d)/stacked: $(d)/stacked.c $($(d)/..)/src/libarena.a
	$(CC) $(CFLAGS_$(@D)) $(CPPFLAGS_$(@D)) $(LDFLAGS_$(@D)) -o $@ $^ $> -lz

$(d)/grow: $(d)/grow.c $($(d)/..)/src/libarena.a
	$(CC) $(CFLAGS_$(@D)) $(CPPFLAGS_$(@D)) $(LDFLAGS_$(@D)) -o $@ $^ $>

.PHONY: $(d)/check

$(d)/check: $(d)/stacked $(d)/grow
	$(@D)/stacked

check: $(d)/check



.PHONY: $(d)/clean

$(d)/clean:
	rm -f $(@D)/stacked

clean: $(d)/clean


# Non-recursive epilogue.
include $(Epilogue)

