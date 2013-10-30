#
# Standard GNU Make non-recursive prologue. Pushes callers $(d) onto a stack,
# for re-entrant Rules.mk processing.
#
sp		:= $(sp).x
dirstack_$(sp)	:= $(d)
d		:= $(dir)


#
# Make sure we have our defaults.
#
include $(d)/mk/Variables.mk


#
# Include source rules
#
dir	:= $(d)/src
include $(dir)/Rules.mk


#
# Include regress rules
#
dir	:= $(d)/regress
include $(dir)/Rules.mk


#
# Standard GNU Make non-recursive epilogue. Restores callers $(d), for
# re-entrant Rules.mk processing.
#
d	:= $(dirstack_$(sp))
sp	:= $(basename $(sp))

