#
# Non-recursive Make prologue. Pushes callers $(d) onto a stack, for re-entrant
# Rules.mk processing.
#
sp		:= $(sp).x
dirstack_$(sp)	:= $(d)
d		:= $(dir)


# Our parent directory
$(d)/..		?= $(shell dirname "$(d)")
$(d)/..		!= dirname "$(d)"

