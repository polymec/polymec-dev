#
# Non-recursive Make epilogue. Restores callers $(d), for re-entrant Rules.mk
# processing.
#
d	:= $(dirstack_$(sp))

sp	?= $(shell basename "$(sp)")
sp	!= basename "$(sp)"

