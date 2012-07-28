add_library(arbi_poisson poisson.c)
set(ARBI_LIBS arbi_poisson;${ARBI_LIBS})
set(ARBI_LIBS arbi_poisson;${ARBI_LIBS} PARENT_SCOPE)

#add_subdirectory(tests)
