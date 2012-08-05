add_library(arbi_poisson poisson.c)
set(ARBI_LIBS arbi_poisson;${ARBI_LIBS} PARENT_SCOPE)
set(ARBI_LIBS arbi_poisson;${ARBI_LIBS})

#add_subdirectory(tests)
