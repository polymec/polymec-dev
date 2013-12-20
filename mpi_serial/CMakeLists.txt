add_library(mpi_serial mpi_serial.c)
set(POLYMEC_LIBS ${POLYMEC_LIBS};mpi_serial PARENT_SCOPE)
set(POLYMEC_LIBS ${POLYMEC_LIBS};mpi_serial)

install(TARGETS mpi_serial DESTINATION lib)
install(FILES mpi.h DESTINATION include/polymec)
