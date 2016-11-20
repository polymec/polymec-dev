add_library(mpi_serial mpi_serial.c)
set_target_properties(mpi_serial PROPERTIES FOLDER Libraries)
set(POLYMEC_LIBRARIES ${POLYMEC_LIBRARIES};mpi_serial PARENT_SCOPE)
set(POLYMEC_LIBRARIES ${POLYMEC_LIBRARIES};mpi_serial)

install(TARGETS mpi_serial DESTINATION lib)
install(FILES mpi.h DESTINATION include/polymec)
