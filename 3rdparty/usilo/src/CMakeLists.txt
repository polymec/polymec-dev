include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_BINARY_DIR})
add_library(usilo H5FDsilo.c alloc.c extface.c silo.c silo_debug.c 
            silo_hdf5.c silo_unknown.c silo_ns.c silo_version.c)
target_link_libraries(usilo ${HDF5_HL_LIB_NAME};${HDF5_LIB_NAME};dl;z)

