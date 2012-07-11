# This function finds all tests embedded in .cc files in the current directory and 
# creates build and test targets for them.
function(find_tests)
  # Create a list of all .cc files in the directory.
  file(GLOB files "*.cc")
  foreach(file ${files}) # For each file...
    # Check to see whether the processed file contains the string 
    # BUILD_TESTS. We only process these files.
    file(READ ${file} file_contents)
    string(REGEX MATCH BUILD_TESTS is_test ${file_contents})

    if (${is_test} MATCHES BUILD_TESTS) # File contains BUILD_TESTS, so it's a test.

      # Preprocess the file with the C preprocessor. This gets rid of all the 
      # crap that isn't relevant to this build. Notice how we must carefully parse 
      # all the arguments to the preprocessor into a ;-delimited list!
      set(cpp_flags "-E -DBUILD_TESTS -DCMAKE_FINDING_TESTS ${CMAKE_CXX_FLAGS}")
      string(REPLACE " " ";" cpp_flags ${cpp_flags})
      execute_process(COMMAND ${CMAKE_CXX_COMPILER} ${cpp_flags} ${file}
                      OUTPUT_VARIABLE pp_file_contents
                      ERROR_VARIABLE pp_errors)

      # Retrieve just the filename with no extension and prepend test_ to it.
      get_filename_component(testname ${file} NAME_WE)
      if (USE_MPI EQUAL 0)
        set(testname test_${testname})
      else ()
        set(testname test_${testname}_MPI)
      endif ()

      # Add the test executable and add compiler flags, libraries, etc.
      add_executable(${testname} ${file})
      set_target_properties(${testname} PROPERTIES COMPILE_FLAGS "-DBUILD_TESTS -DTEST_PASSED=\\\"passed\\\" -DTEST_FAILED=\\\"failed\\\"")

      # Parse test libraries from TEST_LIBRARIES. Use the preprocessed file.
      string(REGEX MATCH "TEST_LIBRARIES\\([a-zA-Z0-9_, ]+\\)" test_libs ${pp_file_contents})
      if (NOT ${test_libs} STREQUAL "")
        # I'm not pretending I completely understand the following line and all the backslashes, but
        # it seems to do okay.
        string(REGEX REPLACE ".*TEST_LIBRARIES(.+)$" "\\1" test_libs ${test_libs})
        string(REPLACE \( "" test_libs ${test_libs})
        string(REPLACE \) "" test_libs ${test_libs})
        string(REPLACE " " "" test_libs ${test_libs})
        string(REPLACE , ";" test_libs ${test_libs})
      endif ()
      if (USE_MPI EQUAL 0)
        target_link_libraries(${testname} ${test_libs} mpi_stubs ${CHARYBDIS_TEST_LIBS})
      else ()
        target_link_libraries(${testname} ${test_libs} ${CHARYBDIS_TEST_LIBS})
      endif ()

      # FIXME: Parse the number of threads from TEST_THREADS

      # Add the test instance(s).
      if (USE_MPI EQUAL 0)
        add_test(run_${testname} ${CMAKE_CURRENT_BINARY_DIR}/${testname})
        # Tests that report "failed" are those that have failed.
        set_tests_properties(run_${testname} PROPERTIES FAIL_REGULAR_EXPRESSION "failed")
      else ()
        # Parse numbers of processes from TEST_MPI_PROCESSES (in the preprocessed file)
        string(REGEX MATCH "TEST_MPI_PROCESSES\\([a-zA-Z0-9 , ]+\\)" test_procs ${pp_file_contents})
        if (NOT ${test_procs} STREQUAL "")
          string(REGEX REPLACE ".*TEST_MPI_PROCESSES(.+)$" "\\1" test_procs ${test_procs})
          string(REPLACE \( "" test_procs ${test_procs})
          string(REPLACE \) "" test_procs ${test_procs})
          string(REPLACE " " "" test_procs ${test_procs})
          string(REPLACE , ";" test_procs ${test_procs})
        else ()
          set(test_procs "1;")
        endif ()

        # Add one test instance for each number of processes.
        foreach(proc ${test_procs}) # For each file...
          add_test(run_${testname}_${proc}_procs ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${proc} ${MPIEXEC_PREFLAGS} ${CMAKE_CURRENT_BINARY_DIR}/${testname} ${MPIEXEC_POSTFLAGS})
          # Tests that report "failed" are those that have failed.
          set_tests_properties(run_${testname}_${proc}_procs PROPERTIES FAIL_REGULAR_EXPRESSION "failed")
        endforeach()  
      endif ()

    endif()

  endforeach()  
endfunction(find_tests)
