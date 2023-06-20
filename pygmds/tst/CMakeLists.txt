#==============================================================================
find_package(Python COMPONENTS Interpreter REQUIRED)

add_test(NAME test_pygmds
        COMMAND ${Python_EXECUTABLE} -m pytest ${CMAKE_CURRENT_SOURCE_DIR} ${TEST_SAMPLES_DIR}
 )
set_tests_properties(test_pygmds
        PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/pygmds:$ENV{PYTHONPATH}")
#==============================================================================
