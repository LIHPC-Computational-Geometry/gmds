#==============================================================================
find_package(Python COMPONENTS Interpreter REQUIRED)

add_test(NAME test_pygmds
        COMMAND pytest ${CMAKE_CURRENT_SOURCE_DIR} ${TEST_SAMPLES_DIR} -v
 )
set_tests_properties(test_pygmds
        PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/pygmds:$ENV{PYTHONPATH}")
#==============================================================================
