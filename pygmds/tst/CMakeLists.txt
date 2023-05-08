#==============================================================================
find_package(Python COMPONENTS Interpreter REQUIRED)

add_test(NAME test_pygmds
        COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/mesh.py ${TEST_SAMPLES_DIR}
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/pygmds
 )
#==============================================================================
