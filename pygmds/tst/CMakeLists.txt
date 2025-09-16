#==============================================================================
message("ENVIRONMENT ENVPYTHONPATH $ENV{PYTHONPATH}")

add_test(NAME test_pygmds_mesh
        COMMAND pytest ${CMAKE_CURRENT_SOURCE_DIR}/test_mesh.py ${TEST_SAMPLES_DIR} -v
 )
add_test(NAME test_pygmds_geometry
        COMMAND pytest ${CMAKE_CURRENT_SOURCE_DIR}/test_geometry.py ${TEST_SAMPLES_DIR} -v
 )
if(ENABLE_BLOCKING)
	add_test(NAME test_pygmds_blocking
        COMMAND pytest ${CMAKE_CURRENT_SOURCE_DIR}/test_blocking.py ${TEST_SAMPLES_DIR} -v
 )
endif ()

set_tests_properties(test_pygmds_mesh
        PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/pygmds:$ENV{PYTHONPATH}")
set_tests_properties(test_pygmds_geometry
        PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/pygmds:$ENV{PYTHONPATH}")
if(ENABLE_BLOCKING)
	set_tests_properties(test_pygmds_blocking
        	PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/pygmds:$ENV{PYTHONPATH}")
endif ()
#==============================================================================
