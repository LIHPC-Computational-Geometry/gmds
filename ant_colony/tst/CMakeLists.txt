include(GoogleTest)
add_executable(GMDS_ANT_COLONY_TEST
	Env_test.h
        main_test.cpp)

target_link_libraries(GMDS_ANT_COLONY_TEST PUBLIC
        ${GMDS_LIB}
        ${LIB_GMDS_IO}
        gtest)

target_include_directories(GMDS_ANT_COLONY_TEST PRIVATE .)

set_target_properties(GMDS_ANT_COLONY_TEST PROPERTIES
        CXX_STANDARD 14
        CXX_STANDARD_REQUIRED YES)
#==============================================================================
gtest_discover_tests(GMDS_ANT_COLONY_TEST
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
#==============================================================================
