add_executable(GMDS_SING_GRAPH_BUILD_TEST
        SingGraphBuilder2DTest.h
        main_test.cpp)

target_link_libraries(GMDS_SING_GRAPH_BUILD_TEST PUBLIC
        ${GMDS_LIB}
        gtest)

target_include_directories(GMDS_SING_GRAPH_BUILD_TEST PRIVATE .)

set_target_properties(GMDS_SING_GRAPH_BUILD_TEST PROPERTIES
        CXX_STANDARD 14
        CXX_STANDARD_REQUIRED YES)
#==============================================================================
gtest_discover_tests(GMDS_SING_GRAPH_BUILD_TEST
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
#==============================================================================
