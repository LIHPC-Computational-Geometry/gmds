add_executable(GMDS_CAD_TEST
        FACManagerTestSuite.h
        GeomLinkerTestSuite.h
        GeomTopologyTestSuite.h
        main_test.cpp)

target_link_libraries(GMDS_CAD_TEST PUBLIC
        ${GMDS_LIB}
        gtest)

target_include_directories(GMDS_CAD_TEST PRIVATE .)

set_target_properties(GMDS_CAD_TEST PROPERTIES
        CXX_STANDARD 14
        CXX_STANDARD_REQUIRED YES)
#==============================================================================
gtest_discover_tests(GMDS_CAD_TEST
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
#==============================================================================
