add_executable(GMDS_ELGMORPHING_TEST
        elgmorphingTestSuite.h
        main_test.cpp
        )

target_link_libraries(GMDS_ELGMORPHING_TEST PUBLIC
        ${GMDS_LIB}
        ${LIB_GMDS_IO}
        ${LIB_GMDS_BLOCKING}
        gtest)

target_include_directories(GMDS_ELGMORPHING_TEST PRIVATE .)

set_target_properties(GMDS_ELGMORPHING_TEST PROPERTIES
        CXX_STANDARD 14
        CXX_STANDARD_REQUIRED YES)
#==============================================================================
gtest_discover_tests(GMDS_ELGMORPHING_TEST
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
#==============================================================================
