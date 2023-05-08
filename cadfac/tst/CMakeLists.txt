add_executable(GMDS_CADFAC_TEST
        FACManagerTestSuite.h
        GeomLinkerTestSuite.h
        GeomTopologyTestSuite.h
        main_test.cpp)
#==============================================================================
target_link_libraries(GMDS_CADFAC_TEST PUBLIC
        ${GMDS_LIB}
        GTest::gtest)
#==============================================================================
gtest_discover_tests(GMDS_CADFAC_TEST
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
#==============================================================================
