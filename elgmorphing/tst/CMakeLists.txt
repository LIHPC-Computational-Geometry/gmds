add_executable(GMDS_ELGMORPHING_TEST
        elgmorphingTestSuite.h
        main_test.cpp
        )
#==============================================================================
target_link_libraries(GMDS_ELGMORPHING_TEST PUBLIC
        ${GMDS_LIB}
        ${LIB_GMDS_IO}
        ${LIB_GMDS_BLOCKING}
        ${LIB_GMDS_SMOOTHY}
        GTest::gtest)
#==============================================================================
gtest_discover_tests(GMDS_ELGMORPHING_TEST
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
#==============================================================================
