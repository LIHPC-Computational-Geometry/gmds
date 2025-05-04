add_executable(GMDS_PADDING_TEST
        PaddingTestSuite.h
        main_test.cpp)
#==============================================================================
target_link_libraries(GMDS_PADDING_TEST PUBLIC
        ${GMDS_LIB}
        ${LIB_GMDS_IG_ALGO}
        GTest::gtest)
#==============================================================================
gtest_discover_tests(GMDS_PADDING_TEST
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
#==============================================================================
