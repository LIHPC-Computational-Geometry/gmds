add_executable(GMDS_QUADFRONT_TEST
        DummyTestSuite.h
        main_test.cpp
)
#==============================================================================
target_link_libraries(GMDS_QUADFRONT_TEST PUBLIC
        ${GMDS_LIB}
        GTest::gtest)
#==============================================================================
gtest_discover_tests(GMDS_QUADFRONT_TEST
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
#==============================================================================
