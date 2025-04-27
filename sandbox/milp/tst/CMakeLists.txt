add_executable(GMDS_MILP_TEST
        DummyTestSuite.h
        main_test.cpp
)
#==============================================================================
target_link_libraries(GMDS_MILP_TEST PUBLIC
        ${GMDS_LIB}
        GTest::gtest)
#==============================================================================
gtest_discover_tests(GMDS_MILP_TEST
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
#==============================================================================
