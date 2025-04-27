add_executable(GMDS_MEDIALAXIS_TEST
        DummyTestSuite.h
        main_test.cpp
)
#==============================================================================
target_link_libraries(GMDS_MEDIALAXIS_TEST PUBLIC
        ${GMDS_LIB}
        GTest::gtest)
#==============================================================================
gtest_discover_tests(GMDS_MEDIALAXIS_TEST
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
#==============================================================================
