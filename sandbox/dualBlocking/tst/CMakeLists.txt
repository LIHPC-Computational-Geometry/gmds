add_executable(GMDS_DUBLO_TEST
        DuBloTestSuite.h
        main_test.cpp)
#==============================================================================
target_link_libraries(GMDS_DUBLO_TEST PUBLIC
        ${GMDS_LIB}
        GTest::gtest)
#==============================================================================
gtest_discover_tests(GMDS_DUBLO_TEST
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
#==============================================================================
