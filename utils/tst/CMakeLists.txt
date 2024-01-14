add_executable(GMDS_UTILS_TEST
        UtilsTestSuite.h
        ParamTestSuite.h
        ArrayTestSuite.h
        RandomGeneratorTest.h
        main_test.cpp)
#==============================================================================
target_link_libraries(GMDS_UTILS_TEST PUBLIC
        ${GMDS_LIB}
        GTest::gtest)
#==============================================================================
gtest_discover_tests(GMDS_UTILS_TEST
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
#==============================================================================
