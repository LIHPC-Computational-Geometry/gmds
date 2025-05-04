add_executable(GMDS_QUALITY_TEST
        QuadEvalTestSuite.h
        HexEvalTestSuite.h
        main_test.cpp)
#==============================================================================
target_link_libraries(GMDS_QUALITY_TEST PUBLIC
        ${GMDS_LIB}
        GTest::gtest)
#==============================================================================
include(GoogleTest)
gtest_discover_tests(GMDS_QUALITY_TEST
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
#==============================================================================
