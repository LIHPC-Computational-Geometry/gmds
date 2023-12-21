add_executable(GMDS_BLOCKINGQUALITY_TEST
        BlockQualityTestSuite.h
        MCTSTestSuite.h
        main_test.cpp)

target_link_libraries(GMDS_BLOCKINGQUALITY_TEST PUBLIC
        ${GMDS_LIB}
        ${LIB_GMDS_IG_ALGO}
        GTest::gtest)

target_include_directories(GMDS_BLOCKINGQUALITY_TEST PRIVATE .)

set_target_properties(GMDS_BLOCKINGQUALITY_TEST PROPERTIES
        CXX_STANDARD 14
        CXX_STANDARD_REQUIRED YES)
#==============================================================================
gtest_discover_tests(GMDS_BLOCKINGQUALITY_TEST
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
#====================================