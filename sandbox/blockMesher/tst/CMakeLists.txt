add_executable(GMDS_BLOCK_MESHER_TEST
        BlockMesherTestSuite.h
        main_test.cpp)
#==============================================================================
target_link_libraries(GMDS_BLOCK_MESHER_TEST PUBLIC
        ${GMDS_LIB}
        GTest::gtest)
#==============================================================================
gtest_discover_tests(GMDS_BLOCK_MESHER_TEST
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
#==============================================================================
