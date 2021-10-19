add_executable(GMDS_BLOCK_MESHER_TEST
        BlockMesherTestSuite.h
        main_test.cpp )

target_link_libraries(GMDS_BLOCK_MESHER_TEST PUBLIC
        ${GMDS_LIB}
        gtest)

target_include_directories(GMDS_BLOCK_MESHER_TEST PRIVATE .)

set_target_properties(GMDS_BLOCK_MESHER_TEST PROPERTIES
        CXX_STANDARD 14
        CXX_STANDARD_REQUIRED YES)
#==============================================================================
gtest_discover_tests(GMDS_BLOCK_MESHER_TEST
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
#==============================================================================
