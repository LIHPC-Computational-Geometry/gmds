add_executable(GMDS_FRAME_TEST
        Frame2dTestSuite.h
        main_test.cpp)
#==============================================================================
target_link_libraries(GMDS_FRAME_TEST PUBLIC
        ${GMDS_LIB}
        GTest::gtest)
#==============================================================================
gtest_discover_tests(GMDS_FRAME_TEST
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
#==============================================================================
