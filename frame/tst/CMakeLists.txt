add_executable(GMDS_FRAME_TEST
        Frame2dTestSuite.h
        main_test.cpp)

target_link_libraries(GMDS_FRAME_TEST PUBLIC
        ${GMDS_LIB}
        gtest)

target_include_directories(GMDS_FRAME_TEST PRIVATE .)

set_target_properties(GMDS_FRAME_TEST PROPERTIES
        CXX_STANDARD 14
        CXX_STANDARD_REQUIRED YES)
#==============================================================================
gtest_discover_tests(GMDS_FRAME_TEST
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
#==============================================================================
