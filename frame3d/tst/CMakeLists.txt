add_executable(GMDS_FRAME3D_TEST
        Frame3dTestSuite.h
        PointGeneratorTestSuite.h
        main_test.cpp )

target_link_libraries(GMDS_FRAME3D_TEST PUBLIC
        ${GMDS_LIB}
        gtest)

target_include_directories(GMDS_FRAME3D_TEST PRIVATE .)

set_target_properties(GMDS_FRAME3D_TEST PROPERTIES
        CXX_STANDARD 14
        CXX_STANDARD_REQUIRED YES)
#==============================================================================
gtest_discover_tests(GMDS_FRAME3D_TEST
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
#==============================================================================
