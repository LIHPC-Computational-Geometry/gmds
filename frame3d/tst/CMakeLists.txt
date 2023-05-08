add_executable(GMDS_FRAME3D_TEST
        Frame3dTestSuite.h
        PointGeneratorTestSuite.h
        main_test.cpp )
#==============================================================================
target_link_libraries(GMDS_FRAME3D_TEST PUBLIC
        ${GMDS_LIB}
        GTest::gtest)
#==============================================================================
include(GoogleTest)
gtest_discover_tests(GMDS_FRAME3D_TEST
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
#==============================================================================
