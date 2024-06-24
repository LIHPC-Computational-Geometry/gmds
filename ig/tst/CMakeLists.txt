add_executable(GMDS_IG_TEST
        BooleanMarkTestSuite.h
        Blocking2DTestSuite.h
        EdgeTestSuite.h
        MeshTestSuite.h
        MeshDoctorTestSuite.h
        CellTestSuite.h
        EdgeContainerTestSuite.h
        NodeTestSuite.h
        NodeContainerTestSuite.h
        main_test.cpp)
#==============================================================================
target_link_libraries(GMDS_IG_TEST PUBLIC
        ${GMDS_LIB}
        ${LIB_GMDS_IO}
        GTest::gtest)
#==============================================================================
gtest_discover_tests(GMDS_IG_TEST
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
#==============================================================================
