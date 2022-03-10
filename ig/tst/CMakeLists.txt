add_executable(GMDS_IG_TEST
        BooleanMarkTestSuite.h
        MeshTestSuite.h
        MeshDoctorTestSuite.h
        main_test.cpp)

target_link_libraries(GMDS_IG_TEST PUBLIC
        ${GMDS_LIB}
        gtest)

target_include_directories(GMDS_IG_TEST
        PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

#sans doute a virer, transit via le compile_feature du GMDSIG
set_target_properties(GMDS_IG_TEST PROPERTIES
        CXX_STANDARD 14
        CXX_STANDARD_REQUIRED YES)
#==============================================================================
gtest_discover_tests(GMDS_IG_TEST
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
#==============================================================================
