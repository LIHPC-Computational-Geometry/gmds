add_executable(GMDS_SERVICES_TEST
        FieldDataPropertyTestSuite.h
        IdleService.h
        MeshDataPropertyTestSuite.h
        ScalarDataPropertyTestSuite.h
        main_test.cpp)

target_link_libraries(GMDS_SERVICES_TEST PUBLIC
        ${GMDS_LIB}
        gtest)

target_include_directories(GMDS_SERVICES_TEST PRIVATE .)

set_target_properties(GMDS_SERVICES_TEST PROPERTIES
        CXX_STANDARD 14
        CXX_STANDARD_REQUIRED YES)
#==============================================================================
gtest_discover_tests(GMDS_SERVICES_TEST
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
#==============================================================================
