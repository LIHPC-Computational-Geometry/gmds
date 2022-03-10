add_executable(GMDS_MEDUSA_TEST
        MedusaTestSuite.h
        main_test.cpp)

target_link_libraries(GMDS_MEDUSA_TEST PUBLIC
        Medusa
        gtest)

target_include_directories(GMDS_MEDUSA_TEST PRIVATE .)

set_target_properties(GMDS_MEDUSA_TEST PROPERTIES
        CXX_STANDARD 14
        CXX_STANDARD_REQUIRED YES)
#==============================================================================
gtest_discover_tests(GMDS_MEDUSA_TEST
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
#==============================================================================
