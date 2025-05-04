add_executable(GMDS_SING_GRAPH_BUILD_TEST
        SingGraphBuilder2DTest.h
        main_test.cpp)
#==============================================================================
target_link_libraries(GMDS_SING_GRAPH_BUILD_TEST PUBLIC
        ${GMDS_LIB}
        GTest::gtest)
#==============================================================================
gtest_discover_tests(GMDS_SING_GRAPH_BUILD_TEST
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
#==============================================================================

# PRE_LINK event to copy glpk dll
if(WIN32)
add_custom_command(TARGET GMDS_SING_GRAPH_BUILD_TEST PRE_LINK       		 
    COMMAND ${CMAKE_COMMAND} -E copy_if_different 							 
        "../../../external/glpk/src/glpk/w64/glpk_4_65.dll"
        $<TARGET_FILE_DIR:GMDS_SING_GRAPH_BUILD_TEST>)
endif ()