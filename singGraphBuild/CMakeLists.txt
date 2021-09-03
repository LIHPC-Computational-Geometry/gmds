#==============================================================================
# LIBRARY DEFINTION (SOURCE FILES)
#==============================================================================
# Nommer tout en GMDS_MODULE_NAME, GMDS_SRC, ... dans les composants
set(GMDS_LIB ${LIB_GMDS_SINGGRAPHBUILD})
set(GMDS_LIB_PREFIX gmds/singGraphBuild)

set(GMDS_INC
        ${CMAKE_BINARY_DIR}/exports/${GMDS_LIB}_export.h
        inc/gmds/singGraphBuild/Tools.h
        inc/gmds/singGraphBuild/SingularityLine.h
        inc/gmds/singGraphBuild/SingularityPoint.h
        inc/gmds/singGraphBuild/SingularityPatch.h
        inc/gmds/singGraphBuild/SingularityGraph.h
        inc/gmds/singGraphBuild/SingularityGraphBuilder.h
        inc/gmds/singGraphBuild/SingularityGraphBuilder2D.h
        )
set(GMDS_SRC
        src/Tools.cpp
        src/SingularityLine.cpp
        src/SingularityPoint.cpp
        src/SingularityPatch.cpp
        src/SingularityGraph.cpp
        src/SingularityGraphBuilder.cpp
        src/SingularityGraphBuilder2D.cpp
        )

#==============================================================================
add_library(${GMDS_LIB} ${GMDS_INC} ${GMDS_SRC})
#==============================================================================
include(GenerateExportHeader)
generate_export_header(${GMDS_LIB} 
    EXPORT_FILE_NAME ${CMAKE_BINARY_DIR}/exports/${GMDS_LIB}_export.h
    EXPORT_MACRO_NAME ${GMDS_LIB}_API)
#==============================================================================
# TARGET DEFINITION
#==============================================================================
include(GNUInstallDirs)

#LIBRARY TO INSTALL
target_link_libraries(${GMDS_LIB} PUBLIC
        ${LIB_GMDS_IG}
        ${LIB_GMDS_IG_ALGO}
        ${LIB_GMDS_IO}
        ${LIB_GMDS_FRAME}
        glpk)

target_compile_features(${GMDS_LIB} PUBLIC cxx_std_14)

# INCLUDE TO INSTALL
target_include_directories(${GMDS_LIB} PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/inc>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_PREFIX}/include>
        )
set_target_properties(${GMDS_LIB} PROPERTIES PUBLIC_HEADER "${GMDS_INC}")

install(TARGETS ${GMDS_LIB}
        EXPORT GMDS_SUITE
        PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${GMDS_LIB_PREFIX}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

#==============================================================================
if(WITH_TEST)
    add_subdirectory(tst)
endif(WITH_TEST)


#==============================================================================
# EXECUTABLE
#==============================================================================
add_executable(singGraphBuilder
        src/main.cpp)
target_link_libraries(singGraphBuilder PRIVATE ${GMDS_LIB})
install(TARGETS singGraphBuilder
        DESTINATION ${CMAKE_INSTALL_BINDIR})
