#==============================================================================
# LIBRARY DEFINTION (SOURCE FILES)
#==============================================================================
# Nommer tout en GMDS_MODULE_NAME, GMDS_SRC, ... dans les composants
set(GMDS_LIB ${LIB_GMDS_IG})
set(GMDS_LIB_PREFIX gmds/ig)

set(GMDS_INC
        ${CMAKE_BINARY_DIR}/exports/${GMDS_LIB}_export.h
        inc/gmds/ig/Cell.h
        inc/gmds/ig/CellGroup.h
        inc/gmds/ig/Edge.h
        inc/gmds/ig/EdgeContainer.h
        inc/gmds/ig/Face.h
        inc/gmds/ig/FaceContainer.h
        inc/gmds/ig/Mesh.h
        inc/gmds/ig/Node.h
        inc/gmds/ig/NodeContainer.h
        inc/gmds/ig/Region.h
        inc/gmds/ig/RegionContainer.h
        inc/gmds/ig/MeshDoctor.h)

set(GMDS_SRC
        src/Cell.cpp
        src/Edge.cpp
        src/EdgeContainer.cpp
        src/Face.cpp
        src/FaceContainer.cpp
        src/Mesh.cpp
        src/MeshDoctor.cpp
        src/Node.cpp
        src/NodeContainer.cpp
        src/Region.cpp
        src/RegionContainer.cpp
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
target_link_libraries(${GMDS_LIB} PUBLIC ${LIB_GMDS_MATH})
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
