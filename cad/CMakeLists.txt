#==============================================================================
# LIBRARY DEFINTION (SOURCE FILES)
#==============================================================================
# Nommer tout en GMDS_MODULE_NAME, GMDS_SRC, ... dans les composants
set(GMDS_LIB ${LIB_GMDS_CAD})
set(GMDS_LIB_PREFIX gmds/cad)

set(GMDS_INC
        ${CMAKE_BINARY_DIR}/exports/${GMDS_LIB}_export.h
        inc/gmds/cad/FACCurve.h
        inc/gmds/cad/FACManager.h
        inc/gmds/cad/FACPoint.h
        inc/gmds/cad/FACSurface.h
        inc/gmds/cad/FACVolume.h
        inc/gmds/cad/GeomCurve.h
        inc/gmds/cad/GeomEntity.h
        inc/gmds/cad/GeomManager.h
        inc/gmds/cad/GeomPoint.h
        inc/gmds/cad/GeomSurface.h
        inc/gmds/cad/GeomMeshLinker.h
        inc/gmds/cad/GeomVolume.h
        inc/gmds/cad/GeomSmoother.h
        )
set(GMDS_SRC
        src/FACCurve.cpp
        src/FACManager.cpp
        src/FACPoint.cpp
        src/FACSurface.cpp
        src/FACVolume.cpp
        src/GeomSurface.cpp
        src/GeomSmoother.cpp
        src/GeomMeshLinker.cpp
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
        ${LIB_GMDS_IO})

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

