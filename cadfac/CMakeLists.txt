#==============================================================================
# LIBRARY DEFINTION (SOURCE FILES)
#==============================================================================
# Nommer tout en GMDS_MODULE_NAME, GMDS_SRC, ... dans les composants
set(GMDS_LIB ${LIB_GMDS_CADFAC})
set(GMDS_LIB_PREFIX gmds/cadfac)

set(GMDS_INC
        ${CMAKE_BINARY_DIR}/exports/${GMDS_LIB}_export.h
        inc/gmds/cadfac/FACCurve.h
        inc/gmds/cadfac/FACManager.h
        inc/gmds/cadfac/FACPoint.h
        inc/gmds/cadfac/FACSurface.h
        inc/gmds/cadfac/FACVolume.h
        )
set(GMDS_SRC
        src/FACCurve.cpp
        src/FACManager.cpp
        src/FACPoint.cpp
        src/FACSurface.cpp
        src/FACVolume.cpp
        )

#==============================================================================
# GTS
find_package(PkgConfig REQUIRED)
pkg_check_modules(PCRE REQUIRED IMPORTED_TARGET libpcre2-8) # required by glib
pkg_search_module(GLIB2 REQUIRED IMPORTED_TARGET glib-2.0) # required by gts
pkg_search_module(GTS REQUIRED IMPORTED_TARGET gts)
#==============================================================================

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
        ${LIB_GMDS_CAD}
        ${LIB_GMDS_IG}
        ${LIB_GMDS_IG_ALGO}
        ${LIB_GMDS_IO}
        ANN
        PRIVATE
        PkgConfig::GTS PkgConfig::GLIB2 PkgConfig::PCRE)

target_compile_features(${GMDS_LIB} PUBLIC cxx_std_17)

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

