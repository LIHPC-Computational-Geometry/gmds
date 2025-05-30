#==============================================================================
# LIBRARY DEFINITION (SOURCE FILES)
#==============================================================================
# Explicitly used the name given in this preamble
set(GMDS_LIB ${LIB_GMDS_QUADFRONT})
set(GMDS_LIB_PREFIX gmds/quadfront)

set(GMDS_INC
        ${CMAKE_BINARY_DIR}/exports/${GMDS_LIB}_export.h
        inc/gmds/quadfront/Quadfront.h
)
set(GMDS_SRC
        src/Quadfront.cpp
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
        ${LIB_GMDS_CAD}
        ${LIB_GMDS_IG})

#==============================================================================
# NOTHING TO UPDATE BELOW
#==============================================================================

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
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/gmds
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/gmds)

#==============================================================================
if(WITH_TEST)
        add_subdirectory(tst)
endif(WITH_TEST)
#==============================================================================
# EXECUTABLE
#==============================================================================
add_executable(quadfront_exec src/main.cpp)
target_link_libraries(quadfront_exec PRIVATE ${GMDS_LIB})
target_compile_features(quadfront_exec PUBLIC cxx_std_14)
install(TARGETS quadfront_exec)