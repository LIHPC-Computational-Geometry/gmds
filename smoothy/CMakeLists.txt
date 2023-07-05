#==============================================================================
# LIBRARY DEFINTION (SOURCE FILES)
#==============================================================================
# Nommer tout en GMDS_MODULE_NAME, GMDS_SRC, ... dans les composants
set(GMDS_LIB ${LIB_GMDS_SMOOTHY})
set(GMDS_LIB_PREFIX gmds/smoothy)

set(GMDS_INC
        ${CMAKE_BINARY_DIR}/exports/${GMDS_LIB}_export.h
        inc/gmds/smoothy/AbstractSmoother.h
        inc/gmds/smoothy/AngleBasedQuadSmoother.h
        inc/gmds/smoothy/LaplacianSmoother.h
        inc/gmds/smoothy/HLBFGSWrapper.h
        inc/gmds/smoothy/EllipticSmoothing.h
        inc/gmds/smoothy/EllipticSmoother2D.h
        )
set(GMDS_SRC
        src/AbstractSmoother.cpp
        src/AngleBasedQuadSmoother.cpp
        src/LaplacianSmoother.cpp
        src/HLBFGSWrapper.cpp
        src/EllipticSmoothing.cpp
        src/EllipticSmoother2D.cpp
        )
#==============================================================================
# HLBFGS include link (old fashion include)
include_directories(${CMAKE_SOURCE_DIR}/external/HLBFGS)
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
        ${LIB_GMDS_CAD}
        ${LIB_GMDS_CADFAC}
        ${LIB_GMDS_IG_ALGO}
        ${LIB_GMDS_IO}
        PRIVATE hlbfgs)
#HLBGFS car remain private as it is forced as STATIC in the HLBGS directory

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
# EXECUTABLES
#==============================================================================
add_executable(smoothy src/smoothy.cpp)
target_link_libraries(smoothy PRIVATE ${GMDS_LIB})
target_compile_features(smoothy PUBLIC cxx_std_14)

add_executable(invertibleMapper2D src/invertibleMapper2D.cpp)
target_link_libraries(invertibleMapper2D PRIVATE ${GMDS_LIB})
target_compile_features(invertibleMapper2D PUBLIC cxx_std_14)


add_executable(morph2d src/morph2d.cpp)
target_link_libraries(morph2d PRIVATE ${GMDS_LIB})
target_compile_features(morph2d PUBLIC cxx_std_14)

