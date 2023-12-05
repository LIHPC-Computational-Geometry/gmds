#==============================================================================
# LIBRARY DEFINTION (SOURCE FILES)
#==============================================================================
# Nommer tout en GMDS_MODULE_NAME, GMDS_SRC, ... dans les composants
set(GMDS_LIB ${LIB_GMDS_MATH})
set(GMDS_LIB_PREFIX gmds/math)
set(GMDS_INC
        ${CMAKE_BINARY_DIR}/exports/${GMDS_LIB}_export.h
        inc/gmds/math/AxisAngleRotation.h
        inc/gmds/math/BezierCurve.h
        inc/gmds/math/BezierTriangle.h
        inc/gmds/math/Chart.h
        inc/gmds/math/Constants.h
        inc/gmds/math/Cross.h
        inc/gmds/math/Cross2D.h
        inc/gmds/math/DiscretizationScheme1D.h
        inc/gmds/math/FE.h
        inc/gmds/math/Hexahedron.h
        inc/gmds/math/Line.h
        inc/gmds/math/Matrix.h
        inc/gmds/math/Numerics.h
        inc/gmds/math/Orientation.h
        inc/gmds/math/Plane.h
        inc/gmds/math/Prism3.h
        inc/gmds/math/Point.h
        inc/gmds/math/Pyramid.h
        inc/gmds/math/Quadrilateral.h
        inc/gmds/math/QualityMeasure.h
        inc/gmds/math/Quaternion.h
        inc/gmds/math/Ray.h
        inc/gmds/math/Segment.h
        inc/gmds/math/SHarmonicL4.h
        inc/gmds/math/SHFrame.h
        inc/gmds/math/Tetrahedron.h
        inc/gmds/math/Triangle.h
        inc/gmds/math/TransfiniteInterpolation.h
        inc/gmds/math/Vector.h
        inc/gmds/math/VectorDyn.h
        inc/gmds/math/StreamComputation.h
        )
set(GMDS_SRC
        src/AxisAngleRotation.cpp
        src/BezierCurve.cpp
        src/BezierTriangle.cpp
        src/Chart.cpp
        src/Cross.cpp
        src/Cross2D.cpp
        src/DiscretizationScheme1D.cpp
        src/FE.cpp
        src/Hexahedron.cpp
        src/Line.cpp
        src/Numerics.cpp
        src/Orientation.cpp
        src/Plane.cpp
        src/Point.cpp
        src/Prism3.cpp
        src/Pyramid.cpp
        src/Quadrilateral.cpp
        src/QualityMeasure.cpp
        src/Quaternion.cpp
        src/Ray.cpp
        src/Segment.cpp
        src/SHarmonicL4.cpp
        src/Tetrahedron.cpp
        src/Triangle.cpp
        src/StreamComputation.cpp
        src/VectorDyn.cpp
        src/Orientation.cpp
        src/Vector.cpp
        src/TransfiniteInterpolation.cpp inc/gmds/math/BezierSurface.h src/BezierSurface.cpp inc/gmds/math/BezierHex.h src/BezierHex.cpp)
#==============================================================================
find_package(Eigen3 3.3 REQUIRED)
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
target_link_libraries(${GMDS_LIB}
        PUBLIC
        ${LIB_GMDS_UTILS}
        Predicates_psm
        Eigen3::Eigen)

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

