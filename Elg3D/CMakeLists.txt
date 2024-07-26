#==============================================================================
# ELG3D CORE MODULE
#==============================================================================
project(elg3d)

add_definitions(-DDLLEXPORT)

#==============================================================================
# KOKKOS
#==============================================================================
find_package(Kokkos REQUIRED)

#==============================================================================
# OPENMP
#==============================================================================
find_package(OpenMP)

if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

#==============================================================================
# EXODUSII
#==============================================================================
if (ELG3D_WITH_EXODUSII)
    find_package(SEACASExodus)
    add_compile_definitions(ELG3D_WITH_EXODUSII)
    # include directory for exodusii and netCDF;
    # it is currently not set by the exodus target when installed using spack
    include_directories(${SEACASExodus_INCLUDE_DIRS})
    include_directories(${SEACASExodus_TPL_INCLUDE_DIRS})
endif()

#==============================================================================
# GTS
#==============================================================================
find_package(PkgConfig REQUIRED)
pkg_search_module(GLIB2 REQUIRED IMPORTED_TARGET glib-2.0)
pkg_search_module(GTS REQUIRED IMPORTED_TARGET gts)

#==============================================================================
# INTERSECT
#==============================================================================
#==============================================================================
# GLPK LIBRARY
#==============================================================================

include_directories(${GLPK_INC})
link_directories(${GLPK_LIB})
#==============================================================================
# CODE STRUCTURE
#==============================================================================

include_directories("${CMAKE_CURRENT_SOURCE_DIR}/inc")

set(ELG3D_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)

set(ELG3D_SRC_FILES

        inc/ELG3D/ALGOCMPNT/AssignCells.h
        inc/ELG3D/ALGOCMPNT/Cavity.h
        inc/ELG3D/ALGOCMPNT/InitData.h
        inc/ELG3D/ALGOCMPNT/MeshExtractor.h
        inc/ELG3D/DATACMPNT/DummyClassElg3D.h
        inc/ELG3D/DATACMPNT/FacetedCurveGeomServices.h
        inc/ELG3D/DATACMPNT/FacetedSurfaceGeomServices.h
        inc/ELG3D/DATACMPNT/FracPres.h
        inc/ELG3D/DATACMPNT/Parameters.h
        inc/ELG3D/ALGOCMPNT/BadPillowDetection.h
        inc/ELG3D/ALGOCMPNT/ExtractGeomModel.h
        inc/ELG3D/ALGOCMPNT/FracPresEnforcement.h
        inc/ELG3D/ALGOCMPNT/InterfaceNodesPos.h
        inc/ELG3D/ALGOCMPNT/InterfaceNodesPosSmoothVF.h
        inc/ELG3D/ALGOCMPNT/ManifoldDetection.h
        inc/ELG3D/DATACMPNT/MaterialAssignment.h
        inc/ELG3D/ALGOCMPNT/MaterialInterfaces.h
        inc/ELG3D/ALGOCMPNT/OptimizationSmooth.h
        inc/ELG3D/ALGOCMPNT/PixelsRepartitioning.h
        inc/ELG3D/ALGOCMPNT/r2d.h
        inc/ELG3D/ALGOCMPNT/r3d.h
        inc/ELG3D/ALGOCMPNT/Refinement.h
        inc/ELG3D/ALGOCMPNT/SubsetProblem.h
        inc/ELG3D/ALGOCMPNT/Tools.h
        src/ALGOCMPNT/AssignCells.cpp
        src/ALGOCMPNT/Cavity.cpp
        src/ALGOCMPNT/InitData.cpp
        src/ALGOCMPNT/MeshExtractor.cpp
        src/DATACMPNT/DummyClassElg3D.cpp
        src/DATACMPNT/FacetedCurveGeomServices.cpp
        src/DATACMPNT/FacetedSurfaceGeomServices.cpp
        src/DATACMPNT/FracPres.cpp
        src/ALGOCMPNT/BadPillowDetection.cpp
        src/ALGOCMPNT/ExtractGeomModel.cpp
        src/ALGOCMPNT/FracPresEnforcement.cpp
        src/ALGOCMPNT/InterfaceNodesPos.cpp
        src/ALGOCMPNT/InterfaceNodesPosSmoothVF.cpp
        src/ALGOCMPNT/ManifoldDetection.cpp
        src/ALGOCMPNT/OptimizationSmooth.cpp
        src/ALGOCMPNT/r2d.c
        src/ALGOCMPNT/r3d.c
        src/ALGOCMPNT/Refinement.cpp
        src/ALGOCMPNT/SubsetProblem.cpp
        src/DATACMPNT/MaterialAssignment.cpp
        src/DATACMPNT/Parameters.cpp
        src/ALGOCMPNT/MaterialInterfaces.cpp
        src/ALGOCMPNT/PixelsRepartitioning.cpp
        inc/ELG3D/ALGOCMPNT/BoundingBoxGeomAssociation.h src/ALGOCMPNT/BoundingBoxGeomAssociation.cpp inc/ELG3D/ALGOCMPNT/MaterialGradientComputation.h src/ALGOCMPNT/MaterialGradientComputation.cpp src/ALGOCMPNT/Tools.cpp inc/ELG3D/ALGOCMPNT/SmartLaplacian.h src/ALGOCMPNT/SmartLaplacian.cpp inc/ELG3D/ALGOCMPNT/Pillow.h src/ALGOCMPNT/Pillow.cpp inc/ELG3D/ALGOCMPNT/MoveToNewPos.h src/ALGOCMPNT/MoveToNewPos.cpp)


#==============================================================================
# TEST MODULE
#==============================================================================

set(ELG3D_TEST ON CACHE BOOL "Elg3D multithread test module")
if (ELG3D_TEST)
  add_subdirectory(test)
endif ()

#==============================================================================
# LIBRARY BUILDING
#==============================================================================

add_executable(run_dummy src/main_dummy.cpp)
add_executable(run_dummy2 src/main_dummy2.cpp)
add_executable(run_elg2d src/main_elg2d.cpp)
add_executable(run_elg3d src/main_elg3d.cpp)

add_executable(run_pixelVF2d src/main_pixelVF2d.cpp)
add_executable(run_pixelVF3d src/main_pixelVF3d.cpp)

add_executable(run_exodusreader2d src/main_exodusreader2d.cpp)
add_executable(run_exodus_pixelVF2d src/main_exodus_pixelVF2d.cpp)
add_executable(run_exodus_compare2d src/main_exodus_compare2d.cpp)



add_executable(run_exodusreader3d src/main_exodusreader3d.cpp)
add_executable(run_exodus_pixelVF3d src/main_exodus_pixelVF3d.cpp)
add_executable(run_exodus_compare3d src/main_exodus_compare3d.cpp)

add_executable(run_elg2d_imr2019 src/main_elg2d_imr2019.cpp)
add_executable(run_elg3d_imr2019 src/main_elg3d_imr2019.cpp)


add_executable(run_compare3d src/main_compare_3d.cpp)

add_library(elg3d ${ELG3D_SRC_FILES})

target_link_libraries(run_dummy elg3d kmds kokkos dl)
target_link_libraries(run_dummy2 elg3d kmds kokkos dl)
target_link_libraries(run_elg2d elg3d kmds kokkos dl)
target_link_libraries(run_elg3d elg3d kmds kokkos dl)
target_link_libraries(run_pixelVF2d elg3d kmds Kokkos::kokkoscore)
target_link_libraries(run_pixelVF3d elg3d kmds Kokkos::kokkoscore)

target_link_libraries(run_exodusreader2d elg3d kmds Kokkos::kokkoscore)
target_link_libraries(run_exodus_pixelVF2d elg3d kmds Kokkos::kokkoscore)
target_link_libraries(run_exodus_compare2d elg3d kmds Kokkos::kokkoscore)
target_link_libraries(run_exodusreader3d elg3d kmds Kokkos::kokkoscore)
target_link_libraries(run_exodus_pixelVF3d elg3d kmds Kokkos::kokkoscore)
target_link_libraries(run_exodus_compare3d elg3d kmds Kokkos::kokkoscore)

target_link_libraries(run_elg3d_imr2019 elg3d kmds Kokkos::kokkoscore)

target_link_libraries(run_elg2d_imr2019 elg3d kmds Kokkos::kokkoscore)

target_link_libraries(run_compare3d elg3d kmds Kokkos::kokkoscore)

target_link_libraries(elg3d
#	PUBLIC ${LIB_GMDS_MATH} ${LIB_GMDS_IG} ${LIB_GMDS_CAD} ${LIB_GMDS_UTILS} kmds Kokkos::kokkoscore Kokkos::kokkoscontainers GCO  gts glib-2.0
        PUBLIC ${LIB_GMDS_MATH} ${LIB_GMDS_IG} ${LIB_GMDS_CAD}
        ${LIB_GMDS_CADFAC} ${LIB_GMDS_UTILS} kmds Kokkos::kokkoscore Kokkos::kokkoscontainers GCO PkgConfig::GTS PkgConfig::GLIB2
	INTERFACE glpk)

if(ELG3D_WITH_EXODUSII)
    target_link_libraries(elg3d PUBLIC exodus)
endif ()

#add_custom_command(
#   TARGET elg3d
#   POST_BUILD
#   COMMAND ${CMAKE_COMMAND} ARGS -E copy  $<TARGET_FILE:elg3d> ${CMAKE_BINARY_DIR}/gen/
#   COMMENT "Copy elg3d lib"
#   VERBATIM
#   )

#target_include_directories(gmds INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>/generated>)



