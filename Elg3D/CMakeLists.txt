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
#include(FindPkgConfig)
find_package(PkgConfig REQUIRED)
#pkg_check_modules(glib-2.0 REQUIRED glib-2.0)
#find_path(GLIB_INCLUDE_DIR glib.h PATH_SUFFIXES glib-2.0)
#include_directories(${GLIB_INCLUDE_DIR})
#find_path(GLIBCONFIG_INCLUDE_DIR glibconfig.h PATH_SUFFIXES ../lib/glib-2.0/include)
#include_directories(${GLIBCONFIG_INCLUDE_DIR})

#pkg_check_modules(gts REQUIRED gts)
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GTS_Cflags}")
#find_path(GTS_INCLUDE_DIR gts.h PATH_SUFFIXES include)
#include_directories(${GTS_INCLUDE_DIR})
#find_path(GTS_LIBRARY_DIR libgts.so PATH_SUFFIXES lib)
#link_directories(${GTS_LIBRARY_DIR})

pkg_search_module(PKG_GTS REQUIRED IMPORTED_TARGET gts)


#if(GTS_LIB_PATH)
#    include_directories(${GTS_INC_PATH})
#    link_directories(${GTS_LIB_PATH})
#    add_compile_definitions(ELG3D_WITH_GTS)
#
#    include_directories(/home/legoffn/travail/utils/spack_20200612/opt/spack/linux-ubuntu16.04-broadwell/gcc-5.4.0/glib-2.64.3-e2s7xjrjzvokp7qjqhalgs5b3ep3cowf/include/glib-2.0)
#    include_directories(/home/legoffn/travail/utils/spack_20200612/opt/spack/linux-ubuntu16.04-broadwell/gcc-5.4.0/glib-2.64.3-e2s7xjrjzvokp7qjqhalgs5b3ep3cowf/lib/glib-2.0/include)
#    link_directories(/home/legoffn/travail/utils/spack_20200612/opt/spack/linux-ubuntu16.04-broadwell/gcc-5.4.0/glib-2.64.3-e2s7xjrjzvokp7qjqhalgs5b3ep3cowf/lib)
#endif()

#==============================================================================
# GLIB
#==============================================================================
#include_directories(/usr/include/glib-2.0 /usr/lib64/glib-2.0/include)

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
        PUBLIC ${LIB_GMDS_MATH} ${LIB_GMDS_IG} ${LIB_GMDS_CAD} ${LIB_GMDS_UTILS} kmds Kokkos::kokkoscore Kokkos::kokkoscontainers GCO PkgConfig::PKG_GTS
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



