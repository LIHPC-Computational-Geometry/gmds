#==============================================================================
# KMDS CORE MODULE
#==============================================================================
project(kmds)

add_definitions(-DDLLEXPORT)

#==============================================================================
# OPENMP
#==============================================================================
find_package(OpenMP)

if(OPENMP_FOUND)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()
#==============================================================================
# CODE STRUCTURE
#==============================================================================
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/inc")

#find_package(PkgConfig)
#pkg_check_modules(KOKKOS IMPORTED_TARGET kokkos)
#target_link_libraries(<lib> PkgConfig::KOKKOS)


#find_package(kokkos REQUIRED)
#include_directories(${KOKKOS_INC_PATH})
#link_directories(${KOKKOS_LIB_PATH})
#link_directories(/home/legoffn/travail/utils/spack_20200612/opt/spack/linux-ubuntu16.04-broadwell/gcc-5.4.0/hwloc-2.2.0-mlb3ymifx2ma5a4vygodmqwq3lznr3dk/lib)

find_package(Kokkos REQUIRED)


set(GMDS_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src/DS ${CMAKE_CURRENT_SOURCE_DIR}/src/Utils ${CMAKE_CURRENT_SOURCE_DIR}/src/IO)

set(KMDS_SRC_FILES 
  inc/KM/DS/CellHandle.h
  inc/KM/DS/Connectivity.h
  inc/KM/DS/EContainer.h
  inc/KM/DS/FContainer.h
  inc/KM/DS/Mesh.h
  inc/KM/DS/NContainer.h
  inc/KM/DS/RContainer.h
  inc/KM/IO/VTKWriter.h
  inc/KM/IO/VTKWriter_def.h
  inc/KM/Utils/VariableManager.h
  inc/KM/Utils/Exception.h
  inc/KM/Utils/GrowingView.h
  inc/KM/Utils/FakeTypes.h
  inc/KM/Utils/KFlags.h
  inc/KM/Utils/KTypes.h
  inc/KM/Utils/InitTools.h
  inc/KM/Utils/Variable.h
  inc/KM/Utils/Variable_def.h
  inc/KM/Utils/VariableManager.h
  src/DS/CellHandle.cpp
  src/DS/Connectivity.cpp
  src/DS/EContainer.cpp
  src/DS/FContainer.cpp
  src/DS/Mesh.cpp
  src/DS/NContainer.cpp
  src/DS/RContainer.cpp
  src/Utils/InitTools.cpp
  src/Utils/Variable.cpp
  src/Utils/VariableManager.cpp
        inc/KM/Utils/Graph.h src/Utils/Graph.cpp)

#==============================================================================
# TEST MODULE
#==============================================================================
if(WITH_TEST)
  add_subdirectory(test)
endif(WITH_TEST)

#==============================================================================
# LIBRARY BUILDING
#==============================================================================

#==============================================================================
# TARGET DEFINITION
#==============================================================================
include(GNUInstallDirs)

add_executable(run_kmds_quads src/main.cpp)
add_executable(run_kmds_hex src/main_hex.cpp)
add_library(kmds ${KMDS_SRC_FILES})

target_compile_features(kmds PUBLIC cxx_std_14)

#target_link_libraries(run_kmds_quads kmds kokkos dl)
#target_link_libraries(run_kmds_hex kmds kokkos dl)
#target_link_libraries(run_kmds_quads kmds kokkos)
#target_link_libraries(run_kmds_hex kmds kokkos)
target_link_libraries(run_kmds_quads kmds Kokkos::kokkoscore)
target_link_libraries(run_kmds_hex kmds Kokkos::kokkoscore)
#target_link_libraries(kmds ${LIB_GMDS_MATH} ${LIB_GMDS_IG} kokkos dl hwloc)
target_link_libraries(kmds ${LIB_GMDS_MATH} ${LIB_GMDS_IG} Kokkos::kokkoscore Kokkos::kokkoscontainers)

#add_custom_command(
#   TARGET kmds
#   POST_BUILD
#   COMMAND ${CMAKE_COMMAND} ARGS -E copy  $<TARGET_FILE:kmds> ${CMAKE_BINARY_DIR}/gen/
#   COMMENT "Copy kmds lib"
#   VERBATIM
#   )


# INCLUDE TO INSTALL
target_include_directories(kmds PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/inc>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_PREFIX}/include>
        )

