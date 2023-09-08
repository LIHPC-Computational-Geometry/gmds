#==============================================================================
# GMDS Project
#==============================================================================
cmake_minimum_required(VERSION 3.10.0)
#==============================================================================
# Operating System specificities (must be put at first)
#==============================================================================
if(APPLE)
    set(CMAKE_OSX_ARCHITECTURES arm64)
endif()
#==============================================================================
project(gmds VERSION 1.0 LANGUAGES CXX)
#==============================================================================
# RPATH handling according to :
# https://gitlab.kitware.com/cmake/community/-/wikis/doc/cmake/RPATH-handling
# use, i.e. don't skip the full RPATH for the build tree
set(CMAKE_SKIP_BUILD_RPATH FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# the RPATH to be used when installing, but only if it's not a system directory
list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
if("${isSystemDir}" STREQUAL "-1")
    set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
endif("${isSystemDir}" STREQUAL "-1")

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
#if (NOT DEFINED CMAKE_INSTALL_RPATH_USE_LINK_PATH)
#  set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
#endif()
#==============================================================================
# Compilation options mainly for testing purpose
#==============================================================================
option(WITH_DOC "Build documentation"                                                   OFF)
option(WITH_TEST "Build test"                                                           OFF)
option(WITH_CODE_COVERAGE "Perform code coverage with avail. test suites"               OFF)
option(WITH_CODE_COVERAGE_HTML_REPORT "Activate html reporting during code covering"    OFF)
option(WITH_PYTHON_API "Build the gmds python API"                                      ON)
option(BUILD_SHARED_LIBS ""                                                             ON)
option(WITH_CGNS "Allows the CGNS export (blocking and aero components)"                ON)

#==============================================================================
if(BUILD_SHARED_LIBS)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/bin")
  set(CMAKE_POSITION_INDEPENDENT_CODE ON)
  # set(CMAKE_CXX_VISIBILITY_PRESET hidden)
elseif(WIN32)
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MT")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MTd")
endif()
#==============================================================================
if(WITH_CGNS)
    add_definitions(-DUSE_CGNS)
endif()
#==============================================================================
if (WITH_CODE_COVERAGE_HTML_REPORT)
    set(WITH_CODE_COVERAGE ON)
endif()
if (WITH_CODE_COVERAGE)
    set(WITH_TEST ON)
endif()
#==============================================================================
# Location of cmake modules
#==============================================================================
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/CMakeModules)
#==============================================================================
# Location of exports headers
#==============================================================================
include_directories(${CMAKE_CURRENT_BINARY_DIR}/exports)
#==============================================================================
#Tester avec External project, qui build en configure et projet à part (pas de redeploiement)
# Fait un projet dont on depend et pas inclus
if(WITH_TEST)
    enable_testing()
    find_package(GTest REQUIRED)
    set(TEST_SAMPLES_DIR "${CMAKE_SOURCE_DIR}/test_samples")

    configure_file(
            unit_test_config.h.in
            unit_test_config.h
            @ONLY)
    include_directories(${CMAKE_CURRENT_BINARY_DIR})
    set(gtest_force_shared_crt ON CACHE BOOL "Always use msvcrt.dll" ) # windows issue with gtest
endif(WITH_TEST)
#==============================================================================
# Commmand to add a new component in GMDS
#==============================================================================
#     COMP_NAME  component name
#     DIR_NAME   the subdirectory where the component code is located
#     LIB_NAME   the final lib name as generated at the end
#     MESSAGE    the message associated to the flag for activating this component
#     ACTIVATE   ON or OFF to indicate that the component must be used or not
#     COVER      ON or OFF, ON meaning that the component must be covered
macro(GMDS_ADD_COMPONENT COMP_NAME DIR_NAME LIB_NAME  MESSAGE ACTIVATE COVER)
    option(ENABLE_${COMP_NAME} ${MESSAGE} ${ACTIVATE} )
    list(APPEND GMDS_ALL_COMPONENTS  ${COMP_NAME})
    if(${COVER})
            list(APPEND GMDS_COVER_COMPONENTS ${COMP_NAME})
    endif()
    if (ENABLE_${COMP_NAME})
        list(APPEND GMDS_ACTIVE_COMPONENTS  ${COMP_NAME})
        list(APPEND GMDS_ACTIVE_LIBS  LIB_GMDS_${COMP_NAME})
        list(APPEND GMDS_ACTIVE_LIB_NAMES  ${LIB_NAME})
        list(APPEND GMDS_ACTIVE_DIRS ${DIR_NAME})
        message(STATUS "Component ${LIB_NAME} is activated")
    else()
        set(LIB_GMDS_${COMP_NAME} )
        message(STATUS "Component ${LIB_NAME} is not activated")
    endif ()

endmacro()


#==============================================================================
include(ExternalProject)
include(GNUInstallDirs)

add_subdirectory(external/ANN)
add_subdirectory(external/Predicates_psm)
add_subdirectory(external/inih)
add_subdirectory(external/OpenNL)
add_subdirectory(external/gco-v3)
find_package(glpk)
add_subdirectory(external/glpk)
add_subdirectory(external/HLBFGS)
#==============================================================================
# OPTIONAL COMPONENTS
#==============================================================================

GMDS_ADD_COMPONENT(
        BLOCKING                       # cmake variable
        blocking                        # src subdirectory name
        GMDSBlocking                    # name of the generated library
        "This module provides a data structure for a blocking representation"  # description
        ON                                 # is activated
        ON                                 # must be covered
)

GMDS_ADD_COMPONENT(
        CLAIRE                       # cmake variable
        claire                        # src subdirectory name
        GMDSClaire                    # name of the generated library
        " "  # description
        ON                                 # is activated
        ON                                 # must be covered
)

GMDS_ADD_COMPONENT(
        BLOCK_MESHER                       # cmake variable
        blockMesher                        # src subdirectory name
        GMDSBlockMesher                    # name of the generated library
        "Provide the mesh of a block structure classified on a CAD geometry"  # description
        ON                                 # is activated
        ON                                 # must be covered
)

GMDS_ADD_COMPONENT(
        FRAME                               # cmake variable
        frame                               # src subdirectory name
        GMDSFrame                           # name of the generated library
        "Frame field generation component"  # description
        ON                                  # is activated
        ON                                 # must be covered
)

GMDS_ADD_COMPONENT(
        FRAME_3D
        frame3d
        GMDSFrame3d
        "3D Frame field generation component"
        ON
        ON                                 # must be covered
)

GMDS_ADD_COMPONENT(
        SINGGRAPHBUILD
        singGraphBuild
        GMDSSingGraphBuild
        "Singularity graph building component"
        ON
        ON                                 # must be covered
)

GMDS_ADD_COMPONENT(
        DUAL_BLOCKING
        dualBlocking
        GMDSDualBlocking
        "Dual blocking component"
        OFF
        OFF                                 # must be covered
)

GMDS_ADD_COMPONENT(
        MEDUSA
        medusa
        GMDSMedusa
        "Visualization component"
        OFF
        OFF                                 # must be covered
)

GMDS_ADD_COMPONENT(
        KMDS
        kmds
        GMDSkmds
        "kokkos-based mesh structure"
        OFF
        OFF                                 # must be covered
)

GMDS_ADD_COMPONENT(
        ELG3D
        Elg3D
        GMDSelg3d
        "overlay-grid algorithm"
        OFF
        OFF                                 # must be covered
)

GMDS_ADD_COMPONENT(
        ELGMORPHING
        elgmorphing
        GMDSelgmorphing
        "morphing the mesh"
        OFF
        OFF                                 # must be covered
)

GMDS_ADD_COMPONENT(
        HYBRIDMESHADAPT
        hybridMeshAdapt
        GMDSHybridMeshAdapt
        "Hybrid remeshing component"
        ON
        ON                                 # must be covered
)

GMDS_ADD_COMPONENT(
        PADDING
        padding
        GMDSPadding
        "Smart Padding/Pillowing algorithm"
        ON
        ON                                # must be covered
)

GMDS_ADD_COMPONENT(
        POLYBLOCK
        polyblock
        GMDSPolyblock
        "Polycube_like algorithms for hexahedral blocking"
        ON
        ON                                # must be covered
)
GMDS_ADD_COMPONENT(
        GEOD_HONEY_COMB
        geodHoneyComb
        GMDSGeodHoneyComb
        "Algorithms to create honeycomb like meshes for spherical geometries"
        ON
        ON                                # must be covered
)

GMDS_ADD_COMPONENT(
        RLBLOCKING                       # cmake variable
        rlBlocking                       # src subdirectory name
        GMDSRlBlocking                    # name of the generated library
        " "  # description
        OFF                                 # is activated
        OFF                                 # must be covered
)



GMDS_ADD_COMPONENT(
        ANT_COLONY
        ant_colony
        GMDSAntColony
        "Selective padding using ACO ants colony"
        ON
        ON                                 # must be covered
)

#==============================================================================
set (GMDS_INCLUDE_DIRS APPEND)
#==============================================================================
# TESTING AND CODE COVERAGE
#
# Up to now, some components are not reliable enough to be added in the full
# code coverage process. It must be improved in the next few months
#==============================================================================
if (WITH_CODE_COVERAGE)
    message(STATUS "============================================")
    message(STATUS "============= COVERAGE MODE  ===============")
    message(STATUS "============================================")
    message(STATUS "CXX Compiler          = ${CMAKE_CXX_COMPILER}")
    message(STATUS "CXX Compiler Version  = ${CMAKE_CXX_COMPILER_VERSION}")
    message(STATUS "CMAKE_CXX_FLAGS       = ${CMAKE_CXX_FLAGS}")
    message(STATUS "CMAKE_CXX_FLAGS_DEBUG = ${CMAKE_CXX_FLAGS_DEBUG}")
    message(STATUS "============================================")
    message(STATUS "============ COVERED COMPONENTS ============")
    message(STATUS "============================================")
    foreach(c IN LISTS GMDS_COVER_COMPONENTS)
        set(ENABLE_${c} ON)
        message(STATUS "COMPONENT  : ${c}")
    endforeach()
    message(STATUS "============================================")
endif(WITH_CODE_COVERAGE)
#==============================================================================
#voir a ajouter une target coverage avec ces options liées à la target
# va m'éviter les variables locales
# target vide si pas besoin ou pas
# forcer mode non optimise plutôt
# Attention ne pas livrer avec code coverage!!
# Faire un build type Coverage, comme Debug ou Release
if (WITH_CODE_COVERAGE)
    include(CodeCoverage)
    set(CMAKE_CXX_FLAGS "-g -O0 --coverage -fprofile-arcs -ftest-coverage")
    set(CMAKE_C_FLAGS "-g -O0 --coverage -fprofile-arcs -ftest-coverage")
    set(COVERAGE_LCOV_EXCLUDES
            '/usr/*'
            '${CMAKE_SOURCE_DIR}/docs/*'
            '${CMAKE_SOURCE_DIR}/test_samples/*'
            '${CMAKE_SOURCE_DIR}/external/*'
            '${CMAKE_SOURCE_DIR}/*/tst/*'
            )

    #==============================================================================
    # Data directory if required for some tests
    set(DATA_DIR ${CMAKE_BINARY_DIR}/test_samples)
    # Directory that contains test suites

    if(WITH_CODE_COVERAGE_HTML_REPORT)
        SETUP_TARGET_FOR_COVERAGE_LCOV(
                NAME code_cover_gmds
                EXECUTABLE ctest
                ${DATA_DIR} .
        )
    else(WITH_CODE_COVERAGE_HTML_REPORT)
        SETUP_TARGET_FOR_COVERAGE_LCOV_SHORT(
                NAME code_cover_gmds
                EXECUTABLE ctest
                ${DATA_DIR}
        )
    endif(WITH_CODE_COVERAGE_HTML_REPORT)
    #==============================================================================
    #voir configure_file qui garde le lien de mise à jour, dépendances conservées
    file(COPY
            ${CMAKE_SOURCE_DIR}/test_samples/
            DESTINATION ${CMAKE_BINARY_DIR}/test_samples/)

endif(WITH_CODE_COVERAGE)


#==============================================================================
# PROJECT STRUCTURE
#==============================================================================
#for each component, we define a library name to use and a directory
set(LIB_GMDS_UTILS GMDSUtils)
add_subdirectory(utils)
set(LIB_GMDS_MATH GMDSMath)
add_subdirectory(math)
set(LIB_GMDS_IG GMDSIg)
add_subdirectory(ig)
set(LIB_GMDS_IO GMDSIo)
add_subdirectory(io)
set(LIB_GMDS_IG_ALGO GMDSIgAlgo)
add_subdirectory(igalgo)
set(LIB_GMDS_QUALITY GMDSQuality)
add_subdirectory(quality)
set(LIB_GMDS_CAD GMDSCad)
add_subdirectory(cad)
set(LIB_GMDS_CADFAC GMDSCadFac)
add_subdirectory(cadfac)
set(LIB_GMDS_SMOOTHY GMDSSmoothy)
add_subdirectory(smoothy)
set(LIB_GMDS_IG_SHEET GMDSSheet)
add_subdirectory(sheet)

if(WITH_PYTHON_API)
    add_subdirectory(pygmds)
endif ()
#==============================================================================
list(LENGTH GMDS_ACTIVE_COMPONENTS count)
math(EXPR count "${count}-1")
foreach(i RANGE ${count})
    list(GET GMDS_ACTIVE_COMPONENTS ${i} name)
    list(GET GMDS_ACTIVE_LIBS ${i} lib)
    list(GET GMDS_ACTIVE_DIRS ${i} dir)
    set(LIB_GMDS_${name} ${lib} )
    add_subdirectory(${dir})
endforeach()
#==============================================================================
# DOCUMENTATION
#==============================================================================
if (WITH_DOC)
    add_subdirectory(docs)
endif (WITH_DOC)
#==============================================================================
# Mise en place pour l'install avec les export
#==============================================================================
install(EXPORT GMDS_SUITE
        DESTINATION
        ${CMAKE_INSTALL_LIBDIR}/cmake/gmds-${PROJECT_VERSION}
        )

set(INCLUDE_INSTALL_DIR include/)
set(LIB_INSTALL_DIR lib/)

include(CMakePackageConfigHelpers)
configure_package_config_file(
        gmdsConfig.cmake.in ${CMAKE_CURRENT_BINARY_DIR}/gmdsConfig.cmake
        INSTALL_DESTINATION ${LIB_INSTALL_DIR}/gmds/cmake
        PATH_VARS INCLUDE_INSTALL_DIR
)

write_basic_package_version_file(
        ${CMAKE_CURRENT_BINARY_DIR}/gmdsConfigVersion.cmake
        VERSION 0.7.2
        COMPATIBILITY SameMajorVersion
)

install(FILES
        ${CMAKE_CURRENT_BINARY_DIR}/gmdsConfig.cmake
#       ${CMAKE_CURRENT_BINARY_DIR}/gmdsConfigVersion.cmake
        DESTINATION ${LIB_INSTALL_DIR}/gmds/cmake
)
