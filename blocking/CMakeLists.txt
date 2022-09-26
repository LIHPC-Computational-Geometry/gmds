#==============================================================================
# LIBRARY DEFINTION (SOURCE FILES)
#==============================================================================
# Nommer tout en GMDS_MODULE_NAME, GMDS_SRC, ... dans les composants
set(GMDS_LIB ${LIB_GMDS_BLOCKING})
set(GMDS_LIB_PREFIX gmds/blocking)

set(GMDS_INC
        ${CMAKE_BINARY_DIR}/exports/${GMDS_LIB}_export.h
        inc/gmds/blocking/Blocking.h
        inc/gmds/blocking/CGNSWriter.h
        inc/gmds/blocking/SheetCollapse.h
        )
set(GMDS_SRC
        src/Blocking.cpp
        src/CGNSWriter.cpp
        src/SheetCollapse.cpp
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

# Deactivate the CGAL warning about not liking CMAKE_BUILD_TYPE=Debug
# because of
# TODO it does not seem to work, the message is still displayed ?
set (CGAL_DO_NOT_WARN_ABOUT_CMAKE_BUILD_TYPE TRUE)

set(CGAL_DATA_DIR "/usr/include/CGAL")

# dependencies
if (WITH_CGAL_QT5)
    find_package(CGAL REQUIRED COMPONENTS Core Qt5)
    add_compile_definitions(CGAL_USE_BASIC_VIEWER)
else()
    find_package(CGAL REQUIRED COMPONENTS Core)
endif(WITH_CGAL_QT5)

message(STATUS "CGAL: ${CGAL_MAJOR_VERSION}CGAL_MINOR_VERSION}")

set(CGNS_DIR /home/simon/soft/cgns/CGNS/build/)

enable_language(C)

find_package(HDF5)
if(HDF5_FOUND)
    set(HDF5_LIB "${HDF5_C_LIBRARIES}")
    list(APPEND EXTERNAL_INCLUDES ${HDF5_INCLUDE_DIRS})
    if(ENABLE_MED)
        find_library(MED_LIB medC)
        if(MED_LIB)
            set_config_option(HAVE_MED "Med")
            list(APPEND EXTERNAL_LIBRARIES ${MED_LIB})
        endif()
    endif()
    find_library(CGNS_LIB cgns HINTS ENV CGNS_ROOT PATH_SUFFIXES lib)
    find_path(CGNS_INC "cgnslib.h" HINTS ENV CGNS_ROOT PATH_SUFFIXES include)
    message(${CGNS_LIB})
    message(${CGNS_INC})
    if(CGNS_LIB AND CGNS_INC)
        set(HAVE_LIBCGNS true)
        list(APPEND EXTERNAL_LIBRARIES ${CGNS_LIB})
        list(APPEND EXTERNAL_INCLUDES ${CGNS_INC})
        if(ENABLE_CGNS_CPEX0045)
            set_config_option(HAVE_LIBCGNS_CPEX0045 "Cgns_CPEX0045")
        endif()
    endif()
    if(MED_LIB OR CGNS_LIB)
        list(APPEND EXTERNAL_LIBRARIES ${HDF5_LIB})
        find_library(SZ_LIB NAMES szlib sz)
        if(SZ_LIB)
            list(APPEND EXTERNAL_LIBRARIES ${SZ_LIB})
        endif()
        if(NOT HAVE_LIBZ) # necessary for non-GUI builds
            find_package(ZLIB)
            if(ZLIB_FOUND)
                #set_config_option(HAVE_LIBZ "Zlib")
                list(APPEND EXTERNAL_LIBRARIES ${ZLIB_LIBRARIES})
            endif()
        endif()
    endif()
else()
    message(STATUS "HDF5 not found")
endif()


#LIBRARY TO INSTALL
target_link_libraries(${GMDS_LIB} PUBLIC
        ${LIB_GMDS_IG}
        #${LIB_GMDS_CAD}
        #${LIB_GMDS_IG_ALGO}
        ${LIB_GMDS_IO}
        #${LIB_GMDS_SMOOTHY})
        ${LIB_GMDS_MATH}
        CGAL::CGAL_Core
        cgns
        )

if (WITH_CGAL_QT5)
    target_link_libraries(${GMDS_LIB} PUBLIC
            CGAL::CGAL_Qt5
            )
endif()

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

#==============================================================================
# EXECUTABLE
#==============================================================================
#add_executable(blockMesher src/main_blockmesher.cpp)
#target_link_libraries(blockMesher PRIVATE ${GMDS_LIB})
#target_compile_features(blockMesher PUBLIC cxx_std_14)

