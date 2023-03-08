#==============================================================================
# LIBRARY DEFINTION (SOURCE FILES)
#==============================================================================
# Nommer tout en GMDS_MODULE_NAME, GMDS_SRC, ... dans les composants
set(GMDS_LIB ${LIB_GMDS_IO})
set(GMDS_LIB_PREFIX gmds/io)



#==============================================================================
# Compilation options mainly for testing purpose
#==============================================================================
option(WITH_LIMA "Activate LIMA format support"        OFF)
option(WITH_MESHB "Activate Meshb.7 format support"    OFF)

set(GMDS_INC
        ${CMAKE_BINARY_DIR}/exports/${GMDS_LIB}_export.h
        inc/gmds/io/IGMeshIOService.h
        inc/gmds/io/IMeshIOService.h
        inc/gmds/io/IReader.h
        inc/gmds/io/IWriter.h
        inc/gmds/io/MeditReader.h
        inc/gmds/io/VTKReader.h
        inc/gmds/io/VTKWriter.h
        inc/gmds/io/GMSHWriter.h
        inc/gmds/io/MdlReader.h)
set(GMDS_SRC
        src/IGMeshIOService.cpp
        src/IMeshIOService.cpp
        src/IReader.cpp
        src/IWriter.cpp
        src/MeditReader.cpp
        src/VTKReader.cpp
        src/VTKWriter.cpp
        src/GMSHWriter.cpp
        src/MdlReader.cpp)

if(WITH_MESHB)
   list(APPEND GMDS_INC inc/gmds/io/MeshBWriter.h)
   list(APPEND GMDS_SRC src/MeshBWriter.cpp)

   find_package(libMeshb)
   if(libMeshb_FOUND)
      MESSAGE(STATUS "Found installed Meshb.7 by the system")
  else()
      MESSAGE(STATUS "Found previously installed Meshb.7")
      list(APPEND CMAKE_PREFIX_PATH /home/ledouxf/cmakebuilds/)
      find_package(libMeshb REQUIRED)
   endif()
endif()


if(WITH_LIMA)
    find_package(Lima REQUIRED)
    list(APPEND GMDS_INC
            inc/gmds/io/LimaWriter.h
            inc/gmds/io/LimaReader.h)
    list(APPEND GMDS_SRC
            src/LimaWriter.cpp
            src/LimaReader.cpp)

endif()
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
if(WITH_LIMA)
    target_link_libraries(${GMDS_LIB} PUBLIC
            ${LIB_GMDS_IG}
            ${LIB_GMDS_SOFIANE}
            Lima::Lima)
elseif(WITH_MESHB)
        target_link_libraries(${GMDS_LIB} PUBLIC
                ${LIB_GMDS_IG}
                ${LIB_GMDS_SOFIANE}
                PRIVATE
                ${libMeshb_LIBRARIES})
else()
   target_link_libraries(${GMDS_LIB} PUBLIC
           ${LIB_GMDS_IG}
           ${LIB_GMDS_SOFIANE})
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
