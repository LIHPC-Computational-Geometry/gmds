#==============================================================================
# LIBRARY DEFINTION (SOURCE FILES)
#==============================================================================
# Nommer tout en GMDS_MODULE_NAME, GMDS_SRC, ... dans les composants
set(GMDS_LIB ${LIB_GMDS_BLOCKING})
set(GMDS_LIB_PREFIX gmds/blocking)
#==============================================================================
set(GMDS_INC
        ${CMAKE_BINARY_DIR}/exports/${GMDS_LIB}_export.h
        inc/gmds/blocking/Blocking.h
        inc/gmds/blocking/CurvedBlocking.h
        inc/gmds/blocking/CurvedBlockingClassifier.h
        inc/gmds/blocking/InputMarkedDarts.h
        inc/gmds/blocking/SheetInsert.h
        inc/gmds/blocking/WriterDartsVTK.h
        inc/gmds/blocking/SheetCollapse.h)

if(WITH_CGNS)
    message(TEST)
    set(GMDS_INC ${GMDS_INC} inc/gmds/blocking/CGNSWriter.h
            )
endif()

set(GMDS_SRC
        src/Blocking.cpp
        src/CurvedBlocking.cpp
        src/CurvedBlockingClassifier.cpp
        src/InputMarkedDarts.cpp
        src/SheetInsert.cpp
        src/WriterDartsVTK.cpp
        src/SheetCollapse.cpp
        )

if(WITH_CGNS)
    set(GMDS_SRC ${GMDS_SRC} src/CGNSWriter.cpp
            )
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

# Deactivate the CGAL warning about not liking CMAKE_BUILD_TYPE=Debug
# because of
# TODO it does not seem to work, the message is still displayed ?
set (CGAL_DO_NOT_WARN_ABOUT_CMAKE_BUILD_TYPE TRUE)

if(WITH_CGNS)
    include (${CMAKE_SOURCE_DIR}/cmake/FindCgns.cmake)
endif()

# dependencies
if (WITH_CGAL_QT5)
    find_package(CGAL REQUIRED COMPONENTS Core Qt5)
    add_compile_definitions(CGAL_USE_BASIC_VIEWER)
else()
    find_package(CGAL REQUIRED COMPONENTS Core)
endif(WITH_CGAL_QT5)

message(STATUS "CGAL: ${CGAL_MAJOR_VERSION}.${CGAL_MINOR_VERSION}")

#LIBRARY TO INSTALL
target_link_libraries(${GMDS_LIB} PUBLIC
        ${LIB_GMDS_IG}
        ${LIB_GMDS_CADFAC}
        #${LIB_GMDS_IG_ALGO}
        ${LIB_GMDS_IO}
        #${LIB_GMDS_SMOOTHY})
        ${LIB_GMDS_MATH}
        CGAL::CGAL_Core
        )

if(WITH_CGNS)
    target_link_libraries(${GMDS_LIB} PUBLIC
            cgns::cgns
            )
endif()
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
if(WITH_CGNS)
    add_executable(blocking
            src/main.cpp)
    target_link_libraries(blocking PUBLIC ${LIB_GMDS_BLOCKING})

    add_executable(CGNSWriter src/main_cgns.cpp)
    target_link_libraries(CGNSWriter PRIVATE ${GMDS_LIB})
    target_compile_features(CGNSWriter PUBLIC cxx_std_14)
    install(TARGETS CGNSWriter)
endif(WITH_CGNS)