#==============================================================================
set(MCTS_DATA_DIR "" CACHE INTERNAL "Data directory where are input files")
#==============================================================================
# LIBRARY DEFINTION (SOURCE FILES)
#==============================================================================
# Nommer tout en GMDS_MODULE_NAME, GMDS_SRC, ... dans les composants
set(GMDS_LIB ${LIB_GMDS_MCTSBLOCK})
set(GMDS_LIB_PREFIX gmds/mctsblock)
#==============================================================================
set(GMDS_INC
        ${CMAKE_BINARY_DIR}/exports/${GMDS_LIB}_export.h
        inc/gmds/mctsblock/Blocking.h
        inc/gmds/mctsblock/BlockingClassifier.h
        inc/gmds/mctsblock/BlockingAction.h
        inc/gmds/mctsblock/BlockingRewardFunction.h
        inc/gmds/mctsblock/BlockingState.h
        inc/gmds/mctsblock/Graph.h
        spam/inc/mcts/IAction.h
        spam/inc/mcts/IRewardFunction.h
        spam/inc/mcts/IState.h
        spam/inc/mcts/MCTSAgent.h
        spam/inc/mcts/MCTSLoop.h
        spam/inc/mcts/MCTSSelectionFunction.h
        spam/inc/mcts/MCTSTree.h
        )

set(GMDS_SRC
        src/Blocking.cpp
        src/BlockingClassifier.cpp
        src/BlockingAction.cpp
        src/BlockingRewardFunction.cpp
        src/BlockingState.cpp
        src/Graph.cpp
        spam/src/MCTSAgent.cpp
        spam/src/MCTSLoop.cpp
        spam/src/MCTSSelectionFunction.cpp
        spam/src/MCTSTree.cpp
        )

configure_file(
        config.h.in
        config.h
        @ONLY)
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

#===============================================================================
# For the MCTS
find_package(nlohmann_json 3.2.0 REQUIRED)

# Deactivate the CGAL warning about not liking CMAKE_BUILD_TYPE=Debug
# because of
# TODO it does not seem to work, the message is still displayed ?
set (CGAL_DO_NOT_WARN_ABOUT_CMAKE_BUILD_TYPE TRUE)

find_package(CGAL REQUIRED COMPONENTS Core)
message(STATUS "CGAL: ${CGAL_MAJOR_VERSION}.${CGAL_MINOR_VERSION}")

#LIBRARY TO INSTALL
target_link_libraries(${GMDS_LIB} PUBLIC
        ${LIB_GMDS_IG}
        ${LIB_GMDS_CADFAC}
        ${LIB_GMDS_IO}
        ${LIB_GMDS_MATH}
        CGAL::CGAL_Core
        nlohmann_json::nlohmann_json
        )

# INCLUDE TO INSTALL
target_include_directories(${GMDS_LIB} PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/inc>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/spam/inc>
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
add_executable(MCTSBlocker src/main_blocker.cpp)
target_link_libraries(MCTSBlocker PRIVATE ${GMDS_LIB})
target_compile_features(MCTSBlocker PUBLIC cxx_std_14)

