#==============================================================================
# LIBRARY DEFINTION (SOURCE FILES)
#==============================================================================
# Nommer tout en GMDS_MODULE_NAME, GMDS_SRC, ... dans les composants
set(GMDS_LIB Medusa)

set(GMDS_INC
        inc/medusa/control/MediatorControl.h
        inc/medusa/control/MediatorDuBloPolicy.h
        inc/medusa/control/MediatorExplorerPolicy.h
        inc/medusa/control/MediatorPolicyItf.h
        inc/medusa/model/MedusaBackEnd.h
        inc/medusa/model/MedusaGrid.h
        inc/medusa/view/GraphicView.h
        inc/medusa/view/TextView.h
        inc/medusa/view/View.h
        inc/medusa/view/MedusaObserver.h
        inc/medusa/control/MouseInteractorStyle.h
        inc/medusa/Medusa.h
        inc/medusa/control/MedusaCommand.h
        inc/medusa/control/MedusaCommandSheet.h
        inc/medusa/control/MedusaCommandDual.h
        inc/medusa/control/MedusaCommandBlock.h
        inc/medusa/control/MedusaCommandBlockSmoothing.h
        )
set(GMDS_SRC
        src/Medusa.cpp
        src/MediatorControl.cpp
        src/MediatorExplorerPolicy.cpp
        src/MediatorDuBloPolicy.cpp
        src/GraphicView.cpp
        src/MedusaBackEnd.cpp
        src/MedusaGrid.cpp
        src/MedusaObserver.cpp
        src/MouseInteractorStyle.cpp
        src/TextView.cpp
        src/View.cpp
        src/MedusaCommandSheet.cpp
        src/MedusaCommandDual.cpp
        src/MedusaCommandBlock.cpp
        src/MedusaCommandBlockSmoothing.cpp)

#==============================================================================
add_library(${GMDS_LIB} ${GMDS_INC} ${GMDS_SRC})
#==============================================================================
# TARGET DEFINITION
#==============================================================================
set(VTK_DIR /home/simon/Dev/VTK/VTK-install/lib)
find_package(VTK)
#find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

#DEPENDENCY TO OTHER LIBS (and so targets)
target_link_libraries(${GMDS_LIB} PUBLIC
        ${LIB_GMDS_IG}
        ${LIB_GMDS_IO}
        ${LIB_GMDS_DUAL_BLOCKING}
        ${VTK_LIBRARIES}
        )
# CPP INTERFACE TO EXPORT
target_include_directories(${GMDS_LIB} PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/inc>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_PREFIX}/include>
        )

target_compile_features(${GMDS_LIB} PUBLIC cxx_std_14)
#==============================================================================
if(WITH_TEST)
    add_subdirectory(tst)
endif(WITH_TEST)

#==============================================================================
# EXECUTABLE
#==============================================================================
add_executable(MedusaExe src/main.cpp)
target_link_libraries(MedusaExe PRIVATE Medusa)
target_compile_features(MedusaExe PUBLIC cxx_std_14)
