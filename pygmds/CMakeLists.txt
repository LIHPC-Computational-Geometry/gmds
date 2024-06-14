cmake_minimum_required(VERSION 3.10)
#==============================================================================

find_package(pybind11 REQUIRED)

if(ENABLE_BLOCKING)
    pybind11_add_module(gmds
            src/binding_math.cpp
            src/binding_mesh.cpp
            src/binding_geometry.cpp
            src/binding_blocking.cpp
            src/gmds_facade.cpp
    )
else ()
    pybind11_add_module(gmds
            src/binding_math.cpp
            src/binding_mesh.cpp
            src/binding_geometry.cpp
            src/gmds_facade.cpp
    )
endif ()

target_link_libraries(gmds PRIVATE
        ${LIB_GMDS_IG}
        ${LIB_GMDS_CADFAC}
        ${LIB_GMDS_IO}
        )
if(ENABLE_BLOCKING)
    target_link_libraries(gmds PRIVATE
            LIB_GMDS_BLOCKING
    )
endif ()

target_compile_definitions(gmds PUBLIC cxx_std_14)
if(ENABLE_BLOCKING)
    target_compile_definitions(gmds PUBLIC ENABLE_BLOCKING=1)
endif ()

install(TARGETS gmds
        COMPONENT python
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        )
#==============================================================================
if(WITH_TEST)
    add_subdirectory(tst)
endif(WITH_TEST)