cmake_minimum_required(VERSION 3.10)
#==============================================================================

find_package(pybind11 REQUIRED)

pybind11_add_module(gmds
        src/binding_math.cpp
        src/binding_mesh.cpp
        src/binding_geometry.cpp
        src/binding_blocking.cpp
        src/gmds_facade.cpp
        )

if(ENABLE_BLOCKING)
    pybind11_add_module(
            src/binding_blocking.cpp
    )
endif ()

target_link_libraries(gmds PRIVATE
        ${LIB_GMDS_IG}
        ${LIB_GMDS_CADFAC}
        ${LIB_GMDS_BLOCKING}
        ${LIB_GMDS_IO}
        )
target_compile_definitions(gmds PUBLIC cxx_std_14)

install(TARGETS gmds
        COMPONENT python
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        )
#==============================================================================
if(WITH_TEST)
    add_subdirectory(tst)
endif(WITH_TEST)