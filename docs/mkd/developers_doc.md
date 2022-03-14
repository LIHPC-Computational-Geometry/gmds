# Developer documentation

Once a component created, we use a [github workflow](git_workflow.md) in order to structure the code development. 
We also intensively use [unit tests](unit_testing.md) to valid our codes and also to perform code coverage when we merge developments.

The creation of a new GMDS component requires to follow the guideline given below.

## Creation of an optional module

### How to write/update a CMakeLists.txt file for a module

Be careful, each time you add a new header file in the **inc** subdirectory or a source file in the **src* subdirectory, a reference to this file must be added in the *CMakeLists.txt* file.

### How to create a new module

GMDS is structured as a set of libraries, each library being defined in a module. We follow the strict rule of "one module gives one library". In order to define a new module, and so library you have to follow the next series of actions. Let us consider that we want to define a new module named *xxx*

1. At the root directory of **gmds**, add a directory, named by your module name (always start with a small letter). Let us call it **xxx**;

2. At the root directory of **gmds**, open the file *CMakeLists.txt* and add the following command in the section entitle *OPTIONAL COMPONENTS*

```cmake
GMDS_ADD_COMPONENT(
        XXX                               # cmake variable
        xxx                               # src subdirectory name
        GMDSXxx                           # name of the generated library
        "description of the component"  # description
        ON                                  # is activated
        ON                                 # must be covered
)
```
**GMDS_ADD_COMPONENT** is a CMake macro for GMDS.
- The first parameter of the macro is the CMake variable that will be used for this component.
- The source code of the component must be in the subdirectory given as the second parameter.
- The third parameter is the name given at the component library when it will be installed. In the project and everywhere in another "CMakeLists.txt" file, the library will be accessible via the marco *${LIB_GMDS_XXX}$*
- The fourth parameter is a Boolean flag indicating that the component must be built (ON) or not (OFF)
- The fifth parameter is also a Boolean flag used to know if the module must be covered by the code coverage procedure that we use when a "push"-like command is performed.

3. Get into the **xxx** subdirectory and create 3 subdirectories named **inc**, **src** and **tst** and a *CMakeLists.txt* file. In the **inc** subdirectory create a subdirectory **gmds** that includes subdirectory **xxx* to achieve the following structure **xxx/inc/gmds/xxx/**

4. Fill the *CMakeLists.txt" file as follows:

```cmake

#==============================================================================
# LIBRARY DEFINTION (SOURCE FILES)
#==============================================================================
# Explicity used the name given in this preamble
set(GMDS_LIB ${LIB_GMDS_XXX})
set(GMDS_LIB_PREFIX gmds/xxx)

set(GMDS_INC
        inc/gmds/xxx/Class1.h
        inc/gmds/xxx/Class2.h
        )
set(GMDS_SRC
        src/Class1.cpp
        src/Class2.cpp
        )
#==============================================================================
add_library(${GMDS_LIB} ${GMDS_INC} ${GMDS_SRC})
#==============================================================================
# TARGET DEFINITION
#==============================================================================
include(GNUInstallDirs)
#LIBRARY TO INSTALL
target_link_libraries(${GMDS_LIB} PUBLIC
        ${LIB_GMDS_CAD}
        ${LIB_GMDS_IG})

#==============================================================================
# NOTHING TO UPDATE BELOW
#==============================================================================

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
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/gmds
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/gmds)

#==============================================================================
if(WITH_TEST)
#    add_subdirectory(tst)
endif(WITH_TEST)
#==============================================================================
```

This file deserves a few comments:
- You must set the two first variable using the name of your component (XXX) and the name of the subdirectory (xxx).
- Each time you create a new class, the corresponding header and src files must be added to the lists *GMDS_INC* and *GMDS_SRC* respectively.
- Note that the *inc* subidrectory is always structured with two subdirectories, i.e. *inc/gmds/xxx/* for getting a homogeneous way of accessing to gmds header files.
- The command **target_link_libraries(${GMDS_LIB} PUBLIC ....)** is very important. It is where you give the dependency of your module to other GDMS module or to external libraries.
- The remainder of the CMakefiles.txt should not be edited. It is used to generate libraries and to define the install procedure

### How to create an executable in a GMDS module
A GMDS module is associated to a library. But sometimes, for debug and test reasons, we can generate an executable. Considering the module *XXX*, you can do it by adding the next lines in the *CMakeLists.txt" file. 
