# Developer documentation

## Setting your development environment
gmds depends on many external components. In order to compile gmds they need to be installed properly in
order to be used with CMake.

Some are optional depending on the gmds functionalities we want to activate;
they include for example:
- `lcov` is used to perform code coverage locally to your computer;
- `py-pybind11` is mandatory for the python API;
- `glpk` is a linear programming solver that we used in some of our *gmds* basic components;
- `googletest` is used for our testing infrastructure;
- `cgal` is required for the blocking component.

On linux systems and macOS we suggest using [spack](https://spack.io/) for installing the dependencies. 
We use this system for our CI workflows. In a nutshell, spack allows you to install a set of libraries in a 
specific directory. You can see it as an equivalent of *Python environment*. In our context,
we will simply install the set of dependencies we need and use them in our CMake build system. As an example,
let's take a look at how to install the basic set of gmds components, plus the blocking component with the Python API and 
with the CGNS writer.

### Installation of dependencies with spack
The following procedure can be found in the [build_spack_gmds.sh](https://github.com/LIHPC-Computational-Geometry/spack_recipes/blob/main/build_spack_gmds.sh)
script. You can run the following commands from an empty directory; you will end up with three directories `spack, spack_recipes and gmds`.
Spack's [prerequisites](https://spack.readthedocs.io/en/v0.20.3/getting_started.html) for the version of your choice 
should be installed.
 
We download a spack release; optionally you can remove the configuration files 
stored in the `.spack` directory located in your home if you have previously used spack and want a fresh start:
```bash
#==========================================
# First get a spack release
git clone --depth=1 -b v0.20.3  https://github.com/spack/spack.git

#==========================================
# can be mandatory if you have already used spack on your computer
# delete the .spack directory in the home of the user in order to 
# have a fresh start
```

The way to build software in spack is described in what is called `recipes`, and ours 
are located in this project:
```bash
#==========================================
# get our recipes
git clone --branch gmds_temp --depth=1 https://github.com/LIHPC-Computational-Geometry/spack_recipes.git
````

Spack can be configured; you can modify the installation directory to shorten the paths at the cost
of disabling the possibility to have several installations of a same package with differing
options/versions (if choosing that any library `toto` installed with spack will be located in  `absolute_path/spack/opt/spack/toto/`.)
```bash
# Optional: modifying the install_tree variable to make it shorter and more human readable;
# the HASH part in install directory names is removed which can lead to collisions.
# The spack/etc/spack/defaults/config.yaml file can be modified by hand
# - in spack version 0.20
#sed -i 's#"{architecture}/{compiler.name}-{compiler.version}/{name}-{version}-{hash}"#"{name}"#g' spack/etc/spack/defaults/config.yaml
```

Spack comes with builtin recipes; still we have to register additional ones: 
- recipes for our own software, in particular gmds itself;
- superseded recipes that replace the builtin ones; we aim to keep this number low and directly contribute 
to the spack upstream project
```bash
# to register our recipes; it assumes that spack_recipes and spack are located at
# the same level. You can use the "spack repo add" commands instead of copying the repos.yaml file
#spack repo add ./spack_recipes/meshing_repo
#spack repo add ./spack_recipes/supersede_repo
cp spack_recipes/config/repos.yaml spack/etc/spack/defaults/repos.yaml
```

We will start using spack itself
```bash
#==========================================
# configure spack using spack commands; it modifies the .spack directory in the user home
source spack/share/spack/setup-env.sh
spack clean -a
```

First register cmake (`cmake` should be available in your `PATH`, 
for example by installing the system package)
```bash
# registering cmake
spack external find cmake
```

And then compilers; you need a compiler that handles both C and CXX. Check that the highest version of the
compiler returned by `spack compiler list` does indeed provide both, using the command `spack compiler info gcc`
```bash
# registering compilers
spack compiler find
# spack uses the highest version of the compiler found by default; if it is incomplete,
# for example the C compiler is installed but not the CXX one the installations will fail.
# Compilers found can be investigated by (see `spack help compiler` commands)
spack compiler list
#spack compiler info gcc
# An undesirable version can be removed by editing ~/.spack/linux/compilers.yaml or using
#spack compiler remove gcc@12
```

Now we get gmds:
```bash
# install for dev purposes
git clone git@github.com:LIHPC-Computational-Geometry/gmds.git
```

And we build gmds from the sources we just cloned, including all the necessary dependencies.
The options given as parameters in this example are :
- `~python~blocking~cgns` deactivates the corresponding variants, which will result in
a lighter gmds with less dependencies
- `dev_path` the path where the gmds sources are located
- `build_type` the build type, `Debug` being most likely what one will need when developing


```bash
# you will probably want build_type=Debug or RelWithDebInfo.
# Choose the variants you need, you can check them using `spack info gmds`.
# The dev_path option does not seem to handle relative paths.

# +mpi should actually be ok, but currently the default openmpi install fails
# It is activated by default in the hdf5 and cgns recipes, so choose not to use it if necessary
#spack install gmds+python+blocking+cgns dev_path=$PWD/gmds build_type=Debug ^cgns~mpi ^hdf5~mpi
spack install gmds~python~blocking~cgns dev_path=$PWD/gmds build_type=Debug
```

### Configuring the IDE
Now in order to develop we need to configure an IDE. This can be done in two ways:
1. Extracting the options from spack. Files named `gmds/spack-*` or `gmds/Buildxxx/spack-*` were created in the gmds directory; we can extract the necessary 
data to configure an IDE, in particular the `CMAKE_PREFIX_PATH` and the cmake options given to gmds.
```bash
# to configure an IDE
# spack created files and directories named gmds/spack-* in the gmds source tree, where the necessary
# options are set up
# get the CMAKE_PREFIX_PATH with ';' as separators
cat gmds/spack-build-env.txt  | grep CMAKE_PREFIX_PATH | awk -F "=" {'print $2'} | awk -F ";" {'print $1'} | sed 's/:/;/g'

# get the cmake options that were explicitly set by spack; add -DWITH_TEST:BOOL=ON
# to activate the tests
cat gmds/spack-configure-args.txt
```

2. Specifying the options by hand; note that you can proceed this way whatever the method used to install
the dependencies:
```cmake
-DWITH_PYTHON_API:BOOL=ON
-DENABLE_BLOCKING:BOOL=ON
-DWITH_CGNS:BOOL=ON
-DWITH_TEST:BOOL=ON
-DCMAKE_PREFIX_PATH=/absolute_path/spack/opt/spack/googletest;/absolute_path/spack/opt/spack/py-pybind11;/absolute_path/spack/opt/spack/cgal;/absolute_path/spack/opt/spack/gmp;/absolute_path/spack/opt/spack/mpfr;/absolute_path/spack/opt/spack/boost;/absolute_path/spack/opt/spack/glpk;/absolute_path/spack/opt/spack/cgns
```

## Creation of an optional module

Once a component is created, we use a [github workflow](git_workflow.md) in order to structure the code development.
We also intensively use [unit tests](unit_test.md) to valid our codes and also to perform code coverage when we merge developments.

The creation of a new GMDS component requires to follow the guideline given below.

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
# LIBRARY DEFINITION (SOURCE FILES)
#==============================================================================
# Explicitly used the name given in this preamble
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
- You must set the two first variables using the name of your component (XXX) and the name of the subdirectory (xxx).
- Each time you create a new class, the corresponding header and src files must be added to the lists *GMDS_INC* and *GMDS_SRC* respectively.
- Note that the *inc* subdirectory is always structured with two subdirectories, i.e. *inc/gmds/xxx/* for getting a homogeneous way of accessing to gmds header files.
- The command **target_link_libraries(${GMDS_LIB} PUBLIC ....)** is very important. It is where you give the dependency of your module to other GDMS module or to external libraries.
- The remainder of the CMakefiles.txt should not be edited. It is used to generate libraries and to define the install procedure

### How to create an executable in a GMDS module
A GMDS module is associated to a library, but an executable can also be produced.
Considering the module *XXX*, you can do so it by adding the next lines at the end 
of the *CMakeLists.txt* file. 

```cmake
#==============================================================================
# EXECUTABLE
#==============================================================================
add_executable(xxx src/main.cpp)
target_link_libraries(xxx PRIVATE ${GMDS_LIB})
target_compile_features(xxx PUBLIC cxx_std_14)
install(TARGETS xxx)
```