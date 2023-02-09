# REM : la recherche de CGNS est bien sur effectuee a partir de CMAKE_PREFIX_PATH mais le contenu des
# repertoires racines identifies par CGNS_DIR et Cgns_ROOT est egalement evalue.
#
# The module defines the following variables:
#  CGNS_FOUND - the system has Cgns
#  CGNS_INCLUDE_DIR - where to find Cgns.hpp
#  CGNS_INCLUDE_DIRS - Cgns includes
#  CGNS_LIBRARY - where to find the Cgns library
#  CGNS_LIBRARIES - additional libraries
#  CGNS_ROOT_DIR - root dir


# message (STATUS "CGNS_INCLUDE_DIR: '${CGNS_INCLUDE_DIR}'")

message (STATUS "=====================================> CGNS_DIR: ${CGNS_DIR} Cgns_ROOT=${Cgns_ROOT}")

if (CGNS_DIR)
	find_path (CGNS_INCLUDE_DIR NAMES cgnsconfig.h cgnslib.h cgnsBuild.defs PATHS ${CGNS_DIR}/include)
elseif (Cgns_ROOT)
	find_path (CGNS_INCLUDE_DIR NAMES cgnsconfig.h cgnslib.h cgnsBuild.defs PATHS ${Cgns_ROOT}/include)
else()
	find_path (CGNS_INCLUDE_DIR NAMES cgnsconfig.h cgnslib.h cgnsBuild.defs)
endif()
	
# message( "CGNS_INCLUDE_DIR: '${CGNS_INCLUDE_DIR}'" )

set (CGNS_INCLUDE_DIRS ${CGNS_INCLUDE_DIR})

if (CGNS_DIR)
	find_library (CGNS_LIBRARY cgns PATHS ${CGNS_DIR} PATH_SUFFIXES shlib64;shlib;lib64;lib)
elseif (Cgns_ROOT)
	find_library (CGNS_LIBRARY cgns PATHS ${Cgns_ROOT} PATH_SUFFIXES shlib64;shlib;lib64;lib)
else()
	find_library (CGNS_LIBRARY cgns)
endif()

if (CGNS_INCLUDE_DIR AND CGNS_LIBRARY)
	message (STATUS "=====================================> CGNS FOUND : CGNS_INCLUDE_DIR=${CGNS_INCLUDE_DIR} CGNS_LIBRARY=${CGNS_LIBRARY}")
	set (CGNS_FOUND TRUE)
else ( )
	message (STATUS "=====================================> CGNS NOT FOUND.")
	unset (CGNS_FOUND)
	return ( )
endif (CGNS_INCLUDE_DIR AND CGNS_LIBRARY)

set (CGNS_LIBRARIES ${CGNS_LIBRARY})

# try to guess root dir from include dir
if (CGNS_INCLUDE_DIR)
  string (REGEX REPLACE "(.*)/include.*" "\\1" CGNS_ROOT_DIR ${CGNS_INCLUDE_DIR})
# try to guess root dir from library dir
elseif (CGNS_LIBRARY)
  string ( REGEX REPLACE "(.*)/lib[/|32|64].*" "\\1" CGNS_ROOT_DIR ${CGNS_LIBRARY})
endif ()

mark_as_advanced (
  CGNS_LIBRARY
  CGNS_LIBRARIES
  CGNS_INCLUDE_DIR
  CGNS_INCLUDE_DIRS
  CGNS_ROOT_DIR
)

set (CGNS_TARGET "cgns::cgns")

if (EXTENSION STREQUAL ".a")
	add_library (${CGNS_TARGET} STATIC IMPORTED)
else()
	add_library (${CGNS_TARGET} SHARED IMPORTED)
endif ()

set_target_properties (cgns::cgns PROPERTIES
		INTERFACE_INCLUDE_DIRECTORIES ${CGNS_INCLUDE_DIR}
		IMPORTED_LOCATION ${CGNS_LIBRARIES}
	)
