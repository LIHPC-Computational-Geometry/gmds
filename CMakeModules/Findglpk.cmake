if (GLPK_INC AND GLPK_LIB)
	message (STATUS "==>GLPK PROVIDED : GLPK_INC=${GLPK_INC} GLPK_LIB=${GLPK_LIB}")
	return ( )
else()
	if(NOT WIN32)
		find_path (GLPK_INC NAMES glpk.h PATH_SUFFIXES include)
		find_path (GLPK_LIB NAMES libglpk.so PATH_SUFFIXES lib)

		if (GLPK_INC AND GLPK_LIB)
			message (STATUS "==>GLPK FOUND : GLPK_INC=${GLPK_INC} GLPK_LIB=${GLPK_LIB}")
			set (GLPK_FOUND TRUE)
		else ( )
			message (STATUS "==> GLPK NOT FOUND.")
			message (STATUS "==> GLPK_INC=${GLPK_INC} GLPK_LIB=${GLPK_LIB}")
			unset (GLPK_FOUND)
			return ( )
		endif (GLPK_INC AND GLPK_LIB)
	endif (NOT WIN32)

#set (TETGEN_TARGET "tetgen::tetgen")
#add_library (${TETGEN_TARGET} UNKNOWN IMPORTED)
#set_target_properties (tetgen::tetgen PROPERTIES
#                       INTERFACE_INCLUDE_DIRECTORIES ${GLPK_INC}
#                       IMPORTED_LOCATION ${TETGEN_LIBRARY})
endif (GLPK_INC AND GLPK_LIB)