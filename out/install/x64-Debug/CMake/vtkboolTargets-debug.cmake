#----------------------------------------------------------------
# Generated CMake target import file for configuration "Debug".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "vp" for configuration "Debug"
set_property(TARGET vp APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(vp PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_DEBUG "CXX"
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/lib/vp.lib"
  )

list(APPEND _IMPORT_CHECK_TARGETS vp )
list(APPEND _IMPORT_CHECK_FILES_FOR_vp "${_IMPORT_PREFIX}/lib/vp.lib" )

# Import target "decomp" for configuration "Debug"
set_property(TARGET decomp APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(decomp PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_DEBUG "CXX"
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/lib/decomp.lib"
  )

list(APPEND _IMPORT_CHECK_TARGETS decomp )
list(APPEND _IMPORT_CHECK_FILES_FOR_decomp "${_IMPORT_PREFIX}/lib/decomp.lib" )

# Import target "merger" for configuration "Debug"
set_property(TARGET merger APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(merger PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_DEBUG "CXX"
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/lib/merger.lib"
  )

list(APPEND _IMPORT_CHECK_TARGETS merger )
list(APPEND _IMPORT_CHECK_FILES_FOR_merger "${_IMPORT_PREFIX}/lib/merger.lib" )

# Import target "vtkbool" for configuration "Debug"
set_property(TARGET vtkbool APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(vtkbool PROPERTIES
  IMPORTED_IMPLIB_DEBUG "${_IMPORT_PREFIX}/lib/vtkbool.lib"
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/bin/vtkbool.dll"
  )

list(APPEND _IMPORT_CHECK_TARGETS vtkbool )
list(APPEND _IMPORT_CHECK_FILES_FOR_vtkbool "${_IMPORT_PREFIX}/lib/vtkbool.lib" "${_IMPORT_PREFIX}/bin/vtkbool.dll" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
