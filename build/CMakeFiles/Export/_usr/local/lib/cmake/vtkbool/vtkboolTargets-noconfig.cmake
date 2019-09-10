#----------------------------------------------------------------
# Generated CMake target import file.
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "vp" for configuration ""
set_property(TARGET vp APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(vp PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_NOCONFIG "CXX"
  IMPORTED_LOCATION_NOCONFIG "/usr/local/lib/libvp.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS vp )
list(APPEND _IMPORT_CHECK_FILES_FOR_vp "/usr/local/lib/libvp.a" )

# Import target "decomp" for configuration ""
set_property(TARGET decomp APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(decomp PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_NOCONFIG "CXX"
  IMPORTED_LOCATION_NOCONFIG "/usr/local/lib/libdecomp.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS decomp )
list(APPEND _IMPORT_CHECK_FILES_FOR_decomp "/usr/local/lib/libdecomp.a" )

# Import target "merger" for configuration ""
set_property(TARGET merger APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(merger PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_NOCONFIG "CXX"
  IMPORTED_LOCATION_NOCONFIG "/usr/local/lib/libmerger.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS merger )
list(APPEND _IMPORT_CHECK_FILES_FOR_merger "/usr/local/lib/libmerger.a" )

# Import target "vtkbool" for configuration ""
set_property(TARGET vtkbool APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(vtkbool PROPERTIES
  IMPORTED_LOCATION_NOCONFIG "/usr/local/lib/libvtkbool.so"
  IMPORTED_SONAME_NOCONFIG "libvtkbool.so"
  )

list(APPEND _IMPORT_CHECK_TARGETS vtkbool )
list(APPEND _IMPORT_CHECK_FILES_FOR_vtkbool "/usr/local/lib/libvtkbool.so" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
