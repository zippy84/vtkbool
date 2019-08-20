# Generated by CMake

if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" LESS 2.5)
   message(FATAL_ERROR "CMake >= 2.6.0 required")
endif()
cmake_policy(PUSH)
cmake_policy(VERSION 2.6)
#----------------------------------------------------------------
# Generated CMake target import file.
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Protect against multiple inclusion, which would fail when already imported targets are added once more.
set(_targetsDefined)
set(_targetsNotDefined)
set(_expectedTargets)
foreach(_expectedTarget vtkbool decomp merger vp)
  list(APPEND _expectedTargets ${_expectedTarget})
  if(NOT TARGET ${_expectedTarget})
    list(APPEND _targetsNotDefined ${_expectedTarget})
  endif()
  if(TARGET ${_expectedTarget})
    list(APPEND _targetsDefined ${_expectedTarget})
  endif()
endforeach()
if("${_targetsDefined}" STREQUAL "${_expectedTargets}")
  unset(_targetsDefined)
  unset(_targetsNotDefined)
  unset(_expectedTargets)
  set(CMAKE_IMPORT_FILE_VERSION)
  cmake_policy(POP)
  return()
endif()
if(NOT "${_targetsDefined}" STREQUAL "")
  message(FATAL_ERROR "Some (but not all) targets in this export set were already defined.\nTargets Defined: ${_targetsDefined}\nTargets not yet defined: ${_targetsNotDefined}\n")
endif()
unset(_targetsDefined)
unset(_targetsNotDefined)
unset(_expectedTargets)


# Create imported target vtkbool
add_library(vtkbool SHARED IMPORTED)

set_target_properties(vtkbool PROPERTIES
  INTERFACE_LINK_LIBRARIES "vtkFiltersSources;vtkCommonComputationalGeometry;vtkCommonCore;vtksys;vtkCommonDataModel;vtkCommonMath;vtkCommonMisc;vtkCommonSystem;vtkCommonTransforms;vtkCommonExecutionModel;vtkFiltersCore;vtkFiltersGeneral;vtkIOLegacy;vtkIOCore;C:/Workarea/vcpkg/installed/x64-windows/share/vtk/../../lib/double-conversion.lib;\$<\$<NOT:\$<CONFIG:DEBUG>>:C:/Workarea/vcpkg/installed/x64-windows/share/vtk/../../lib/lz4.lib>;\$<\$<CONFIG:DEBUG>:C:/Workarea/vcpkg/installed/x64-windows/share/vtk/../../debug/lib/lz4d.lib>;C:/Workarea/vcpkg/installed/x64-windows/share/vtk/../../lib/lzma.lib;C:/Workarea/vcpkg/installed/x64-windows/share/vtk/../../lib/zlib.lib;C:/Workarea/vcpkg/installed/x64-windows/share/vtk/../../debug/lib/zlibd.lib;vtkFiltersExtraction;vtkFiltersStatistics;vtkImagingFourier;vtkImagingCore;vtkFiltersGeometry;vtkFiltersModeling;vtkRenderingFreeType;vtkRenderingCore;vtkCommonColor;C:/Workarea/vcpkg/installed/x64-windows/share/vtk/../../lib/freetype.lib;C:/Workarea/vcpkg/installed/x64-windows/share/vtk/../../debug/lib/freetyped.lib;vtkWrappingPythonCore;C:/Workarea/vcpkg/installed/x64-windows/share/vtk/../../lib/python37.lib;C:/Workarea/vcpkg/installed/x64-windows/share/vtk/../../debug/lib/python37_d.lib;vtkPythonInterpreter;vtkWrappingTools;decomp;merger"
)

# Create imported target decomp
add_library(decomp STATIC IMPORTED)

set_target_properties(decomp PROPERTIES
  INTERFACE_LINK_LIBRARIES "vp"
)

# Create imported target merger
add_library(merger STATIC IMPORTED)

set_target_properties(merger PROPERTIES
  INTERFACE_LINK_LIBRARIES "vtkFiltersSources;vtkCommonComputationalGeometry;vtkCommonCore;vtksys;vtkCommonDataModel;vtkCommonMath;vtkCommonMisc;vtkCommonSystem;vtkCommonTransforms;vtkCommonExecutionModel;vtkFiltersCore;vtkFiltersGeneral;vtkIOLegacy;vtkIOCore;C:/Workarea/vcpkg/installed/x64-windows/share/vtk/../../lib/double-conversion.lib;\$<\$<NOT:\$<CONFIG:DEBUG>>:C:/Workarea/vcpkg/installed/x64-windows/share/vtk/../../lib/lz4.lib>;\$<\$<CONFIG:DEBUG>:C:/Workarea/vcpkg/installed/x64-windows/share/vtk/../../debug/lib/lz4d.lib>;C:/Workarea/vcpkg/installed/x64-windows/share/vtk/../../lib/lzma.lib;C:/Workarea/vcpkg/installed/x64-windows/share/vtk/../../lib/zlib.lib;C:/Workarea/vcpkg/installed/x64-windows/share/vtk/../../debug/lib/zlibd.lib;vtkFiltersExtraction;vtkFiltersStatistics;vtkImagingFourier;vtkImagingCore;vtkFiltersGeometry;vtkFiltersModeling;vtkRenderingFreeType;vtkRenderingCore;vtkCommonColor;C:/Workarea/vcpkg/installed/x64-windows/share/vtk/../../lib/freetype.lib;C:/Workarea/vcpkg/installed/x64-windows/share/vtk/../../debug/lib/freetyped.lib;vtkWrappingPythonCore;C:/Workarea/vcpkg/installed/x64-windows/share/vtk/../../lib/python37.lib;C:/Workarea/vcpkg/installed/x64-windows/share/vtk/../../debug/lib/python37_d.lib;vtkPythonInterpreter;vtkWrappingTools"
)

# Create imported target vp
add_library(vp STATIC IMPORTED)

# Import target "vtkbool" for configuration "Debug"
set_property(TARGET vtkbool APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(vtkbool PROPERTIES
  IMPORTED_IMPLIB_DEBUG "C:/Workarea/vtkbool/out/build/x64-Debug/vtkbool.lib"
  IMPORTED_LOCATION_DEBUG "C:/Workarea/vtkbool/out/build/x64-Debug/vtkbool.dll"
  )

# Import target "decomp" for configuration "Debug"
set_property(TARGET decomp APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(decomp PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_DEBUG "CXX"
  IMPORTED_LOCATION_DEBUG "C:/Workarea/vtkbool/out/build/x64-Debug/libs/decomp_build/decomp.lib"
  )

# Import target "merger" for configuration "Debug"
set_property(TARGET merger APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(merger PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_DEBUG "CXX"
  IMPORTED_LOCATION_DEBUG "C:/Workarea/vtkbool/out/build/x64-Debug/libs/merger_build/merger.lib"
  )

# Import target "vp" for configuration "Debug"
set_property(TARGET vp APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(vp PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_DEBUG "CXX"
  IMPORTED_LOCATION_DEBUG "C:/Workarea/vtkbool/out/build/x64-Debug/libs/decomp_build/vp_build/vp.lib"
  )

# This file does not depend on other imported targets which have
# been exported from the same project but in a separate export set.

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
cmake_policy(POP)
