# - Config file for the VTKBOOL package
# It defines the following variables
#  VTKBOOL_INCLUDE_DIRS - include directories for MeshFix
#  VTKBOOL_LIBRARIES    - libraries to link against

# Compute paths
get_filename_component(VTKBOOL_CMAKE_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)
set(VTKBOOL_INCLUDE_DIRS "/home/cheriewilson/dev/vtkbool;/home/cheriewilson/dev/vtkbool/build;/home/cheriewilson/dev/vtkbool/libs/vp;/home/cheriewilson/dev/vtkbool/build/libs/vp")
include("${VTKBOOL_CMAKE_DIR}/vtkboolTargets.cmake")
set(VTKBOOL_LIBRARIES vtkbool)
set(VTKBOOL_LIBRARIES_DIR "${VTKBOOL_CMAKE_DIR}")

