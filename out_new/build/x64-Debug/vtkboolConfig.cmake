# - Config file for the VTKBOOL package
# It defines the following variables
#  VTKBOOL_INCLUDE_DIRS - include directories for MeshFix
#  VTKBOOL_LIBRARIES    - libraries to link against

# Compute paths
get_filename_component(VTKBOOL_CMAKE_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)
set(VTKBOOL_INCLUDE_DIRS "C:/Workarea/vtkbool;C:/Workarea/vtkbool/out/build/x64-Debug")

set(VTKBOOL_LIBRARIES VTKBOOL)

