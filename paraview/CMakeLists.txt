paraview_add_plugin(PolyDataBooleanFilter
  VERSION "${vtkbool_VERSION}"
  MODULES vtkBool
  MODULE_FILES "${CMAKE_CURRENT_SOURCE_DIR}/module/vtk.module")
