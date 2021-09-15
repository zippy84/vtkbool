set(srcs
    ../../vtkPolyDataBooleanFilter.cxx
    ../../vtkPolyDataContactFilter.cxx
    ../../Utilities.cxx)

set(headers
    ../../vtkPolyDataBooleanFilter.h)

vtk_module_add_module(vtkBool
    FORCE_STATIC
    SOURCES ${srcs}
    HEADERS ${headers})

paraview_add_server_manager_xmls(
  XMLS vtkPolyDataBooleanFilter.xml)
