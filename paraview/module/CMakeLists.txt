set(srcs
    ../../Utilities.cxx
    ../../Optimize.cxx
    ../../Contact.cxx
    ../../Merger.cxx
    ../../vtkPolyDataBooleanFilter.cxx
)

set(headers
    ../../vtkPolyDataBooleanFilter.h
)

vtk_module_add_module(vtkbool
    FORCE_STATIC
    SOURCES ${srcs}
    HEADERS ${headers}
)

paraview_add_server_manager_xmls(
    XMLS vtkPolyDataBooleanFilter.xml
)
