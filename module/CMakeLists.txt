set(srcs
    ../vtkPolyDataBooleanFilter.cxx
    ../vtkPolyDataContactFilter.cxx
    ../Utilities.cxx)

set(headers
    ../vtkPolyDataBooleanFilter.h)

vtk_module_add_module(vtkBool
    SOURCES ${srcs}
    HEADERS ${headers}
)