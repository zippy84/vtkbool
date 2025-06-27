set(srcs
    ../Utilities.cxx
    ../Optimize.cxx
    ../Contact.cxx
    ../Merger.cxx
    ../vtkPolyDataBooleanFilter.cxx
)

set(headers
    ../vtkPolyDataBooleanFilter.h
)

vtk_module_add_module(vtkbool
    SOURCES ${srcs}
    HEADERS ${headers}
)
