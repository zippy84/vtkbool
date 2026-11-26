set(srcs
    ../../Utilities.cxx
    ../../Contact2.cxx
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

find_package(TBB QUIET)
if(TARGET TBB::tbb)
    vtk_module_link(vtkbool PRIVATE TBB::tbb)
endif()
