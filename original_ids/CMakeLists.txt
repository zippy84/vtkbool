cmake_minimum_required(VERSION 2.6)
project(vtkbool)

set(CMAKE_CXX_FLAGS "-std=c++0x -Wno-deprecated")

set(vtk_libs
    vtkCommon
    vtkFiltering
    vtkIO
    vtkGraphics)

find_package(VTK REQUIRED)

if(VTK_FOUND)
    include(${VTK_USE_FILE})
    include_directories(..)
    add_executable(ids ids.cxx ../vtkPolyDataBooleanFilter.cxx ../vtkPolyDataContactFilter.cxx ../GeomHelper.cxx)
    target_link_libraries(ids ${vtk_libs})
endif()
