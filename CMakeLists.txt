cmake_minimum_required(VERSION 2.6)
project(vtkbool)

set(CMAKE_CXX_FLAGS "-std=c++0x -Wno-deprecated")

option(BUILD_PARAVIEW_PLUGIN "" On)
option(BUILD_AS_DEBUG "" Off)
option(BUILD_TESTS "" Off)

set(vtk_libs
    vtkCommon
    vtkFiltering
    vtkIO
    vtkGraphics)

set(srcs
    vtkPolyDataBooleanFilter.cxx
    vtkPolyDataContactFilter.cxx
    GeomHelper.cxx)

if(BUILD_AS_DEBUG)
    add_definitions(-DDEBUG)
endif()

if(BUILD_PARAVIEW_PLUGIN)
    set(BUILD_TESTS Off CACHE BOOL "" FORCE)

    add_subdirectory(paraview_plugin ${CMAKE_BINARY_DIR}/paraview_plugin)

else()

    unset(VTK_DIR CACHE)

    find_package(VTK REQUIRED)

    if(VTK_FOUND)
        include(${VTK_USE_FILE})

        add_library(${PROJECT_NAME} SHARED ${srcs})
        target_link_libraries(${PROJECT_NAME} ${vtk_libs})

        if(BUILD_TESTS)

            add_executable(tests tests.cxx)
            target_link_libraries(tests ${PROJECT_NAME} ${vtk_libs})

        endif()

    endif(VTK_FOUND)

endif()
