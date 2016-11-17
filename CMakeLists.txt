cmake_minimum_required(VERSION 2.8)
project(vtkbool)

set(CMAKE_CXX_FLAGS "-std=c++11 -Wno-deprecated -fPIC")

option(BUILD_PARAVIEW_PLUGIN "" On)
option(BUILD_AS_DEBUG "" Off)
option(BUILD_TESTING "" Off)

set(srcs
    vtkPolyDataBooleanFilter.cxx
    vtkPolyDataContactFilter.cxx
    Utilities.cxx
    Decomposer.cxx)

include_directories(".")

if(BUILD_AS_DEBUG)
    add_definitions(-DDEBUG)
endif()

#add_definitions(-D_DEBUG)

if(BUILD_PARAVIEW_PLUGIN)
    set(BUILD_TESTING Off CACHE BOOL "" FORCE)

    add_subdirectory(paraview_plugin ${CMAKE_BINARY_DIR}/paraview_plugin)

else()

    #unset(VTK_DIR CACHE)

    find_package(VTK REQUIRED COMPONENTS vtkFiltersSources vtkIOLegacy vtkFiltersExtraction vtkFiltersGeometry vtkFiltersModeling vtkRenderingFreeType NO_MODULE)

    if(VTK_FOUND)
        include(${VTK_USE_FILE})

        add_library(${PROJECT_NAME} SHARED ${srcs})
        target_link_libraries(${PROJECT_NAME} ${VTK_LIBRARIES})

        if(BUILD_TESTING)
            enable_testing()

            add_executable(testing testing/testing.cxx)
            target_link_libraries(testing ${PROJECT_NAME} ${VTK_LIBRARIES})

            foreach(i RANGE 13)
                add_test(Test_${i} testing ${i})
            endforeach()

        endif()

        if(VTK_WRAP_PYTHON)
            find_package(PythonLibs 2.7 REQUIRED)

            if(PYTHONLIBS_FOUND)

                include_directories(${PYTHON_INCLUDE_DIRS} ${PROJECT_SOURCE_DIR})

                include(vtkWrapPython)

                vtk_wrap_python3(${PROJECT_NAME}Python ${PROJECT_NAME}Python_srcs vtkPolyDataBooleanFilter.cxx)

                add_library(${PROJECT_NAME}PythonD ${${PROJECT_NAME}Python_srcs})
                target_link_libraries(${PROJECT_NAME}PythonD ${PROJECT_NAME})

                add_library(${PROJECT_NAME}Python MODULE ${PROJECT_NAME}PythonInit.cxx)
                target_link_libraries(${PROJECT_NAME}Python ${PROJECT_NAME}PythonD)

            endif()
        endif()

    endif(VTK_FOUND)

endif()
