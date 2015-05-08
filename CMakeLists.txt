cmake_minimum_required(VERSION 2.6)
project(vtkbool)

set(CMAKE_CXX_FLAGS "-std=c++0x -Wno-deprecated -fPIC")

option(BUILD_PARAVIEW_PLUGIN "" On)
option(BUILD_AS_DEBUG "" Off)
option(BUILD_TESTING "" Off)

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
    set(BUILD_TESTING Off CACHE BOOL "" FORCE)

    add_subdirectory(paraview_plugin ${CMAKE_BINARY_DIR}/paraview_plugin)

else()

    unset(VTK_DIR CACHE)

    find_package(VTK REQUIRED)

    if(VTK_FOUND)
        include(${VTK_USE_FILE})

        add_library(${PROJECT_NAME} SHARED ${srcs})
        target_link_libraries(${PROJECT_NAME} ${vtk_libs})

        if(BUILD_TESTING)
            enable_testing()

            add_executable(testing testing/testing.cxx)
            target_link_libraries(testing ${PROJECT_NAME} ${vtk_libs})

            add_test(Test_0 testing 0)
            add_test(Test_1 testing 1)
            add_test(Test_2 testing 2)
            add_test(Test_3 testing 3)
            add_test(Test_4 testing 4)
            add_test(Test_5 testing 5)
            add_test(Test_6 testing 6)
            add_test(Test_7 testing 7)
            add_test(Test_8 testing 8)
            add_test(Test_9 testing 9)

        endif()

        string(FIND "${VTK_LANGUAGES}" "PYTHON" VTK_PYTHON)

        if (NOT ${VTK_PYTHON} EQUAL -1)

            find_package(PythonLibs REQUIRED)

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
