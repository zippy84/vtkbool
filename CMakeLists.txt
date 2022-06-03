# Copyright 2012-2022 Ronald Römer
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

cmake_minimum_required(VERSION 3.12 FATAL_ERROR)
project(vtkbool VERSION 3.0)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_OUTPUT_EXTENSION_REPLACE 1)

if(MSVC)
    add_compile_options(/EHsc)
    add_definitions(-D_SCL_SECURE_NO_WARNINGS)
    set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
else()
    add_compile_options(-Wall -Wextra -fPIC -Wconversion)
endif()

option(VTKBOOL_PARAVIEW "" OFF)
option(VTKBOOL_DEBUG "" OFF)
option(VTKBOOL_COVERAGE "" OFF)

mark_as_advanced(VTKBOOL_DEBUG)
mark_as_advanced(VTKBOOL_COVERAGE)

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND VTKBOOL_COVERAGE)
    set(CMAKE_C_FLAGS_PROFILE --coverage)
    set(CMAKE_CXX_FLAGS_PROFILE --coverage)
    add_link_options("--coverage")
endif()

include_directories(".")

if(VTKBOOL_DEBUG)
    add_definitions(-DDEBUG)
endif()

if(VTKBOOL_PARAVIEW)
    find_package(ParaView REQUIRED)

    if(ParaView_FOUND)
        paraview_plugin_scan(
            PLUGIN_FILES "${CMAKE_CURRENT_SOURCE_DIR}/paraview/paraview.plugin"
            PROVIDES_PLUGINS plugins
            ENABLE_BY_DEFAULT ON)

        set(BUILD_SHARED_LIBS ON)

        include(GNUInstallDirs)

        paraview_plugin_build(
            PLUGINS ${plugins})

    endif()

else()

    find_package(VTK REQUIRED COMPONENTS FiltersSources IOLegacy FiltersExtraction FiltersGeometry FiltersModeling FiltersFlowPaths OPTIONAL_COMPONENTS WrappingPythonCore NO_MODULE)

    if(VTK_FOUND)

        if(VTK_VERSION VERSION_LESS "9.0.0")
            message(FATAL_ERROR "vtkbool requires VTK 9.0.0 or newer.")
        endif()

        vtk_module_scan(
            MODULE_FILES "${CMAKE_CURRENT_SOURCE_DIR}/module/vtk.module"
            REQUEST_MODULES vtkBool
            PROVIDES_MODULES modules
            ENABLE_TESTS ON)

        set(BUILD_SHARED_LIBS ON)

        vtk_module_build(MODULES ${modules})

        if (VTK_WrappingPythonCore_FOUND)

            vtk_module_wrap_python(
                MODULES ${modules}
                PYTHON_PACKAGE "vtkbool"
                BUILD_STATIC OFF
                INSTALL_HEADERS OFF)

            include(CTest)

            include(GNUInstallDirs)

            vtk_module_python_default_destination(python_default_destination)

            add_test(NAME "import_vtkbool"
                COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/testing/test_python.py)

            set_property(TEST "import_vtkbool" APPEND PROPERTY ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/${python_default_destination}/vtkbool")

            add_test(NAME "generate_frieze"
                COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/testing/generate_frieze.py
                WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/testing)

            set_property(TEST "generate_frieze" APPEND PROPERTY ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/${python_default_destination}/vtkbool")

            add_test(NAME "test_filter"
                COMMAND ${Python3_EXECUTABLE} -m pytest
                WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/testing)

            set_property(TEST "test_filter" APPEND PROPERTY ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/${python_default_destination}/vtkbool")

            add_executable(test_merger testing/test_merger.cxx)
            target_link_libraries(test_merger PRIVATE vtkBool ${VTK_LIBRARIES})

            vtk_module_autoinit(
                TARGETS test_merger
                MODULES ${VTK_LIBRARIES}
            )

            add_test(NAME "test_merger"
                COMMAND test_merger
                WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/testing)

        endif()

    endif()

endif()
