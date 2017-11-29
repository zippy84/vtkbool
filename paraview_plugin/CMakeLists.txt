#   Copyright 2012-2018 Ronald Römer
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.

unset(VTK_DIR CACHE)

find_package(ParaView REQUIRED)

if(ParaView_FOUND)
    include(${PARAVIEW_USE_FILE})

    set(VTK_DIR ${ParaView_DIR}/VTK CACHE STRING "" FORCE)

    include_directories("..")

    set(srcs
        ../vtkPolyDataBooleanFilter.cxx
        ../vtkPolyDataContactFilter.cxx
        ../Utilities.cxx)
        
    add_subdirectory(../dev dev_build)
    include_directories(../dev/merger ../dev/decomp ../dev/vp_new_try)

    add_paraview_plugin(
        PolyDataBooleanFilter
        "2.0"
        SERVER_MANAGER_XML vtkPolyDataBooleanFilter.xml
        SERVER_MANAGER_SOURCES ../vtkPolyDataBooleanFilter.cxx
        SERVER_SOURCES ${srcs}
        #GUI_RESOURCE_FILES vtkPolyDataBooleanFilterGUI.xml
    )
    
    target_link_libraries(PolyDataBooleanFilter PRIVATE decomp merger)

endif(ParaView_FOUND)
