find_package(ParaView REQUIRED)

if(ParaView_FOUND)
    include(${PARAVIEW_USE_FILE})

    set(VTK_DIR ${ParaView_DIR}/VTK CACHE STRING "" FORCE)

    include_directories("..")

    set(srcs
        ../vtkPolyDataBooleanFilter.cxx
        ../vtkPolyDataContactFilter.cxx
        ../Utilities.cxx
        ../Decomposer.cxx
        ../Merger.cxx
        ../VisiblePoly.cxx)

    add_paraview_plugin(
        PolyDataBooleanFilter
        "1.5"
        SERVER_MANAGER_XML vtkPolyDataBooleanFilter.xml
        SERVER_MANAGER_SOURCES ../vtkPolyDataBooleanFilter.cxx
        SERVER_SOURCES ${srcs}
        #GUI_RESOURCE_FILES vtkPolyDataBooleanFilterGUI.xml
    )

endif(ParaView_FOUND)
