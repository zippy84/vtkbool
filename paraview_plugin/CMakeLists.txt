find_package(ParaView REQUIRED)

if(ParaView_FOUND)
    include(${PARAVIEW_USE_FILE})

    set(VTK_DIR ${ParaView_DIR}/VTK CACHE STRING "" FORCE)

    include_directories("..")

    add_paraview_plugin(
        PolyDataBooleanFilter
        "1.1"
        SERVER_MANAGER_XML vtkPolyDataBooleanFilter.xml
        SERVER_MANAGER_SOURCES ../vtkPolyDataBooleanFilter.cxx
        SERVER_SOURCES ../vtkPolyDataContactFilter.cxx ../GeomHelper.cxx
        #GUI_RESOURCE_FILES vtkPolyDataBooleanFilterGUI.xml
    )

endif(ParaView_FOUND)
