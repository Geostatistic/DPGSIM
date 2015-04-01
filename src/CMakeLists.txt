cmake_minimum_required(VERSION 2.8)

set(PROJECT_NAME lpm-dpgsim)

project(${PROJECT_NAME})

include_directories(${CMAKE_SOURCE_DIR}/src)
include_directories(${CMAKE_BINARY_DIR}/src)


foreach(dir ${PROJ_DIRS})
  add_subdirectory(${dir})
endforeach()

# Process all headers with Q_OBJECT
QT5_WRAP_CPP(lpm_dpgsim_MOC
  report.h
)

QT5_WRAP_UI(lpm_dpgsim_FORMS
  report.ui
)

## File to compile
add_library(lpm_dpgsim SHARED
  library_init.cpp
  dpgsim.cpp
  dpgsim.h
  common.h
  dpg_grid_path.h
  dpg_grid_path.cpp
  report.h
  report.cpp
  
  dpg_sequential_simulation.h

  ${lpm_dpgsim_MOC}
  ${lpm_dpgsim_FORMS}
)

target_link_libraries(lpm_dpgsim
  ${AR2GEMS_LIBRARIES}
  ${Qt5Core_LIBRARIES}
    ${Qt5Gui_LIBRARIES}
    ${Qt5Gui_LIBRARIES}
    ${Qt5Xml_LIBRARIES}
    ${Qt5Svg_LIBRARIES}
    ${Qt5Designer_LIBRARIES}
    ${Qt5Widgets_LIBRARIES}
    ${Qt5PrintSupport_LIBRARIES}
    ${Qt5Multimedia_LIBRARIES}
    ${Qt5UiTools_LIBRARIES}
  ${VTK_LIBRARIES}    
)

source_group(Qt FILES ${lpm_dpgsim_MOC} )

INSTALL(TARGETS lpm_dpgsim DESTINATION "${INSTALL_GEOSTAT_LIB_DIR}")
