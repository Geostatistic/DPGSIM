cmake_minimum_required(VERSION 2.8)

set(PROJECT_NAME lpm-ufrgs-dpgsim)

project(${PROJECT_NAME})

######################################################################
# Please set the variables with your preferences. The defaults should
# work out of the box.
######################################################################


set(CMAKE_BUILD_TYPE Release) # Debug, Release, RelWithDebInfo or MinSizeRel

if( NOT CMAKE_BUILD_TYPE )
  set( CMAKE_BUILD_TYPE Release CACHE STRING
       "Choose the type of build, options are: Debug Release"
       FORCE )
endif()

set(ALLOW_BUILD_PACKAGE 0 CACHE TYPE BOOL)

find_package(Qt5Core REQUIRED)
include_directories(${Qt5Core_INCLUDE_DIRS})
find_package(Qt5Gui REQUIRED)
include_directories(${Qt5Gui_INCLUDE_DIRS})
find_package(Qt5Xml REQUIRED)
include_directories(${Qt5Xml_INCLUDE_DIRS})
find_package(Qt5Svg REQUIRED)
include_directories(${Qt5Svg_INCLUDE_DIRS})
find_package(Qt5Designer REQUIRED)
include_directories(${Qt5Designer_INCLUDE_DIRS})
find_package(Qt5Widgets REQUIRED)
include_directories(${Qt5Widgets_INCLUDE_DIRS})
find_package(Qt5PrintSupport REQUIRED)
include_directories(${Qt5PrintSupport_INCLUDE_DIRS})
find_package(Qt5Multimedia REQUIRED)
include_directories(${Qt5Multimedia_INCLUDE_DIRS})
find_package(Qt5UiTools REQUIRED)
include_directories(${Qt5UiTools_INCLUDE_DIRS})

# SET RELEVANT VARS FOR PACKAGE
if (ALLOW_BUILD_PACKAGE)
	if (APPLE)
	elseif (UNIX)
	    SET(PACKAGE_TYPE "DEB")    # DEB or RPM
	endif (APPLE)
	
	SET(PACKAGE_MAINTAINER "Pericles Lopes Machado <pericles.raskolnikoff@gmail.com>")
	SET(PACKAGE_VENDOR "LPM/UFRGS")
	    
	SET(PACKAGE_VERSION_MAJOR "1")
	SET(PACKAGE_VERSION_MINOR "0")
	SET(PACKAGE_VERSION_PATCH "0")
	    
	SET(PACKAGE_SYSTEM_NAME "x64")
endif (ALLOW_BUILD_PACKAGE)

if (APPLE)
     # ...
elseif (UNIX)
    set(AR2GEMS_SOURCE_PATH "/home/balda/workspace/freelancejob/AlexSanto/ar2gems" CACHE TYPE STRING)
    set(AR2GEMS_BUILD_PATH "/home/balda/workspace/freelancejob/AlexSanto/ar2gems_build" CACHE TYPE STRING)
elseif (MSVC)
    set(BOOST_ROOT "d:/release-V1.49.0-VC2010" CACHE TYPE STRING)
    set(VTK_DIR "D:/workspace/VTK-master/build" CACHE TYPE STRING)
    set(AR2GEMS_SOURCE_PATH "d:/workspace/ar2gems" CACHE TYPE STRING)
    set(AR2GEMS_BUILD_PATH "d:/workspace/ar2gems_build" CACHE TYPE STRING)
endif (APPLE)

if (APPLE)
     # ...
elseif (UNIX)
    set(INSTALL_LIB_DIR "/usr/bin/plugins/designer")
    SET(INSTALL_GEOSTAT_LIB_DIR "/usr/bin/plugins/Geostat") # don't set ABS path it will fail at stage of package creation (ABS path forbiden
    SET(INSTALL_MAIN_LIB_DIR    "/usr/bin/") # don't set ABS path it will fail at stage of package creation (ABS path forbiden)
    
    SET(PACKAGE_DIR         "/usr/bin/plugins/") # don't set ABS path it will fail at stage of package creation (ABS path forbiden)
	
	SET(MY_CPACK_MUI_WELCOMEFINISHPAGE_BITMAP
    "${CMAKE_CURRENT_SOURCE_DIR}/images/lpm.bmp")
	SET(MY_CPACK_MUI_UNWELCOMEFINISHPAGE_BITMAP
    "${CMAKE_CURRENT_SOURCE_DIR}/images/lpm.bmp")

elseif (MSVC)

	
	set(CMAKE_INSTALL_PREFIX "${CMAKE_CURRENT_BINARY_DIR}") # this is needed for proper headers path install. The WIn doesn't support abs path in install()
    
    SET(INSTALL_LIB_DIR         "../LPMStore/plugins/designer") # don't set ABS path it will fail at stage of package creation (ABS path forbiden)
    SET(INSTALL_GEOSTAT_LIB_DIR "../LPMStore/plugins/Geostat") # don't set ABS path it will fail at stage of package creation (ABS path forbiden)
	SET(INSTALL_MAIN_LIB_DIR    "../LPMStore") # don't set ABS path it will fail at stage of package creation (ABS path forbiden)
    
    SET(PACKAGE_DIR         "../../../../LPMStore/plugins/") # don't set ABS path it will fail at stage of package creation (ABS path forbiden)
	SET(PLUGIN_BINARY_DIR "${PACKAGE_DIR}")
	SET(MY_CPACK_MUI_WELCOMEFINISHPAGE_BITMAP
    "${CMAKE_CURRENT_SOURCE_DIR}/images/lpm.bmp")
	SET(MY_CPACK_MUI_UNWELCOMEFINISHPAGE_BITMAP
    "${CMAKE_CURRENT_SOURCE_DIR}/images/lpm.bmp")
	
	
	
endif (APPLE)



set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules")


### AR2GEMS SETTINGS ###################################################
set(AR2GEMS_DIR "/home/pericles/sgems")


find_package(AR2GEMS REQUIRED)

FindAR2GEMS()
if (AR2GEMS_FOUND)
	include_directories(${AR2GEMS_INCLUDE_DIRS}/ar2gems)
	include_directories(${AR2GEMS_INCLUDE_DIRS}/ar2GsTL)
    include_directories(${AR2GEMS_INCLUDE_DIRS})
    message(STATUS "AR2GEMS found")
endif ()

if (AR2GEMS_BUILD_PATH)
    #include automatic generated ui_*.h headers
	include_directories(${AR2GEMS_BUILD_PATH}/ar2gems)
endif ()

include_directories()
    


### Enabling Boost ###################################################
FIND_PACKAGE( Boost REQUIRED )
#set(BOOST_ROOT "/usr/include/boost")
# Uncomment if Boost isn't in your path:
set(Boost_NO_SYSTEM_PATHS TRUE)
FIND_PACKAGE( Boost REQUIRED )
if(Boost_FOUND)
  message(STATUS "Boost found in ${BOOST_ROOT}")
  include_directories(${Boost_INCLUDE_DIRS})
endif()

### Enabling FLEX and BISON #############################################
find_package(BISON)
find_package(FLEX)

if(CMAKE_COMPILER_IS_GNUCXX)
elseif(MSVC)
  set(BISON_EXECUTABLE "win_bison.exe")
  set(FLEX_EXECUTABLE "win_flex.exe")
endif()



######################################################################
############### DON'T NEED CHANGE ANYTHING HERE ######################
######################################################################

if(CMAKE_COMPILER_IS_GNUCXX)
  set(CMAKE_CXX_FLAGS_RELEASE "-std=c++11 -O3 -march=native -fPIC")
  set(CMAKE_CXX_FLAGS_DEBUG "-std=c++11 -ggdb -march=native -fPIC")
elseif(MSVC)
  # customize flags for Visual Studio
endif()

set(BUILD_SHARED_LIBS TRUE)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")



### Enabling Qt4 #####################################################

### Enabling VTK #####################################################
#find_package(VTK REQUIRED)
find_package(VTK REQUIRED COMPONENTS
  vtkCommonCore
  vtkFiltersSources
  vtkRenderingCore
  vtkRenderingOpenGL
  vtkRenderingVolumeOpenGL
  vtkInteractionStyle
  vtkChartsCore
  vtkRenderingFreeTypeOpenGL
  vtkRenderingContextOpenGL
  vtkIOXML
  vtkRenderingQt
  vtkGUISupportQt
  vtkViewsQt
  vtkViewsContext2D
  vtkInteractionWidgets
)
  
if(VTK_FOUND)
  message(STATUS "VTK found in ${VTK_DIR}")
  include(${VTK_USE_FILE})
endif()

### Enabling Python ##################################################
find_package(PythonLibs 2.7 REQUIRED)
if(PYTHONLIBS_FOUND)
  message(STATUS "Python found in ${PYTHON_INCLUDE_DIR}")
  include_directories(${PYTHON_INCLUDE_DIR})
endif()



### UI files in binary directory #####################################
include_directories(${CMAKE_CURRENT_BINARY_DIR})

if ( MSVC )
	#set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${AR2GEMS_DIR}/lib/x64/release/plugins/Geostat )
	set(CMAKE_LIBRARY_OUTPUT_DIRECTORY_DEBUG ${AR2GEMS_DIR}/lib/x64/debug/plugins/Geostat )
	set(CMAKE_LIBRARY_OUTPUT_DIRECTORY_RELEASE ${AR2GEMS_DIR}/lib/x64/release/plugins/Geostat )
	set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ${AR2GEMS_DIR}/lib/x64/debug/plugins/Geostat )
	set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${AR2GEMS_DIR}/lib/x64/release/plugins/Geostat )	


elseif(UNIX)
	 set(INSTALL_LIB_DIR "/usr/bin/plugins/designer")
endif ( MSVC )

if (MSVC)
  add_definitions(-DPLUGINS_LPM_UFRGS_EXPORT)
  add_definitions(-DCLEVER_WIN32)
  add_definitions(-DCLEVER_WIN)
  add_definitions(-DCLEVER_MSVC)
  
  INSTALL(
	FILES 
	"${CMAKE_CURRENT_SOURCE_DIR}/src/lpm_ufrgs_dpgsim.ui"
	
	DESTINATION 
	"${INSTALL_GEOSTAT_LIB_DIR}")

elseif(UNIX)
  add_definitions(-DPLUGINS_LPM_UFRGS_EXPORT)
  

  INSTALL(
	FILES 
	"${CMAKE_CURRENT_SOURCE_DIR}/src/lpm_ufrgs_dpgsim.ui"
	
	DESTINATION 
	"${INSTALL_GEOSTAT_LIB_DIR}")

endif(MSVC)

add_subdirectory(src)

# CREATE PACKAGE    
if (ALLOW_BUILD_PACKAGE)
    # COMMON    
    SET(CPACK_DEBIAN_PACKAGE_MAINTAINER ${PACKAGE_MAINTAINER})
    
    SET(CPACK_PACKAGE_DESCRIPTION_FILE ${PACKAGE_DESCRIPTION_FILE})
    SET(CPACK_RESOURCE_FILE_LICENSE ${PACKAGE_LICENSE_FILE})
    
    SET(CPACK_PACKAGE_VERSION_MAJOR ${PACKAGE_VERSION_MAJOR})
    SET(CPACK_PACKAGE_VERSION_MINOR ${PACKAGE_VERSION_MINOR})
    SET(CPACK_PACKAGE_VERSION_PATCH ${PACKAGE_VERSION_PATCH})
    set(CPACK_PACKAGE_VENDOR ${PACKAGE_VENDOR})
		
		
    SET(CPACK_SYSTEM_NAME ${PACKAGE_SYSTEM_NAME})
        	                          
	if (APPLE)
	    # ...
	elseif (UNIX)
	    # DEB SPECIFIC
	    if (PACKAGE_TYPE STREQUAL "DEB")
	        SET(CPACK_GENERATOR "DEB")
	        SET(CPACK_DEBIAN_PACKAGE_SHLIBDEPS ON) # autogenerate dependency information (package d-shlibs should be installed)
	        SET(CPACK_DEBIAN_PACKAGE_DEPENDS "ar2gems (>= 1.0.0)")
	    endif (PACKAGE_TYPE STREQUAL "DEB")
	    
	    # RPM SPECIFIC
	    if (PACKAGE_TYPE STREQUAL "RPM")
	        SET(CPACK_GENERATOR "RPM")
	    
	    endif(PACKAGE_TYPE STREQUAL "RPM")
	              
	elseif (MSVC)

		SET(CPACK_PACKAGE_ICON "${CMAKE_CURRENT_SOURCE_DIR}\\\\images\\\\lpm.bmp")
		set(CPACK_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
	    set(CPACK_GENERATOR "NSIS")
	    
		set(CPACK_PACKAGE_FILE_NAME "LPMStore-dpgsim")
		set(CPACK_NSIS_INSTALLER_OUTFILE_NAME "LPMStore-dpgsim")
    	 
        #pass NATIVE PATHS to NSIS template
		string(REPLACE "/" "\\\\" WINNATIVE_PACKAGE_DIR ${PACKAGE_DIR})

		message(${WINNATIVE_PACKAGE_DIR})
	    set(CPACK_NSIS_AR2GEMS_PLUGIN_BINARY_DIR ${WINNATIVE_PACKAGE_DIR})
            
	endif (APPLE)
    
    INCLUDE(CPack)
        	
endif (ALLOW_BUILD_PACKAGE)
  
  

