set ( SRC_FILES
	SANSAddFiles.cpp
	SANSBackgroundCorrectionSettings.cpp
	SANSBackgroundCorrectionWidget.cpp
	SANSConstants.cpp
	SANSDiagnostics.cpp
	SANSEventSlicing.cpp
	SANSPlotSpecial.cpp
	SANSRunWindow.cpp
)

# Include files aren't required, but this makes them appear in Visual Studio
# IMPORTANT: Include files are required in the MOC_FILES set. Scroll down to find it.
set ( INC_FILES
	SANSAddFiles.h
	SANSBackgroundCorrectionSettings.h
	SANSBackgroundCorrectionWidget.h
	SANSConstants.h
	SANSDiagnostics.h
	SANSEventSlicing.h
	SANSPlotSpecial.h
	SANSRunWindow.h
)

set ( MOC_FILES 
    SANSBackgroundCorrectionWidget.h
    SANSAddFiles.h
    SANSPlotSpecial.h
    SANSRunWindow.h
    SANSEventSlicing.h
    SANSDiagnostics.h
)

set ( UI_FILES
    SANSBackgroundCorrectionWidget.ui
    SANSPlotSpecial.ui
    SANSRunWindow.ui
    SANSEventSlicing.ui
)

qt4_wrap_cpp ( MOCCED_FILES ${MOC_FILES} OPTIONS ${extra_options} )

qt4_wrap_ui ( UI_HDRS ${UI_FILES} )
include_directories ( ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR} )

set ( ALL_SRC ${SRC_FILES} ${MOCCED_FILES})

# Use a precompiled header where they are supported
enable_precompiled_headers( PrecompiledHeader.h  ALL_SRC )
add_library ( MantidScientificInterfacesISISSANS ${ALL_SRC} ${INC_FILES} ${UI_HDRS} )

# Set the name of the generated library
set_target_properties ( MantidScientificInterfacesISISSANS PROPERTIES COMPILE_DEFINITIONS IN_MANTIDQT_ISISSANS )
if ( MSVC_IDE )
#  # Add to the 'ScientificInterfaces' group in VS
  set_property ( TARGET MantidScientificInterfacesISISSANS PROPERTY FOLDER "ScientificInterfaces" )
endif()

if (OSX_VERSION VERSION_GREATER 10.8)
  set_target_properties ( MantidScientificInterfacesISISSANS PROPERTIES INSTALL_RPATH "@loader_path/../../../Contents/MacOS")
endif ()

find_Package(Threads)

target_link_libraries ( MantidScientificInterfacesISISSANS LINK_PRIVATE ${TCMALLOC_LIBRARIES_LINKTIME}
            MantidQtWidgetsCommon
            ${CORE_MANTIDLIBS}
            ${QT_LIBRARIES}
            ${QWT_LIBRARIES}
            ${POCO_LIBRARIES}
            ${Boost_LIBRARIES}
            ${JSONCPP_LIBRARIES}
            ${CMAKE_THREAD_LIBS_INIT})

###########################################################################
# Installation settings
###########################################################################

install ( TARGETS MantidScientificInterfacesISISSANS ${SYSTEM_PACKAGE_TARGET} DESTINATION ${PLUGINS_DIR}/qtplugins/mantid )