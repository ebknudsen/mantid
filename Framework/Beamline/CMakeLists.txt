set ( SRC_FILES
	src/ComponentInfo.cpp
	src/DetectorInfo.cpp
	src/SpectrumInfo.cpp
)

set ( INC_FILES
	inc/MantidBeamline/ComponentInfo.h
	inc/MantidBeamline/DetectorInfo.h
	inc/MantidBeamline/SpectrumInfo.h
)

set ( TEST_FILES
	ComponentInfoTest.h
	DetectorInfoTest.h
	SpectrumInfoTest.h
)

if (COVERALLS)
  foreach( loop_var ${SRC_FILES} ${INC_FILES})
    set_property(GLOBAL APPEND PROPERTY COVERAGE_SRCS "${CMAKE_CURRENT_SOURCE_DIR}/${loop_var}")
  endforeach(loop_var)
endif()

if(UNITY_BUILD)
  include(UnityBuild)
  enable_unity_build(Beamline SRC_FILES SRC_UNITY_IGNORE_FILES 10)
endif(UNITY_BUILD)

# Add the target for this directory
add_library ( Beamline ${SRC_FILES} ${INC_FILES} )
# Set the name of the generated library
set_target_properties ( Beamline PROPERTIES OUTPUT_NAME MantidBeamline
                                 COMPILE_DEFINITIONS IN_MANTID_BEAMLINE
                                 INSTALL_RPATH "@loader_path/../MacOS" )
target_include_directories( Beamline SYSTEM PUBLIC ${EIGEN3_INCLUDE_DIR} )

# Add to the 'Framework' group in VS
set_property ( TARGET Beamline PROPERTY FOLDER "MantidFramework" )

target_link_libraries ( Beamline LINK_PRIVATE ${TCMALLOC_LIBRARIES_LINKTIME} 
                        ${GSL_LIBRARIES} ${MANTIDLIBS} )

# Add the unit tests directory
add_subdirectory ( test )

###########################################################################
# Installation settings
###########################################################################

install ( TARGETS Beamline ${SYSTEM_PACKAGE_TARGET} DESTINATION ${LIB_DIR} )
