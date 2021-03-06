# GLOBs should be replaced with explicit listings
set ( SRC_FILES
    # Old TMP convertToMD code prepared for deprecation:
    # src/CreateMDFitWorkspace.cpp
    #end TMP
    src/AccumulateMD.cpp
    src/AndMD.cpp
    src/BaseConvertToDiffractionMDWorkspace.cpp
    src/BinMD.cpp
    src/BinaryOperationMD.cpp
    src/BooleanBinaryOperationMD.cpp
    src/BoxControllerSettingsAlgorithm.cpp
    src/CalculateCoverageDGS.cpp
    src/CentroidPeaksMD.cpp
    src/CentroidPeaksMD2.cpp
    src/ChangeQConvention.cpp
    src/CloneMDWorkspace.cpp
    src/CompactMD.cpp
    src/CompareMDWorkspaces.cpp
    src/ConvToMDBase.cpp
    src/ConvToMDEventsWS.cpp
    src/ConvToMDHistoWS.cpp
    src/ConvToMDSelector.cpp
    src/ConvertCWPDMDToSpectra.cpp
    src/ConvertCWSDExpToMomentum.cpp
    src/ConvertCWSDMDtoHKL.cpp
    src/ConvertMDHistoToMatrixWorkspace.cpp
    src/ConvertSpiceDataToRealSpace.cpp
    src/ConvertToDetectorFaceMD.cpp
    src/ConvertToDiffractionMDWorkspace.cpp
    src/ConvertToDiffractionMDWorkspace2.cpp
    src/ConvertToDiffractionMDWorkspace3.cpp
    src/ConvertToMD.cpp
    src/ConvertToMDMinMaxGlobal.cpp
    src/ConvertToMDMinMaxLocal.cpp
    src/ConvertToMDParent.cpp
    src/ConvertToReflectometryQ.cpp
    src/CreateMD.cpp
    src/CreateMDHistoWorkspace.cpp
    src/CreateMDWorkspace.cpp
    src/CutMD.cpp
    src/DisplayNormalizationSetter.cpp
    src/DivideMD.cpp
    src/EqualToMD.cpp
    src/EvaluateMDFunction.cpp
    src/ExponentialMD.cpp
    src/FakeMDEventData.cpp
    src/FindPeaksMD.cpp
    src/FitMD.cpp
    src/GetSpiceDataRawCountsFromMD.cpp
    src/GreaterThanMD.cpp
    src/IDynamicRebinning.cpp
    src/ImportMDEventWorkspace.cpp
    src/ImportMDHistoWorkspace.cpp
    src/ImportMDHistoWorkspaceBase.cpp
    src/Integrate3DEvents.cpp
    src/IntegrateEllipsoids.cpp
    src/IntegrateEllipsoidsTwoStep.cpp
    src/IntegrateFlux.cpp
    src/IntegrateMDHistoWorkspace.cpp
    src/IntegratePeaksMD.cpp
    src/IntegratePeaksMD2.cpp
    src/IntegratePeaksMDHKL.cpp
    src/IntegratePeaksCWSD.cpp
    src/InvalidParameter.cpp
    src/InvalidParameterParser.cpp
    src/LessThanMD.cpp
    src/LoadMD.cpp
    src/LoadSQW.cpp
    src/LoadSQW2.cpp
    src/LogarithmMD.cpp
    src/MDEventWSWrapper.cpp
    src/MDNormDirectSC.cpp
    src/MDNormSCD.cpp
    src/MDTransfAxisNames.cpp
    src/MDTransfFactory.cpp
    src/MDTransfModQ.cpp
    src/MDTransfNoQ.cpp
    src/MDTransfQ3D.cpp
    src/MDWSDescription.cpp
    src/MDWSTransform.cpp
    src/MaskMD.cpp
    src/MergeMD.cpp
    src/MergeMDFiles.cpp
    src/MinusMD.cpp
    src/MultiplyMD.cpp
    src/NotMD.cpp
    src/OneStepMDEW.cpp
    src/OrMD.cpp
    src/PlusMD.cpp
    src/PowerMD.cpp
    src/PreprocessDetectorsToMD.cpp
    src/Quantification/CachedExperimentInfo.cpp
    src/Quantification/FitResolutionConvolvedModel.cpp
    src/Quantification/ForegroundModel.cpp
    src/Quantification/ForegroundModelFactory.cpp
    src/Quantification/MDResolutionConvolution.cpp
    src/Quantification/MDResolutionConvolutionFactory.cpp
    src/Quantification/Models/MullerAnsatz.cpp
    src/Quantification/Models/QCoordinate.cpp
    src/Quantification/Models/Strontium122.cpp
    src/Quantification/Resolution/ModeratorChopperResolution.cpp
    src/Quantification/Resolution/TobyFitBMatrix.cpp
    src/Quantification/Resolution/TobyFitResolutionModel.cpp
    src/Quantification/Resolution/TobyFitYVector.cpp
    src/Quantification/ResolutionConvolvedCrossSection.cpp
    src/Quantification/SimulateResolutionConvolvedModel.cpp
    src/QueryMDWorkspace.cpp
    src/ReflectometryTransformKiKf.cpp
    src/ReflectometryTransformP.cpp
    src/ReflectometryTransformQxQz.cpp
    src/ReplicateMD.cpp
    src/SaveIsawQvector.cpp
    src/SaveMD.cpp
    src/SaveMD2.cpp
    src/SaveZODS.cpp
    src/SetMDFrame.cpp
    src/SetMDUsingMask.cpp
    src/SliceMD.cpp
    src/SlicingAlgorithm.cpp
    src/SmoothMD.cpp
    src/ThresholdMD.cpp
    src/TransformMD.cpp
    src/TransposeMD.cpp
    src/UnaryOperationMD.cpp
    src/UnitsConversionHelper.cpp
    src/UserFunctionMD.cpp
    src/WeightedMeanMD.cpp
    src/XorMD.cpp
    )

#set ( SRC_UNITY_IGNORE_FILES src/IDynamicRebinning.cpp
#)

set ( INC_FILES
    inc/MantidMDAlgorithms/AccumulateMD.h
    inc/MantidMDAlgorithms/AndMD.h
    inc/MantidMDAlgorithms/BaseConvertToDiffractionMDWorkspace.h
    inc/MantidMDAlgorithms/BinMD.h
    inc/MantidMDAlgorithms/BinaryOperationMD.h
    inc/MantidMDAlgorithms/BooleanBinaryOperationMD.h
    inc/MantidMDAlgorithms/BoxControllerSettingsAlgorithm.h
    inc/MantidMDAlgorithms/CalculateCoverageDGS.h
    inc/MantidMDAlgorithms/CentroidPeaksMD.h
    inc/MantidMDAlgorithms/CentroidPeaksMD2.h
    inc/MantidMDAlgorithms/ChangeQConvention.h
    inc/MantidMDAlgorithms/CloneMDWorkspace.h
    inc/MantidMDAlgorithms/CompactMD.h
    inc/MantidMDAlgorithms/CompareMDWorkspaces.h
    inc/MantidMDAlgorithms/ConvToMDBase.h
    inc/MantidMDAlgorithms/ConvertCWPDMDToSpectra.h
    inc/MantidMDAlgorithms/ConvertCWSDExpToMomentum.h
    inc/MantidMDAlgorithms/ConvertCWSDMDtoHKL.h
    inc/MantidMDAlgorithms/ConvertMDHistoToMatrixWorkspace.h
    inc/MantidMDAlgorithms/ConvertSpiceDataToRealSpace.h
    inc/MantidMDAlgorithms/ConvertToDetectorFaceMD.h
    inc/MantidMDAlgorithms/ConvertToDiffractionMDWorkspace.h
    inc/MantidMDAlgorithms/ConvertToDiffractionMDWorkspace2.h
    inc/MantidMDAlgorithms/ConvertToDiffractionMDWorkspace3.h
    inc/MantidMDAlgorithms/ConvertToMD.h
    inc/MantidMDAlgorithms/ConvertToMDMinMaxGlobal.h
    inc/MantidMDAlgorithms/ConvertToMDMinMaxLocal.h
    inc/MantidMDAlgorithms/ConvertToMDParent.h
    inc/MantidMDAlgorithms/ConvertToReflectometryQ.h
    inc/MantidMDAlgorithms/CreateMD.h
    inc/MantidMDAlgorithms/CreateMDHistoWorkspace.h
    inc/MantidMDAlgorithms/CreateMDWorkspace.h
    inc/MantidMDAlgorithms/CutMD.h
    inc/MantidMDAlgorithms/DisplayNormalizationSetter.h
    inc/MantidMDAlgorithms/DivideMD.h
    inc/MantidMDAlgorithms/DllConfig.h
    inc/MantidMDAlgorithms/EqualToMD.h
    inc/MantidMDAlgorithms/EvaluateMDFunction.h
    inc/MantidMDAlgorithms/ExponentialMD.h
    inc/MantidMDAlgorithms/FakeMDEventData.h
    inc/MantidMDAlgorithms/FindPeaksMD.h
    inc/MantidMDAlgorithms/FitMD.h
    inc/MantidMDAlgorithms/GSLFunctions.h
    inc/MantidMDAlgorithms/GetSpiceDataRawCountsFromMD.h
    inc/MantidMDAlgorithms/GreaterThanMD.h
    inc/MantidMDAlgorithms/IDynamicRebinning.h
    inc/MantidMDAlgorithms/ImportMDEventWorkspace.h
    inc/MantidMDAlgorithms/ImportMDHistoWorkspace.h
    inc/MantidMDAlgorithms/ImportMDHistoWorkspaceBase.h
    inc/MantidMDAlgorithms/Integrate3DEvents.h
    inc/MantidMDAlgorithms/IntegrateEllipsoids.h
    inc/MantidMDAlgorithms/IntegrateEllipsoidsTwoStep.h
    inc/MantidMDAlgorithms/IntegrateFlux.h
    inc/MantidMDAlgorithms/IntegrateMDHistoWorkspace.h
    inc/MantidMDAlgorithms/IntegratePeaksMD.h
    inc/MantidMDAlgorithms/IntegratePeaksMD2.h
    inc/MantidMDAlgorithms/IntegratePeaksMDHKL.h
    inc/MantidMDAlgorithms/IntegratePeaksCWSD.h
    inc/MantidMDAlgorithms/InvalidParameter.h
    inc/MantidMDAlgorithms/InvalidParameterParser.h
    inc/MantidMDAlgorithms/LessThanMD.h
    inc/MantidMDAlgorithms/LoadMD.h
    inc/MantidMDAlgorithms/LoadSQW.h
    inc/MantidMDAlgorithms/LoadSQW2.h
    inc/MantidMDAlgorithms/LogarithmMD.h
    inc/MantidMDAlgorithms/MDEventWSWrapper.h
    inc/MantidMDAlgorithms/MDNormDirectSC.h
    inc/MantidMDAlgorithms/MDNormSCD.h
    inc/MantidMDAlgorithms/MDTransfAxisNames.h
    inc/MantidMDAlgorithms/MDTransfFactory.h
    inc/MantidMDAlgorithms/MDTransfInterface.h
    inc/MantidMDAlgorithms/MDTransfModQ.h
    inc/MantidMDAlgorithms/MDTransfNoQ.h
    inc/MantidMDAlgorithms/MDTransfQ3D.h
    inc/MantidMDAlgorithms/MDWSDescription.h
    inc/MantidMDAlgorithms/MDWSTransform.h
    inc/MantidMDAlgorithms/MaskMD.h
    inc/MantidMDAlgorithms/MergeMD.h
    inc/MantidMDAlgorithms/MergeMDFiles.h
    inc/MantidMDAlgorithms/MinusMD.h
    inc/MantidMDAlgorithms/MultiplyMD.h
    inc/MantidMDAlgorithms/NotMD.h
    inc/MantidMDAlgorithms/OneStepMDEW.h
    inc/MantidMDAlgorithms/OrMD.h
    inc/MantidMDAlgorithms/PlusMD.h
    inc/MantidMDAlgorithms/PowerMD.h
    inc/MantidMDAlgorithms/PreprocessDetectorsToMD.h
    inc/MantidMDAlgorithms/Quantification/CachedExperimentInfo.h
    inc/MantidMDAlgorithms/Quantification/FitResolutionConvolvedModel.h
    inc/MantidMDAlgorithms/Quantification/ForegroundModel.h
    inc/MantidMDAlgorithms/Quantification/ForegroundModelFactory.h
    inc/MantidMDAlgorithms/Quantification/MDResolutionConvolution.h
    inc/MantidMDAlgorithms/Quantification/MDResolutionConvolutionFactory.h
    inc/MantidMDAlgorithms/Quantification/Models/MullerAnsatz.h
    inc/MantidMDAlgorithms/Quantification/Models/QCoordinate.h
    inc/MantidMDAlgorithms/Quantification/Models/Strontium122.h
    inc/MantidMDAlgorithms/Quantification/Resolution/ModeratorChopperResolution.h
    inc/MantidMDAlgorithms/Quantification/Resolution/TobyFitBMatrix.h
    inc/MantidMDAlgorithms/Quantification/Resolution/TobyFitResolutionModel.h
    inc/MantidMDAlgorithms/Quantification/Resolution/TobyFitYVector.h
    inc/MantidMDAlgorithms/Quantification/ResolutionConvolvedCrossSection.h
    inc/MantidMDAlgorithms/Quantification/SimulateResolutionConvolvedModel.h
    inc/MantidMDAlgorithms/QueryMDWorkspace.h
    inc/MantidMDAlgorithms/ReflectometryTransformKiKf.h
    inc/MantidMDAlgorithms/ReflectometryTransformP.h
    inc/MantidMDAlgorithms/ReflectometryTransformQxQz.h
    inc/MantidMDAlgorithms/ReplicateMD.h
    inc/MantidMDAlgorithms/SaveIsawQvector.h
    inc/MantidMDAlgorithms/SaveMD.h
    inc/MantidMDAlgorithms/SaveMD2.h
    inc/MantidMDAlgorithms/SaveZODS.h
    inc/MantidMDAlgorithms/SetMDFrame.h
    inc/MantidMDAlgorithms/SetMDUsingMask.h
    inc/MantidMDAlgorithms/SliceMD.h
    inc/MantidMDAlgorithms/SlicingAlgorithm.h
    inc/MantidMDAlgorithms/SmoothMD.h
    inc/MantidMDAlgorithms/ThresholdMD.h
    inc/MantidMDAlgorithms/TransformMD.h
    inc/MantidMDAlgorithms/TransposeMD.h
    inc/MantidMDAlgorithms/UnaryOperationMD.h
    inc/MantidMDAlgorithms/UnitsConversionHelper.h
    inc/MantidMDAlgorithms/Vector3DParameter.h
    inc/MantidMDAlgorithms/Vector3DParameterParser.h
    inc/MantidMDAlgorithms/WeightedMeanMD.h
    inc/MantidMDAlgorithms/XorMD.h
    )

# Test files. Other source files required.
set ( TEST_FILES
    #
    # these tests are as they test verify different parts of the CPR algorithms
    #CreateMDFitWorkspaceTest.h
    AccumulateMDTest.h
    AndMDTest.h
    BooleanBinaryOperationMDTest.h
    BoxControllerSettingsAlgorithmTest.h
    CachedExperimentInfoTest.h
    CalculateCoverageDGSTest.h
    CentroidPeaksMD2Test.h
    CentroidPeaksMDTest.h
    ChangeQConventionTest.h
    CloneMDWorkspaceTest.h
    CompactMDTest.h
    CompareMDWorkspacesTest.h
    ConvertCWPDMDToSpectraTest.h
    ConvertCWSDExpToMomentumTest.h
    ConvertCWSDMDtoHKLTest.h
    ConvertEventsToMDTest.h
    ConvertMDHistoToMatrixWorkspaceTest.h
    ConvertSpiceDataToRealSpaceTest.h
    ConvertToDetectorFaceMDTest.h
    ConvertToDiffractionMDWorkspaceTest.h
    ConvertToDiffractionMDWorkspace2Test.h
    ConvertToDiffractionMDWorkspace3Test.h
    ConvertToMDComponentsTest.h
    ConvertToMDMinMaxGlobalTest.h
    ConvertToMDMinMaxLocalTest.h
    ConvertToMDTest.h
    ConvertToQ3DdETest.h
    ConvertToReflectometryQTest.h
    CreateMDHistoWorkspaceTest.h
    CreateMDTest.h
    CreateMDWorkspaceTest.h
    CutMDTest.h
    DisplayNormalizationSetterTest.h
    DivideMDTest.h
    EqualToMDTest.h
    EvaluateMDFunctionTest.h
    ExponentialMDTest.h
    FakeMDEventDataTest.h
    FindPeaksMDTest.h
    FitMDTest.h
    FitResolutionConvolvedModelTest.h
    ForegroundModelTest.h
    GetSpiceDataRawCountsFromMDTest.h
    GreaterThanMDTest.h
    ImportMDEventWorkspaceTest.h
    ImportMDHistoWorkspaceTest.h
    Integrate3DEventsTest.h
    IntegrateEllipsoidsTest.h
    IntegrateEllipsoidsTwoStepTest.h
    IntegrateFluxTest.h
    IntegrateMDHistoWorkspaceTest.h
    IntegratePeaksMD2Test.h
    IntegratePeaksMDHKLTest.h
    IntegratePeaksMDTest.h
    IntegratePeaksCWSDTest.h
    InvalidParameterParserTest.h
    InvalidParameterTest.h
    LessThanMDTest.h
    LoadMDTest.h
    LoadSQWTest.h
    LoadSQW2Test.h
    LogarithmMDTest.h
    MDEventWSWrapperTest.h
    MDNormDirectSCTest.h
    MDNormSCDTest.h
    MDResolutionConvolutionFactoryTest.h
    MDTransfAxisNamesTest.h
    MDTransfFactoryTest.h
    MDTransfModQTest.h
    MDTransfQ3DTest.h
    MDWSDescriptionTest.h
    MDWSTransfTest.h
    MaskMDTest.h
    MergeMDFilesTest.h
    MergeMDTest.h
    MinusMDTest.h
    ModeratorChopperResolutionTest.h
    MullerAnsatzTest.h
    MultiplyMDTest.h
    NotMDTest.h
    OneStepMDEWTest.h
    OrMDTest.h
    PlusMDTest.h
    PowerMDTest.h
    PreprocessDetectorsToMDTest.h
    QueryMDWorkspaceTest.h
    ReflectometryTransformKiKfTest.h
    ReflectometryTransformPTest.h
    ReflectometryTransformQxQzTest.h
    ReplicateMDTest.h
    ResolutionConvolvedCrossSectionTest.h
    SaveIsawQvectorTest.h
    SaveMD2Test.h
    SaveMDTest.h
    SaveZODSTest.h
    SetMDFrameTest.h
    SetMDUsingMaskTest.h
    SimulateResolutionConvolvedModelTest.h
    SliceMDTest.h
    SlicingAlgorithmTest.h
    SmoothMDTest.h
    Strontium122Test.h
    ThresholdMDTest.h
    TobyFitBMatrixTest.h
    TobyFitResolutionModelTest.h
    TobyFitYVectorTest.h
    TransformMDTest.h
    TransposeMDTest.h
    UnaryOperationMDTest.h
    UnitsConversionHelperTest.h
    WeightedMeanMDTest.h
    XorMDTest.h
    )

set ( GMOCK_TEST_FILES
    BinaryOperationMDTest.h
    BinMDTest.h
    )

if (COVERALLS)
    foreach( loop_var ${SRC_FILES} ${INC_FILES})
        set_property(GLOBAL APPEND PROPERTY COVERAGE_SRCS "${CMAKE_CURRENT_SOURCE_DIR}/${loop_var}")
    endforeach(loop_var)
endif()

if(UNITY_BUILD)
    include(UnityBuild)
    enable_unity_build(MDAlgorithms SRC_FILES SRC_UNITY_IGNORE_FILES 10)
endif(UNITY_BUILD)

# Use a precompiled header where they are supported
enable_precompiled_headers( inc/MantidMDAlgorithms/PrecompiledHeader.h SRC_FILES )
# Add the target for this directory
add_library ( MDAlgorithms ${SRC_FILES} ${INC_FILES})
# Set the name of the generated library
set_target_properties ( MDAlgorithms PROPERTIES OUTPUT_NAME MantidMDAlgorithms 
    COMPILE_DEFINITIONS IN_MANTID_MDALGORITHMS
)

if (OSX_VERSION VERSION_GREATER 10.8)
    set_target_properties ( MDAlgorithms PROPERTIES INSTALL_RPATH "@loader_path/../Contents/MacOS")
endif () 

# Add to the 'Framework' group in VS
set_property ( TARGET MDAlgorithms PROPERTY FOLDER "MantidFramework" )

include_directories ( inc )
target_link_libraries ( MDAlgorithms LINK_PRIVATE ${TCMALLOC_LIBRARIES_LINKTIME} ${MANTIDLIBS} ${GSL_LIBRARIES} ${NEXUS_LIBRARIES}  )

# Add the unit tests directory
add_subdirectory ( test )

###########################################################################
# Installation settings
###########################################################################

install ( TARGETS MDAlgorithms ${SYSTEM_PACKAGE_TARGET} DESTINATION ${PLUGINS_DIR} )
