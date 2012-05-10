#ifndef CHECKWORKSPACESMATCHTEST_H_
#define CHECKWORKSPACESMATCHTEST_H_

#include <cxxtest/TestSuite.h>

#include "MantidAlgorithms/CheckWorkspacesMatch.h"
#include "MantidAPI/SpectraDetectorMap.h"
#include "MantidAPI/NumericAxis.h"
#include "MantidKernel/UnitFactory.h"
#include "MantidDataObjects/EventWorkspace.h"
#include "MantidDataObjects/PeaksWorkspace.h"
#include "MantidMDEvents/MDHistoWorkspace.h"
#include "MantidTestHelpers/WorkspaceCreationHelper.h"
#include "MantidTestHelpers/MDEventsTestHelper.h"
#include "MantidDataObjects/Workspace2D.h"
#include "MantidAlgorithms/CreatePeaksWorkspace.h"
#include "MantidAPI/MatrixWorkspace.h"
#include "MantidGeometry/Instrument.h"
#include "MantidMDEvents/MDBoxBase.h"

using namespace Mantid::Algorithms;
using namespace Mantid::API;
using namespace Mantid::DataObjects;
using namespace Mantid::Geometry;
using namespace Mantid::MDEvents;

class CheckWorkspacesMatchTest : public CxxTest::TestSuite
{
public:
  // This pair of boilerplate methods prevent the suite being created statically
  // This means the constructor isn't called when running other tests
  static CheckWorkspacesMatchTest *createSuite() { return new CheckWorkspacesMatchTest(); }
  static void destroySuite( CheckWorkspacesMatchTest *suite ) { delete suite; }

  CheckWorkspacesMatchTest() : ws1(WorkspaceCreationHelper::Create2DWorkspace123(2,2))
  {
  }
  
  void testName()
  {
    TS_ASSERT_EQUALS( checker.name(), "CheckWorkspacesMatch" );
  }

  void testVersion()
  {
    TS_ASSERT_EQUALS( checker.version(), 1 );
  }

  void testInit()
  {
    TS_ASSERT_THROWS_NOTHING( checker.initialize() );
    TS_ASSERT( checker.isInitialized() );
  }
  
  void testMatches()
  {
    if ( !checker.isInitialized() ) checker.initialize();
    
    MatrixWorkspace_sptr ws = WorkspaceCreationHelper::Create2DWorkspaceBinned(10,100);
    // A workspace had better match itself!
    TS_ASSERT_THROWS_NOTHING( checker.setProperty("Workspace1",ws) );
    TS_ASSERT_THROWS_NOTHING( checker.setProperty("Workspace2",ws) );
    
    TS_ASSERT( checker.execute() );
    TS_ASSERT_EQUALS( checker.getPropertyValue("Result"), checker.successString() );
    // Same, using the Mantid::API::equals() function
    TS_ASSERT( Mantid::API::equals(ws, ws) );
  }
  
  void testPeaks_matches()
  {
    if ( !checker.isInitialized() ) checker.initialize();

    std::string outWS1Name("CreatePeaks1WorkspaceTest_OutputWS");
    std::string outWS2Name("CreatePeaks2WorkspaceTest_OutputWS");

    Workspace2D_sptr instws = WorkspaceCreationHelper::create2DWorkspaceWithFullInstrument(2, 10);

    CreatePeaksWorkspace alg;
    TS_ASSERT_THROWS_NOTHING( alg.initialize() )
    TS_ASSERT( alg.isInitialized() )
    TS_ASSERT_THROWS_NOTHING( alg.setProperty("InstrumentWorkspace", boost::dynamic_pointer_cast<MatrixWorkspace>(instws)) );
    TS_ASSERT_THROWS_NOTHING( alg.setPropertyValue("OutputWorkspace", outWS1Name) );
    TS_ASSERT_THROWS_NOTHING( alg.setProperty("NumberOfPeaks", 13) );
    TS_ASSERT_THROWS_NOTHING( alg.execute(); )
    TS_ASSERT( alg.isExecuted() );

    TS_ASSERT_THROWS_NOTHING( alg.initialize() )
    TS_ASSERT( alg.isInitialized() )
    TS_ASSERT_THROWS_NOTHING( alg.setProperty("InstrumentWorkspace", boost::dynamic_pointer_cast<MatrixWorkspace>(instws)) );
    TS_ASSERT_THROWS_NOTHING( alg.setPropertyValue("OutputWorkspace", outWS2Name) );
    TS_ASSERT_THROWS_NOTHING( alg.setProperty("NumberOfPeaks", 13) );
    TS_ASSERT_THROWS_NOTHING( alg.execute(); )
    TS_ASSERT( alg.isExecuted() );

    // Retrieve the workspace from data service.
    PeaksWorkspace_sptr pws1, pws2;
    TS_ASSERT_THROWS_NOTHING( pws1 = boost::dynamic_pointer_cast<PeaksWorkspace>(
        AnalysisDataService::Instance().retrieve(outWS1Name) ) );
    TS_ASSERT_THROWS_NOTHING( pws2 = boost::dynamic_pointer_cast<PeaksWorkspace>(
        AnalysisDataService::Instance().retrieve(outWS2Name) ) );
    TS_ASSERT_THROWS_NOTHING( checker.setProperty("Workspace1",boost::dynamic_pointer_cast<Workspace>(pws1)) );
    TS_ASSERT_THROWS_NOTHING( checker.setProperty("Workspace2",boost::dynamic_pointer_cast<Workspace>(pws2)) );
    TS_ASSERT( checker.execute() );
    TS_ASSERT_EQUALS( checker.getPropertyValue("Result"), checker.successString() );
  }

  void testPeaks_extrapeak()
  {
    if ( !checker.isInitialized() ) checker.initialize();

    std::string outWS3Name("CreatePeaks3WorkspaceTest_OutputWS");
    std::string outWS4Name("CreatePeaks4WorkspaceTest_OutputWS");

    Workspace2D_sptr instws = WorkspaceCreationHelper::create2DWorkspaceWithFullInstrument(2, 10);

    CreatePeaksWorkspace alg;
    TS_ASSERT_THROWS_NOTHING( alg.initialize() )
    TS_ASSERT( alg.isInitialized() )
    TS_ASSERT_THROWS_NOTHING( alg.setProperty("InstrumentWorkspace", boost::dynamic_pointer_cast<MatrixWorkspace>(instws)) );
    TS_ASSERT_THROWS_NOTHING( alg.setPropertyValue("OutputWorkspace", outWS3Name) );
    TS_ASSERT_THROWS_NOTHING( alg.setProperty("NumberOfPeaks", 13) );
    TS_ASSERT_THROWS_NOTHING( alg.execute(); )
    TS_ASSERT( alg.isExecuted() );

    TS_ASSERT_THROWS_NOTHING( alg.initialize() )
    TS_ASSERT( alg.isInitialized() )
    TS_ASSERT_THROWS_NOTHING( alg.setProperty("InstrumentWorkspace", boost::dynamic_pointer_cast<MatrixWorkspace>(instws)) );
    TS_ASSERT_THROWS_NOTHING( alg.setPropertyValue("OutputWorkspace", outWS4Name) );
    TS_ASSERT_THROWS_NOTHING( alg.setProperty("NumberOfPeaks", 14) );
    TS_ASSERT_THROWS_NOTHING( alg.execute(); )
    TS_ASSERT( alg.isExecuted() );

    // Retrieve the workspace from data service.
    PeaksWorkspace_sptr pws1, pws2;
    TS_ASSERT_THROWS_NOTHING( pws1 = boost::dynamic_pointer_cast<PeaksWorkspace>(
        AnalysisDataService::Instance().retrieve(outWS3Name) ) );
    TS_ASSERT_THROWS_NOTHING( pws2 = boost::dynamic_pointer_cast<PeaksWorkspace>(
        AnalysisDataService::Instance().retrieve(outWS4Name) ) );
    TS_ASSERT_EQUALS(pws1->getNumberPeaks(), 13);
    TS_ASSERT_EQUALS(pws2->getNumberPeaks(), 14);
    TS_ASSERT_THROWS_NOTHING( checker.setProperty("Workspace1",boost::dynamic_pointer_cast<Workspace>(pws1)) );
    TS_ASSERT_THROWS_NOTHING( checker.setProperty("Workspace2",boost::dynamic_pointer_cast<Workspace>(pws2)) );
    TS_ASSERT( checker.execute() );
    TS_ASSERT_DIFFERS( checker.getPropertyValue("Result"), checker.successString() );
  }

  void testEvent_matches()
  {
    if ( !checker.isInitialized() ) checker.initialize();

    EventWorkspace_sptr ews1 = WorkspaceCreationHelper::CreateEventWorkspace(10,20,30);
    EventWorkspace_sptr ews2 = WorkspaceCreationHelper::CreateEventWorkspace(10,20,30);
    TS_ASSERT_THROWS_NOTHING( checker.setProperty("Workspace1",boost::dynamic_pointer_cast<MatrixWorkspace>(ews1)) );
    TS_ASSERT_THROWS_NOTHING( checker.setProperty("Workspace2",boost::dynamic_pointer_cast<MatrixWorkspace>(ews2)) );
    TS_ASSERT( checker.execute() );
    TS_ASSERT_EQUALS( checker.getPropertyValue("Result"), checker.successString() );

    // Same, using the Mantid::API::equals() function
    TS_ASSERT( Mantid::API::equals(ews1, ews2) );
  }

  void testEvent_different_type()
  {
    if ( !checker.isInitialized() ) checker.initialize();

    EventWorkspace_sptr ews2 = WorkspaceCreationHelper::CreateEventWorkspace(10,20,30);
    TS_ASSERT_THROWS_NOTHING( checker.setProperty("Workspace1",ws1) );
    TS_ASSERT_THROWS_NOTHING( checker.setProperty("Workspace2",boost::dynamic_pointer_cast<MatrixWorkspace>(ews2)) );
    TS_ASSERT( checker.execute() );
    TS_ASSERT_DIFFERS( checker.getPropertyValue("Result"), checker.successString() );

    // Same, using the Mantid::API::equals() function
    TS_ASSERT( !Mantid::API::equals(ws1, ews2) );

  }

  void testEvent_different_number_histograms()
  {
    if ( !checker.isInitialized() ) checker.initialize();

    EventWorkspace_sptr ews1 = WorkspaceCreationHelper::CreateEventWorkspace(10,20,30);
    EventWorkspace_sptr ews2 = WorkspaceCreationHelper::CreateEventWorkspace(15,20,30);
    TS_ASSERT_THROWS_NOTHING( checker.setProperty("Workspace1",boost::dynamic_pointer_cast<MatrixWorkspace>(ews1)) );
    TS_ASSERT_THROWS_NOTHING( checker.setProperty("Workspace2",boost::dynamic_pointer_cast<MatrixWorkspace>(ews2)) );
    TS_ASSERT( checker.execute() );
    TS_ASSERT_DIFFERS( checker.getPropertyValue("Result"), checker.successString() );
    // Same, using the !Mantid::API::equals() function
    TS_ASSERT( (!Mantid::API::equals(ews1, ews2)) );
  }

  void testEvent_differentEventLists()
  {
    if ( !checker.isInitialized() ) checker.initialize();
    EventWorkspace_sptr ews1 = WorkspaceCreationHelper::CreateEventWorkspace(10,20,30);
    EventWorkspace_sptr ews2 = WorkspaceCreationHelper::CreateEventWorkspace(10,20,30, 0.0, 1.0, 2);
    TS_ASSERT_THROWS_NOTHING( checker.setProperty("Workspace1",boost::dynamic_pointer_cast<MatrixWorkspace>(ews1)) );
    TS_ASSERT_THROWS_NOTHING( checker.setProperty("Workspace2",boost::dynamic_pointer_cast<MatrixWorkspace>(ews2)) );
    TS_ASSERT( checker.execute() );
    TS_ASSERT_DIFFERS( checker.getPropertyValue("Result"), checker.successString() );
    // Same, using the !Mantid::API::equals() function
    TS_ASSERT( (!Mantid::API::equals(ews1, ews2)) );
  }

  void testEvent_differentBinBoundaries()
  {
    if ( !checker.isInitialized() ) checker.initialize();
    EventWorkspace_sptr ews1 = WorkspaceCreationHelper::CreateEventWorkspace(10,20,30, 15.0, 10.0);
    EventWorkspace_sptr ews2 = WorkspaceCreationHelper::CreateEventWorkspace(10,20,30, 5.0, 10.0);
    TS_ASSERT_THROWS_NOTHING( checker.setProperty("Workspace1",boost::dynamic_pointer_cast<MatrixWorkspace>(ews1)) );
    TS_ASSERT_THROWS_NOTHING( checker.setProperty("Workspace2",boost::dynamic_pointer_cast<MatrixWorkspace>(ews2)) );
    TS_ASSERT( checker.execute() );
    TS_ASSERT_DIFFERS( checker.getPropertyValue("Result"), checker.successString() );
    // Same, using the !Mantid::API::equals() function
    TS_ASSERT( (!Mantid::API::equals(ews1, ews2)) );
  }

  void testMDEvents_matches()
  {
    if ( !checker.isInitialized() ) checker.initialize();
    MDEventWorkspace3Lean::sptr mdews1 = MDEventsTestHelper::makeFileBackedMDEW("mdev1", false);
    MDEventWorkspace3Lean::sptr mdews2 = MDEventsTestHelper::makeFileBackedMDEW("mdev2", false);
    TS_ASSERT_THROWS_NOTHING( checker.setProperty("Workspace1", boost::dynamic_pointer_cast<IMDWorkspace>(mdews1)) );
    TS_ASSERT_THROWS_NOTHING( checker.setProperty("Workspace2", boost::dynamic_pointer_cast<IMDWorkspace>(mdews2)) );
    TS_ASSERT( checker.execute() );
    TS_ASSERT_EQUALS( checker.getPropertyValue("Result"), checker.successString() );
  }

  void testMDEvents_different_eventtypes()
  {
    if ( !checker.isInitialized() ) checker.initialize();
    MDEventWorkspace3Lean::sptr mdews1 = MDEventsTestHelper::makeAnyMDEW<MDLeanEvent<3>, 3>(2, 0.0, 10.0, 1000, "A");
    MDEventWorkspace3::sptr mdews2 = MDEventsTestHelper::makeAnyMDEW<MDEvent<3>, 3>(2, 0.0, 10.0, 1000, "B");
    TS_ASSERT_THROWS_NOTHING( checker.setProperty("Workspace1", boost::dynamic_pointer_cast<IMDWorkspace>(mdews1)) );
    TS_ASSERT_THROWS_NOTHING( checker.setProperty("Workspace2", boost::dynamic_pointer_cast<IMDWorkspace>(mdews2)) );
    TS_ASSERT( checker.execute() );
    TS_ASSERT_DIFFERS( checker.getPropertyValue("Result"), checker.successString() );
  }

  void testMDEvents_different_dims()
  {
    if ( !checker.isInitialized() ) checker.initialize();
    MDEventWorkspace4Lean::sptr mdews1 = MDEventsTestHelper::makeMDEW<4>(5, -10.0, 10.0, 1);
    MDEventWorkspace3Lean::sptr mdews2 = MDEventsTestHelper::makeMDEW<3>(5, -10.0, 10.0, 1);
    TS_ASSERT_THROWS_NOTHING( checker.setProperty("Workspace1", boost::dynamic_pointer_cast<IMDWorkspace>(mdews1)) );
    TS_ASSERT_THROWS_NOTHING( checker.setProperty("Workspace2", boost::dynamic_pointer_cast<IMDWorkspace>(mdews2)) );
    TS_ASSERT( checker.execute() );
    TS_ASSERT_DIFFERS( checker.getPropertyValue("Result"), checker.successString() );
  }

  void testMDEvents_different_dimnames()
  {
    if ( !checker.isInitialized() ) checker.initialize();
    MDEventWorkspace3Lean::sptr mdews1 = MDEventsTestHelper::makeAnyMDEW<MDLeanEvent<3>, 3>(2, 0.0, 10.0, 1000, "A");
    MDEventWorkspace3Lean::sptr mdews2 = MDEventsTestHelper::makeAnyMDEW<MDLeanEvent<3>, 3>(2, 0.0, 10.0, 1000, "B", "X%d");
    TS_ASSERT_THROWS_NOTHING( checker.setProperty("Workspace1", boost::dynamic_pointer_cast<IMDWorkspace>(mdews1)) );
    TS_ASSERT_THROWS_NOTHING( checker.setProperty("Workspace2", boost::dynamic_pointer_cast<IMDWorkspace>(mdews2)) );
    TS_ASSERT( checker.execute() );
    TS_ASSERT_DIFFERS( checker.getPropertyValue("Result"), checker.successString() );
  }

  void testMDEvents_different_dimmin()
  {
    if ( !checker.isInitialized() ) checker.initialize();
    MDEventWorkspace3Lean::sptr mdews1 = MDEventsTestHelper::makeAnyMDEW<MDLeanEvent<3>, 3>(2, 0.0, 10.0, 1000, "A");
    MDEventWorkspace3Lean::sptr mdews2 = MDEventsTestHelper::makeAnyMDEW<MDLeanEvent<3>, 3>(2, 1.0, 10.0, 1000, "B");
    TS_ASSERT_THROWS_NOTHING( checker.setProperty("Workspace1", boost::dynamic_pointer_cast<IMDWorkspace>(mdews1)) );
    TS_ASSERT_THROWS_NOTHING( checker.setProperty("Workspace2", boost::dynamic_pointer_cast<IMDWorkspace>(mdews2)) );
    TS_ASSERT( checker.execute() );
    TS_ASSERT_DIFFERS( checker.getPropertyValue("Result"), checker.successString() );
  }

  void testMDEvents_different_numdata()
  {
    if ( !checker.isInitialized() ) checker.initialize();
    MDEventWorkspace3Lean::sptr mdews1 = MDEventsTestHelper::makeAnyMDEW<MDLeanEvent<3>, 3>(2, 0.0, 10.0, 1000, "A");
    MDEventWorkspace3Lean::sptr mdews2 = MDEventsTestHelper::makeAnyMDEW<MDLeanEvent<3>, 3>(2, 0.0, 10.0, 5000, "B");
    TS_ASSERT_THROWS_NOTHING( checker.setProperty("Workspace1", boost::dynamic_pointer_cast<IMDWorkspace>(mdews1)) );
    TS_ASSERT_THROWS_NOTHING( checker.setProperty("Workspace2", boost::dynamic_pointer_cast<IMDWorkspace>(mdews2)) );
    TS_ASSERT( checker.execute() );
    TS_ASSERT_DIFFERS( checker.getPropertyValue("Result"), checker.successString() );
  }

  void testMDEvents_different_data()
  {
    if ( !checker.isInitialized() ) checker.initialize();
    MDEventWorkspace3Lean::sptr mdews1 = MDEventsTestHelper::makeAnyMDEW<MDLeanEvent<3>, 3>(2, 0.0, 10.0, 1000, "A");
    MDEventWorkspace3Lean::sptr mdews2 = MDEventsTestHelper::makeAnyMDEW<MDLeanEvent<3>, 3>(2, 0.0, 10.0, 1000, "B");
    MDBoxBase<MDLeanEvent<3>, 3> *parentBox = dynamic_cast<MDBoxBase<MDLeanEvent<3>, 3> *>(mdews2->getBox());
    std::vector<MDBoxBase<MDLeanEvent<3>, 3> *> boxes;
    parentBox->getBoxes(boxes, 1000, true);
    MDBox<MDLeanEvent<3>, 3> *box = dynamic_cast<MDBox<MDLeanEvent<3>, 3> *>(boxes[0]);
    std::vector<MDLeanEvent<3> > &events = box->getEvents();
    const float offset = static_cast<const float>(0.1);
    events[0].setSignal(events[0].getSignal() + offset);
    TS_ASSERT_THROWS_NOTHING( checker.setProperty("Workspace1", boost::dynamic_pointer_cast<IMDWorkspace>(mdews1)) );
    TS_ASSERT_THROWS_NOTHING( checker.setProperty("Workspace2", boost::dynamic_pointer_cast<IMDWorkspace>(mdews2)) );
    TS_ASSERT( checker.execute() );
    TS_ASSERT_DIFFERS( checker.getPropertyValue("Result"), checker.successString() );
  }

  void testMDEvents_different_error()
  {
    if ( !checker.isInitialized() ) checker.initialize();
    MDEventWorkspace3Lean::sptr mdews1 = MDEventsTestHelper::makeAnyMDEW<MDLeanEvent<3>, 3>(2, 0.0, 10.0, 1000, "A");
    MDEventWorkspace3Lean::sptr mdews2 = MDEventsTestHelper::makeAnyMDEW<MDLeanEvent<3>, 3>(2, 0.0, 10.0, 1000, "B");
    MDBoxBase<MDLeanEvent<3>, 3> *parentBox = dynamic_cast<MDBoxBase<MDLeanEvent<3>, 3> *>(mdews2->getBox());
    std::vector<MDBoxBase<MDLeanEvent<3>, 3> *> boxes;
    parentBox->getBoxes(boxes, 1000, true);
    MDBox<MDLeanEvent<3>, 3> *box = dynamic_cast<MDBox<MDLeanEvent<3>, 3> *>(boxes[0]);
    std::vector<MDLeanEvent<3> > &events = box->getEvents();
    const float offset = static_cast<const float>(0.1);
    events[0].setErrorSquared(events[0].getErrorSquared() + offset);
    TS_ASSERT_THROWS_NOTHING( checker.setProperty("Workspace1", boost::dynamic_pointer_cast<IMDWorkspace>(mdews1)) );
    TS_ASSERT_THROWS_NOTHING( checker.setProperty("Workspace2", boost::dynamic_pointer_cast<IMDWorkspace>(mdews2)) );
    TS_ASSERT( checker.execute() );
    TS_ASSERT_DIFFERS( checker.getPropertyValue("Result"), checker.successString() );
  }

  void testMDHisto_matches()
  {
    if ( !checker.isInitialized() ) checker.initialize();
    MDHistoWorkspace_sptr mdhws1 = MDEventsTestHelper::makeFakeMDHistoWorkspace(5.0, 4);
    MDHistoWorkspace_sptr mdhws2 = MDEventsTestHelper::makeFakeMDHistoWorkspace(5.0, 4);
    TS_ASSERT_THROWS_NOTHING( checker.setProperty("Workspace1", boost::dynamic_pointer_cast<IMDWorkspace>(mdhws1)) );
    TS_ASSERT_THROWS_NOTHING( checker.setProperty("Workspace2", boost::dynamic_pointer_cast<IMDWorkspace>(mdhws2)) );
    TS_ASSERT( checker.execute() );
    TS_ASSERT_EQUALS( checker.getPropertyValue("Result"), checker.successString() );
  }

  void testMDHist_different_dims()
  {
    if ( !checker.isInitialized() ) checker.initialize();
    MDHistoWorkspace_sptr mdhws1 = MDEventsTestHelper::makeFakeMDHistoWorkspace(5.0, 4);
    MDHistoWorkspace_sptr mdhws2 = MDEventsTestHelper::makeFakeMDHistoWorkspace(5.0, 3);
    TS_ASSERT_THROWS_NOTHING( checker.setProperty("Workspace1", boost::dynamic_pointer_cast<IMDWorkspace>(mdhws1)) );
    TS_ASSERT_THROWS_NOTHING( checker.setProperty("Workspace2", boost::dynamic_pointer_cast<IMDWorkspace>(mdhws2)) );
    TS_ASSERT( checker.execute() );
    TS_ASSERT_DIFFERS( checker.getPropertyValue("Result"), checker.successString() );
  }

  void testMDHist_different_dimnames()
  {
    if ( !checker.isInitialized() ) checker.initialize();
    MDHistoWorkspace_sptr mdhws1 = MDEventsTestHelper::makeFakeMDHistoWorkspace(5.0, 3);
    const int dims = 3;
    std::size_t numBins[dims] = {10, 10, 10};
    Mantid::coord_t min[dims] = {0.0, 0.0, 0.0};
    Mantid::coord_t max[dims] = {10.0, 10.0, 10.0};
    std::vector<std::string> names;
    names.push_back("h");
    names.push_back("k");
    names.push_back("l");
    MDHistoWorkspace_sptr mdhws2 = MDEventsTestHelper::makeFakeMDHistoWorkspaceGeneral(3, 5.0, 1.0, numBins, min, max, names);
    TS_ASSERT_THROWS_NOTHING( checker.setProperty("Workspace1", boost::dynamic_pointer_cast<IMDWorkspace>(mdhws1)) );
    TS_ASSERT_THROWS_NOTHING( checker.setProperty("Workspace2", boost::dynamic_pointer_cast<IMDWorkspace>(mdhws2)) );
    TS_ASSERT( checker.execute() );
    TS_ASSERT_DIFFERS( checker.getPropertyValue("Result"), checker.successString() );
  }

  void testMDHist_different_dimbins()
  {
    if ( !checker.isInitialized() ) checker.initialize();
    MDHistoWorkspace_sptr mdhws1 = MDEventsTestHelper::makeFakeMDHistoWorkspace(5.0, 3);
    MDHistoWorkspace_sptr mdhws2 = MDEventsTestHelper::makeFakeMDHistoWorkspace(5.0, 3, 5);
    TS_ASSERT_THROWS_NOTHING( checker.setProperty("Workspace1", boost::dynamic_pointer_cast<IMDWorkspace>(mdhws1)) );
    TS_ASSERT_THROWS_NOTHING( checker.setProperty("Workspace2", boost::dynamic_pointer_cast<IMDWorkspace>(mdhws2)) );
    TS_ASSERT( checker.execute() );
    TS_ASSERT_DIFFERS( checker.getPropertyValue("Result"), checker.successString() );
  }

  void testMDHist_different_dimmax()
  {
    if ( !checker.isInitialized() ) checker.initialize();
    MDHistoWorkspace_sptr mdhws1 = MDEventsTestHelper::makeFakeMDHistoWorkspace(5.0, 3);
    Mantid::coord_t max = static_cast<Mantid::coord_t>(10.1);
    MDHistoWorkspace_sptr mdhws2 = MDEventsTestHelper::makeFakeMDHistoWorkspace(5.0, 3, 10, max);
    TS_ASSERT_THROWS_NOTHING( checker.setProperty("Workspace1", boost::dynamic_pointer_cast<IMDWorkspace>(mdhws1)) );
    TS_ASSERT_THROWS_NOTHING( checker.setProperty("Workspace2", boost::dynamic_pointer_cast<IMDWorkspace>(mdhws2)) );
    TS_ASSERT( checker.execute() );
    TS_ASSERT_DIFFERS( checker.getPropertyValue("Result"), checker.successString() );
  }

  void testMDHist_different_data()
  {
    if ( !checker.isInitialized() ) checker.initialize();
    MDHistoWorkspace_sptr mdhws1 = MDEventsTestHelper::makeFakeMDHistoWorkspace(5.0, 4);
    MDHistoWorkspace_sptr mdhws2 = MDEventsTestHelper::makeFakeMDHistoWorkspace(5.1, 4);
    TS_ASSERT_THROWS_NOTHING( checker.setProperty("Workspace1", boost::dynamic_pointer_cast<IMDWorkspace>(mdhws1)) );
    TS_ASSERT_THROWS_NOTHING( checker.setProperty("Workspace2", boost::dynamic_pointer_cast<IMDWorkspace>(mdhws2)) );
    TS_ASSERT_THROWS_NOTHING( checker.setProperty("Tolerance", 1.0e-5 ) );
    TS_ASSERT( checker.execute() );
    TS_ASSERT_DIFFERS( checker.getPropertyValue("Result"), checker.successString() );
  }

  void testMDHist_different_error()
  {
    if ( !checker.isInitialized() ) checker.initialize();
    MDHistoWorkspace_sptr mdhws1 = MDEventsTestHelper::makeFakeMDHistoWorkspace(5.0, 4);
    MDHistoWorkspace_sptr mdhws2 = MDEventsTestHelper::makeFakeMDHistoWorkspace(5.0, 4, 10, 10.0, 1.1);
    TS_ASSERT_THROWS_NOTHING( checker.setProperty("Workspace1", boost::dynamic_pointer_cast<IMDWorkspace>(mdhws1)) );
    TS_ASSERT_THROWS_NOTHING( checker.setProperty("Workspace2", boost::dynamic_pointer_cast<IMDWorkspace>(mdhws2)) );
    TS_ASSERT_THROWS_NOTHING( checker.setProperty("Tolerance", 1.0e-5 ) );
    TS_ASSERT( checker.execute() );
    TS_ASSERT_DIFFERS( checker.getPropertyValue("Result"), checker.successString() );
  }

  void testDifferentSize()
  {
    if ( !checker.isInitialized() ) checker.initialize();

    Mantid::API::MatrixWorkspace_sptr ws2 = WorkspaceCreationHelper::Create1DWorkspaceFib(2);
    
    TS_ASSERT_THROWS_NOTHING( checker.setProperty("Workspace1",ws1) );
    TS_ASSERT_THROWS_NOTHING( checker.setProperty("Workspace2",ws2) );
    
    TS_ASSERT( checker.execute() );
    TS_ASSERT_EQUALS( checker.getPropertyValue("Result"), "Size mismatch" );
    
    // Same, using the !Mantid::API::equals() function
    TS_ASSERT( (!Mantid::API::equals(ws1, ws2)) );
  }
  
  void testHistNotHist()
  {
    if ( !checker.isInitialized() ) checker.initialize();

    Mantid::API::MatrixWorkspace_sptr ws2 = WorkspaceCreationHelper::Create2DWorkspace123(2,2,true);
    
    TS_ASSERT_THROWS_NOTHING( checker.setProperty("Workspace1",ws1) );
    TS_ASSERT_THROWS_NOTHING( checker.setProperty("Workspace2",ws2) );
    
    TS_ASSERT( checker.execute() );
    TS_ASSERT_EQUALS( checker.getPropertyValue("Result"), "Histogram/point-like mismatch" );

    // Same, using the !Mantid::API::equals() function
    TS_ASSERT( (!Mantid::API::equals(ws1, ws2)) );
  }

  void testDistNonDist()
  {
    if ( !checker.isInitialized() ) checker.initialize();

    Mantid::API::MatrixWorkspace_sptr ws2 = WorkspaceCreationHelper::Create2DWorkspace123(2,2);
    ws2->isDistribution(true);
    
    TS_ASSERT_THROWS_NOTHING( checker.setProperty("Workspace1",ws1) );
    TS_ASSERT_THROWS_NOTHING( checker.setProperty("Workspace2",ws2) );
    
    TS_ASSERT( checker.execute() );
    TS_ASSERT_EQUALS( checker.getPropertyValue("Result"), "Distribution flag mismatch" );

    // Same, using the !Mantid::API::equals() function
    TS_ASSERT( (!Mantid::API::equals(ws1, ws2)) );
  }

  void testDifferentAxisType()
  {
    if ( !checker.isInitialized() ) checker.initialize();

    Mantid::API::MatrixWorkspace_sptr ws2 = WorkspaceCreationHelper::Create2DWorkspace123(2,2);
    Mantid::API::Axis* const newAxis = new Mantid::API::NumericAxis(2);
    ws2->replaceAxis(1,newAxis);
    
    TS_ASSERT_THROWS_NOTHING( checker.setProperty("Workspace1",ws1) );
    TS_ASSERT_THROWS_NOTHING( checker.setProperty("Workspace2",ws2) );
    
    TS_ASSERT( checker.execute() )
    TS_ASSERT_EQUALS( checker.getPropertyValue("Result"), "Axis 1 type mismatch" );

    // Same, using the !Mantid::API::equals() function
    TS_ASSERT( (!Mantid::API::equals(ws1, ws2)) );
  }

  void testDifferentAxisTitles()
  {
    if ( !checker.isInitialized() ) checker.initialize();

    Mantid::API::MatrixWorkspace_sptr ws2 = WorkspaceCreationHelper::Create2DWorkspace123(2,2);
    ws2->getAxis(0)->title() = "blah";
    
    TS_ASSERT_THROWS_NOTHING( checker.setProperty("Workspace1",ws1) );
    TS_ASSERT_THROWS_NOTHING( checker.setProperty("Workspace2",ws2) );
    
    TS_ASSERT( checker.execute() );
    TS_ASSERT_EQUALS( checker.getPropertyValue("Result"), "Axis 0 title mismatch" );

    // Same, using the !Mantid::API::equals() function
    TS_ASSERT( (!Mantid::API::equals(ws1, ws2)) );
  }

  void testDifferentAxisUnit()
  {
    if ( !checker.isInitialized() ) checker.initialize();

    Mantid::API::MatrixWorkspace_sptr ws2 = WorkspaceCreationHelper::Create2DWorkspace123(2,2);
    ws2->getAxis(0)->unit() = Mantid::Kernel::UnitFactory::Instance().create("Wavelength");
    
    TS_ASSERT_THROWS_NOTHING( checker.setProperty("Workspace1",ws1) );
    TS_ASSERT_THROWS_NOTHING( checker.setProperty("Workspace2",ws2) );
    
    TS_ASSERT( checker.execute() );
    TS_ASSERT_EQUALS( checker.getPropertyValue("Result"), "Axis 0 unit mismatch" );
    // Same, using the !Mantid::API::equals() function
    TS_ASSERT( (!Mantid::API::equals(ws1, ws2)) );
  }

  void testDifferentAxisValues()
  {
    if ( !checker.isInitialized() ) checker.initialize();

    Mantid::API::MatrixWorkspace_sptr ws2 = WorkspaceCreationHelper::Create2DWorkspace123(2,2);
    ws2->getAxis(1)->setValue(1,99);
    
    TS_ASSERT_THROWS_NOTHING( checker.setProperty("Workspace1",ws1) );
    TS_ASSERT_THROWS_NOTHING( checker.setProperty("Workspace2",ws2) );
    
    TS_ASSERT( checker.execute() )
    TS_ASSERT_EQUALS( checker.getPropertyValue("Result"), "Axis 1 values mismatch" ) ;
    // Same, using the !Mantid::API::equals() function
    TS_ASSERT( (!Mantid::API::equals(ws1, ws2)) );
  }

  void testDifferentYUnit()
  {
    if ( !checker.isInitialized() ) checker.initialize();

    Mantid::API::MatrixWorkspace_sptr ws2 = WorkspaceCreationHelper::Create2DWorkspace123(2,2);
    ws2->setYUnit("blah");
    
    TS_ASSERT_THROWS_NOTHING( checker.setProperty("Workspace1",ws1) );
    TS_ASSERT_THROWS_NOTHING( checker.setProperty("Workspace2",ws2) );
    
    TS_ASSERT( checker.execute() );
    TS_ASSERT_EQUALS( checker.getPropertyValue("Result"), "YUnit mismatch" );
    // Same, using the !Mantid::API::equals() function
    TS_ASSERT( (!Mantid::API::equals(ws1, ws2)) );
  }

  void testDifferentSpectraMap()
  {
    if ( !checker.isInitialized() ) checker.initialize();

    Mantid::API::MatrixWorkspace_sptr ws2 = WorkspaceCreationHelper::Create2DWorkspace123(2,2);
    ws2->getSpectrum(0)->setSpectrumNo(1234);
    TS_ASSERT_THROWS_NOTHING( checker.setProperty("Workspace1",ws1) );
    TS_ASSERT_THROWS_NOTHING( checker.setProperty("Workspace2",ws2) );
    TS_ASSERT( checker.execute() );
    TS_ASSERT_EQUALS( checker.getPropertyValue("Result"), "Spectrum number mismatch" );

    ws2 = WorkspaceCreationHelper::Create2DWorkspace123(2,2);
    ws2->getSpectrum(0)->setDetectorID(99);
    ws2->getSpectrum(1)->setDetectorID(98);
    TS_ASSERT_THROWS_NOTHING( checker.setProperty("Workspace1",ws1) );
    TS_ASSERT_THROWS_NOTHING( checker.setProperty("Workspace2",ws2) );
    TS_ASSERT( checker.execute() );
    TS_ASSERT_EQUALS( checker.getPropertyValue("Result"), "Detector IDs mismatch" );

    // Same, using the !Mantid::API::equals() function
    TS_ASSERT( (!Mantid::API::equals(ws1, ws2)) );
  }

  void testDifferentInstruments()
  {
    if ( !checker.isInitialized() ) checker.initialize();

    Mantid::API::MatrixWorkspace_sptr ws2 = WorkspaceCreationHelper::Create2DWorkspace123(2,2);
    Mantid::Geometry::Instrument_sptr instrument(new Mantid::Geometry::Instrument("different"));
    ws2->setInstrument(instrument);
    
    TS_ASSERT_THROWS_NOTHING( checker.setProperty("Workspace1",ws1) );
    TS_ASSERT_THROWS_NOTHING( checker.setProperty("Workspace2",ws2) );
    
    TS_ASSERT( checker.execute() );
    TS_ASSERT_EQUALS( checker.getPropertyValue("Result"), "Instrument name mismatch" );
    // Same, using the !Mantid::API::equals() function
    TS_ASSERT( (!Mantid::API::equals(ws1, ws2)) );
  }

  void testDifferentParameterMaps()
  {
    if ( !checker.isInitialized() ) checker.initialize();

    Mantid::API::MatrixWorkspace_sptr ws2 = WorkspaceCreationHelper::Create2DWorkspace123(2,2);
    ws2->instrumentParameters().addBool(new Mantid::Geometry::Component, "myParam", true);
    
    TS_ASSERT_THROWS_NOTHING( checker.setProperty("Workspace1",ws1) );
    TS_ASSERT_THROWS_NOTHING( checker.setProperty("Workspace2",ws2) );
    
    TS_ASSERT( checker.execute() );
    TS_ASSERT_EQUALS( checker.getPropertyValue("Result"), "Instrument ParameterMap mismatch" );
    // Same, using the !Mantid::API::equals() function
    TS_ASSERT( (!Mantid::API::equals(ws1, ws2)) );
  }

  void testDifferentMasking()
  {
    if ( !checker.isInitialized() ) checker.initialize();

    Mantid::API::MatrixWorkspace_sptr ws2 = WorkspaceCreationHelper::Create2DWorkspace123(2,2);
    ws2->maskBin(0,0);
    ws2->dataY(0)[0] = 2;
    ws2->dataE(0)[0] = 3;
    
    TS_ASSERT_THROWS_NOTHING( checker.setProperty("Workspace1",ws1) );
    TS_ASSERT_THROWS_NOTHING( checker.setProperty("Workspace2",ws2) );
    
    TS_ASSERT( checker.execute() );
    TS_ASSERT_EQUALS( checker.getPropertyValue("Result"), "Masking mismatch" );
    
    Mantid::API::MatrixWorkspace_sptr ws3 = WorkspaceCreationHelper::Create2DWorkspace123(2,2);
    ws3->maskBin(0,1);
    ws3->dataY(0)[1] = 2;
    ws3->dataE(0)[1] = 3;

    TS_ASSERT_THROWS_NOTHING( checker.setProperty("Workspace1",ws3) );
    TS_ASSERT_THROWS_NOTHING( checker.setProperty("Workspace2",ws2) );
    
    TS_ASSERT( checker.execute() );
    TS_ASSERT_EQUALS( checker.getPropertyValue("Result"), "Masking mismatch" );
    // Same, using the !Mantid::API::equals() function
    TS_ASSERT( (!Mantid::API::equals(ws1, ws2)) );
  }
  
  void testDifferentSampleName()
  {
    if ( !checker.isInitialized() ) checker.initialize();
    checker.setProperty("CheckSample",true);

    Mantid::API::MatrixWorkspace_sptr ws2 = WorkspaceCreationHelper::Create2DWorkspace123(2,2);
    ws2->mutableSample().setName("different");
    
    TS_ASSERT_THROWS_NOTHING( checker.setProperty("Workspace1",ws1) );
    TS_ASSERT_THROWS_NOTHING( checker.setProperty("Workspace2",ws2) );
    
    TS_ASSERT( checker.execute() );
    TS_ASSERT_EQUALS( checker.getPropertyValue("Result"), "Sample name mismatch" );
  }

  void testDifferentProtonCharge()
  {
    if ( !checker.isInitialized() ) checker.initialize();
    checker.setProperty("CheckSample",true);

    Mantid::API::MatrixWorkspace_sptr ws2 = WorkspaceCreationHelper::Create2DWorkspace123(2,2);
    ws2->mutableRun().setProtonCharge(99.99);
    
    TS_ASSERT_THROWS_NOTHING( checker.setProperty("Workspace1",ws1) );
    TS_ASSERT_THROWS_NOTHING( checker.setProperty("Workspace2",ws2) );
    
    TS_ASSERT( checker.execute() );
    TS_ASSERT_EQUALS( checker.getPropertyValue("Result"), "Proton charge mismatch" );
  }

  void testDifferentLogs()
  {
    if ( !checker.isInitialized() ) checker.initialize();
    checker.setProperty("CheckSample",true);
    
    Mantid::API::MatrixWorkspace_sptr ws2 = WorkspaceCreationHelper::Create2DWorkspace123(2,2);
    ws2->mutableRun().addLogData(new Mantid::Kernel::PropertyWithValue<int>("Prop1",99));
    
    TS_ASSERT_THROWS_NOTHING( checker.setProperty("Workspace1",ws1) );
    TS_ASSERT_THROWS_NOTHING( checker.setProperty("Workspace2",ws2) );
    
    TS_ASSERT( checker.execute() );
    TS_ASSERT_EQUALS( checker.getPropertyValue("Result"), "Different numbers of logs" );
    
    Mantid::API::MatrixWorkspace_sptr ws3 = WorkspaceCreationHelper::Create2DWorkspace123(2,2);
    ws3->mutableRun().addLogData(new Mantid::Kernel::PropertyWithValue<int>("Prop2",99));
    
    TS_ASSERT_THROWS_NOTHING( checker.setProperty("Workspace1",ws2) );
    TS_ASSERT_THROWS_NOTHING( checker.setProperty("Workspace2",ws3) );
    
    TS_ASSERT( checker.execute() );
    TS_ASSERT_EQUALS( checker.getPropertyValue("Result"), "Log name mismatch" );
    
    Mantid::API::MatrixWorkspace_sptr ws4 = WorkspaceCreationHelper::Create2DWorkspace123(2,2);
    ws4->mutableRun().addLogData(new Mantid::Kernel::PropertyWithValue<int>("Prop1",100));
    
    TS_ASSERT_THROWS_NOTHING( checker.setProperty("Workspace1",ws2) );
    TS_ASSERT_THROWS_NOTHING( checker.setProperty("Workspace2",ws4) );
    
    TS_ASSERT( checker.execute() );
    TS_ASSERT_EQUALS( checker.getPropertyValue("Result"), "Log value mismatch" );
  }

  void test_Input_With_Two_Groups_That_Are_The_Same_Matches()
  {
    // Create a group
    const std::string groupName("TestGroup");
    WorkspaceGroup_sptr group = WorkspaceCreationHelper::CreateWorkspaceGroup(2, 2, 2, groupName);

    doGroupTest(groupName, groupName, Mantid::Algorithms::CheckWorkspacesMatch::successString());

    cleanupGroup(group);
  }

  void test_Input_With_Two_Groups_That_Are_Different_Sizes_Fails()
  {
    // Create a group
    const std::string groupOneName("TestGroupOne");
    WorkspaceGroup_sptr groupOne = WorkspaceCreationHelper::CreateWorkspaceGroup(2, 2, 2, groupOneName);
    const std::string groupTwoName("TestGroupTwo");
    WorkspaceGroup_sptr groupTwo = WorkspaceCreationHelper::CreateWorkspaceGroup(3, 2, 2, groupTwoName);

    doGroupTest(groupOneName, groupTwoName, "GroupWorkspaces size mismatch.", std::map<std::string,std::string>(), true);

    cleanupGroup(groupOne);
    cleanupGroup(groupTwo);
  }

  void test_Input_With_A_Group_And_A_Single_Workspace_Gives_Type_Mismatch()
  {
    const std::string groupName("CheckWorkspacesMatch_TestGroup");
    WorkspaceGroup_sptr group = WorkspaceCreationHelper::CreateWorkspaceGroup(2, 2, 2, groupName);
    Mantid::API::MatrixWorkspace_sptr ws2 = WorkspaceCreationHelper::Create2DWorkspace123(2,2);
    const std::string wsName("CheckWorkspacesMatch_TestWS");
    Mantid::API::AnalysisDataService::Instance().add(wsName,ws2);

    doGroupTest(groupName, wsName, "Type mismatch. One workspace is a group, the other is not.");

    // Cleanup
    cleanupGroup(group);
    Mantid::API::AnalysisDataService::Instance().remove(wsName);
  }

  void test_Input_With_Two_Groups_When_Single_Item_Checking_Is_Disabled()
  {
    // Create a group
    const std::string groupOneName("TestGroupOne");
    WorkspaceGroup_sptr groupOne = WorkspaceCreationHelper::CreateWorkspaceGroup(2, 2, 2, groupOneName);
    const std::string groupTwoName("TestGroupTwo");
    WorkspaceGroup_sptr groupTwo = WorkspaceCreationHelper::CreateWorkspaceGroup(2, 2, 2, groupTwoName);
    Mantid::API::AnalysisDataServiceImpl& dataStore =  Mantid::API::AnalysisDataService::Instance();
    // Extract the zeroth element of groupTwo and add a spurious log
    MatrixWorkspace_sptr zero = boost::dynamic_pointer_cast<MatrixWorkspace>(dataStore.retrieve(groupTwo->getNames()[0]));
    TS_ASSERT(zero);
    using Mantid::Kernel::PropertyWithValue;
    zero->mutableRun().addProperty(new PropertyWithValue<double>("ExtraLog", 10));

    std::map<std::string,std::string> otherProps;
    otherProps.insert(std::make_pair("CheckSample", "1"));

    doGroupTest(groupOneName, groupTwoName, "Different numbers of logs. Inputs=[TestGroupOne_0,TestGroupTwo_0]", otherProps);

    // Cleanup
    cleanupGroup(groupOne);
    cleanupGroup(groupTwo);
  }


private:

  void doGroupTest(const std::string & inputWSOne, const std::string & inputWSTwo,
                   const std::string & expectedResult,
                   const std::map<std::string,std::string> & otherProps = std::map<std::string,std::string>(),
                   bool expectFail = false
                   )
  {
    Mantid::Algorithms::CheckWorkspacesMatch matcher;
    matcher.initialize();
    matcher.setPropertyValue("Workspace1", inputWSOne);
    matcher.setPropertyValue("Workspace2", inputWSTwo);
    std::map<std::string,std::string>::const_iterator iend = otherProps.end();
    std::map<std::string,std::string>::const_iterator itr = otherProps.begin();
    for(; itr != iend; ++itr)
    {
      matcher.setPropertyValue(itr->first, itr->second);
    }

    TS_ASSERT_THROWS_NOTHING(matcher.execute());
    if (expectFail)
    {
      TS_ASSERT_EQUALS(matcher.isExecuted(), false);
      return;
    }
    TS_ASSERT_EQUALS(matcher.isExecuted(), true);
    TS_ASSERT_EQUALS(matcher.getPropertyValue("Result"), expectedResult);
  }

  void cleanupGroup(const WorkspaceGroup_sptr group)
  {
    group->deepRemoveAll();
    const std::string name = group->getName();
    Mantid::API::AnalysisDataService::Instance().remove(name);
  }

private:
  Mantid::Algorithms::CheckWorkspacesMatch checker;
  const Mantid::API::MatrixWorkspace_sptr ws1;
};

#endif /*CHECKWORKSPACESMATCHTEST_H_*/
