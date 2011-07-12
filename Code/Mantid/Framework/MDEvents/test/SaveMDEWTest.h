#ifndef MANTID_MDEVENTS_SAVEMDEWTEST_H_
#define MANTID_MDEVENTS_SAVEMDEWTEST_H_

#include "MantidKernel/System.h"
#include "MantidKernel/Timer.h"
#include "MantidMDEvents/MDEventFactory.h"
#include "MantidMDEvents/SaveMDEW.h"
#include "MantidTestHelpers/MDEventsTestHelper.h"
#include <cxxtest/TestSuite.h>
#include <iomanip>
#include <iostream>
#include "MantidKernel/CPUTimer.h"
#include "MantidTestHelpers/AlgorithmHelper.h"

using namespace Mantid::MDEvents;
using namespace Mantid::API;
using Mantid::Kernel::CPUTimer;

class SaveMDEWTest : public CxxTest::TestSuite
{
public:

    
  void test_Init()
  {
    SaveMDEW alg;
    TS_ASSERT_THROWS_NOTHING( alg.initialize() )
    TS_ASSERT( alg.isInitialized() )
  }
  
  void test_exec()
  {
    // Make a 1D MDEventWorkspace
    MDEventWorkspace1::sptr ws = MDEventsTestHelper::makeMDEW<1>(10, 0.0, 10.0, 23);

    // Recurse split so that it has lots more boxes, recursively
    MDEventsTestHelper::recurseSplit<1>( dynamic_cast<MDGridBox<MDEvent<1>,1>*>(ws->getBox()), 0, 4);

    // Add some points
    MDEventsTestHelper::feedMDBox(ws->getBox(), 1, 9e3, 1e-3);

    AnalysisDataService::Instance().addOrReplace("SaveMDEWTest_ws", ws);

//    std::ostringstream mess;
//    mess << num << ", " << x << ", " << y << ", " << z << ", " << radius;
//    AlgorithmHelper::runAlgorithm("FakeMDEventData", 4,
//        "InputWorkspace", "SaveMDEWTest_ws", "UniformParams", "1000000");


    ws->refreshCache();

    // There are this many boxes, so this is the max ID.
    TS_ASSERT_EQUALS( ws->getBoxController()->getMaxId(), 111111);

    IMDEventWorkspace_sptr iws = ws;

    CPUTimer tim;

    SaveMDEW alg;
    TS_ASSERT_THROWS_NOTHING( alg.initialize() )
    TS_ASSERT( alg.isInitialized() )
    TS_ASSERT_THROWS_NOTHING( alg.setPropertyValue("InputWorkspace", "SaveMDEWTest_ws") );
    TS_ASSERT_THROWS_NOTHING( alg.setPropertyValue("Filename", "SaveMDEWTest.nxs") );
    TS_ASSERT_THROWS_NOTHING( alg.execute(); );
    TS_ASSERT( alg.isExecuted() );

    std::cout << tim << " to save " << ws->getBoxController()->getMaxId() << " boxes." << std::endl;
  }
  

};


#endif /* MANTID_MDEVENTS_SAVEMDEWTEST_H_ */

