#ifndef MANTIDQTCUSTOMINTERFACES_PROJECTSAVEMODELTEST_H
#define MANTIDQTCUSTOMINTERFACES_PROJECTSAVEMODELTEST_H

#include <cxxtest/TestSuite.h>

#include "MantidQtAPI/IProjectSerialisable.h"
#include "MantidQtCustomInterfaces/ProjectSaveModel.h"
#include "MantidTestHelpers/WorkspaceCreationHelper.h"
#include "ProjectSaveMockObjects.h"

using namespace MantidQt::API;
using namespace MantidQt::CustomInterfaces;

//=====================================================================================
// Functional tests
//=====================================================================================
class ProjectSaveModelTest : public CxxTest::TestSuite {

public:
  void setUp() override {
    auto ws1 = WorkspaceCreationHelper::Create1DWorkspaceRand(10);
    WorkspaceCreationHelper::storeWS("ws1", ws1);
    auto ws2 = WorkspaceCreationHelper::Create1DWorkspaceRand(10);
    WorkspaceCreationHelper::storeWS("ws2", ws2);
  }

  void tearDown() override {
    WorkspaceCreationHelper::removeWS("ws1");
    WorkspaceCreationHelper::removeWS("ws2");
  }

  void testConstructNoWorkspacesNoWindows() {
    tearDown(); // remove workspaces setup by default
    std::vector<MantidQt::API::IProjectSerialisable*> windows;
    TS_ASSERT_THROWS_NOTHING( ProjectSaveModel model(windows) );
  }

  void testConstructOneWorkspaceNoWindows() {
    std::vector<MantidQt::API::IProjectSerialisable*> windows;

    ProjectSaveModel model(windows);
    TS_ASSERT(!model.hasWindows("ws1"));
    TS_ASSERT_EQUALS(model.getWindows("ws1").size(), 0);
  }

  void testGetWindowsForWorkspaceNoWindows() {
    std::vector<MantidQt::API::IProjectSerialisable*> windows;

    ProjectSaveModel model(windows);
    TS_ASSERT(!model.hasWindows("ws1"));
    TS_ASSERT_EQUALS(model.getWindows("ws1").size(), 0);
  }

  void testGetWindowsForWorkspaceOneWindow() {
    std::vector<MantidQt::API::IProjectSerialisable*> windows;
    WindowStub win1("window1", {"ws1"});
    windows.push_back(&win1);

    ProjectSaveModel model(windows);
    TS_ASSERT(model.hasWindows("ws1"));
    TS_ASSERT_EQUALS(model.getWindows("ws1").size(), 1);
  }

  void testGetWindowsForWorkspaceTwoWindows() {
    std::vector<MantidQt::API::IProjectSerialisable*> windows;
    WindowStub win1("window1", {"ws1"});
    WindowStub win2("window2", {"ws1"});
    windows.push_back(&win1);
    windows.push_back(&win2);

    ProjectSaveModel model(windows);
    TS_ASSERT(model.hasWindows("ws1"));
    TS_ASSERT_EQUALS(model.getWindows("ws1").size(), 2);
  }

  void testGetWindowsForTwoWorkspacesAndTwoWindows() {
    std::vector<MantidQt::API::IProjectSerialisable*> windows;
    WindowStub win1("window1", {"ws1"});
    WindowStub win2("window2", {"ws2"});
    windows.push_back(&win1);
    windows.push_back(&win2);

    ProjectSaveModel model(windows);
    TS_ASSERT(model.hasWindows("ws1"));
    TS_ASSERT_EQUALS(model.getWindows("ws1").size(), 1);
    TS_ASSERT(model.hasWindows("ws2"));
    TS_ASSERT_EQUALS(model.getWindows("ws2").size(), 1);
  }

  void testGetWorkspaceNames() {
    std::vector<MantidQt::API::IProjectSerialisable*> windows;

    ProjectSaveModel model(windows);
    TS_ASSERT(!model.hasWindows("ws1"));
    TS_ASSERT(!model.hasWindows("ws2"));

    auto names = model.getWorkspaceNames();
    TS_ASSERT_EQUALS(names.size(), 2);
    TS_ASSERT_EQUALS(names[0], "ws1");
    TS_ASSERT_EQUALS(names[1], "ws2");
  }

};

#endif // MANTIDQTCUSTOMINTERFACES_PROJECTSAVEMODELTEST_H
