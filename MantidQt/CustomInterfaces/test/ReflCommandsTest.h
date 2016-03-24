#ifndef MANTID_CUSTOMINTERFACES_REFLCOMMANDSTEST_H
#define MANTID_CUSTOMINTERFACES_REFLCOMMANDSTEST_H

#include <cxxtest/TestSuite.h>
#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "MantidAPI/FrameworkManager.h"
#include "MantidQtCustomInterfaces/Reflectometry/IReflTablePresenter.h"
#include "MantidQtCustomInterfaces/Reflectometry/ReflAppendRowCommand.h"
#include "MantidQtCustomInterfaces/Reflectometry/ReflClearSelectedCommand.h"
#include "MantidQtCustomInterfaces/Reflectometry/ReflCopySelectedCommand.h"
#include "MantidQtCustomInterfaces/Reflectometry/ReflCutSelectedCommand.h"
#include "MantidQtCustomInterfaces/Reflectometry/ReflDeleteRowCommand.h"
#include "MantidQtCustomInterfaces/Reflectometry/ReflExpandCommand.h"
#include "MantidQtCustomInterfaces/Reflectometry/ReflExportTableCommand.h"
#include "MantidQtCustomInterfaces/Reflectometry/ReflGroupRowsCommand.h"
#include "MantidQtCustomInterfaces/Reflectometry/ReflImportTableCommand.h"
#include "MantidQtCustomInterfaces/Reflectometry/ReflNewTableCommand.h"
#include "MantidQtCustomInterfaces/Reflectometry/ReflOpenTableCommand.h"
#include "MantidQtCustomInterfaces/Reflectometry/ReflOptionsCommand.h"
#include "MantidQtCustomInterfaces/Reflectometry/ReflPasteSelectedCommand.h"
#include "MantidQtCustomInterfaces/Reflectometry/ReflPlotGroupCommand.h"
#include "MantidQtCustomInterfaces/Reflectometry/ReflPlotRowCommand.h"
#include "MantidQtCustomInterfaces/Reflectometry/ReflPrependRowCommand.h"
#include "MantidQtCustomInterfaces/Reflectometry/ReflProcessCommand.h"
#include "MantidQtCustomInterfaces/Reflectometry/ReflSaveTableAsCommand.h"
#include "MantidQtCustomInterfaces/Reflectometry/ReflSaveTableCommand.h"
#include "MantidQtCustomInterfaces/Reflectometry/ReflTableViewPresenter.h"

using namespace MantidQt::CustomInterfaces;
using namespace Mantid::API;
using namespace testing;

class MockIReflTablePresenter : public IReflTablePresenter {

public:
  MOCK_METHOD1(notify, void(IReflTablePresenter::Flag));

private:
  // Calls we don't care about
  const std::map<std::string, QVariant> &options() const { return m_options; };
  std::vector<ReflCommandBase_uptr> publishCommands() {
    return std::vector<ReflCommandBase_uptr>();
  };
  std::set<std::string> getTableList() const {
    return std::set<std::string>();
  };
  // Calls we don't care about
  void setOptions(const std::map<std::string, QVariant> &options){};
  void transfer(const std::vector<std::map<std::string, std::string>> &runs){};
  void setInstrumentList(const std::vector<std::string> &instruments,
                         const std::string &defaultInstrument){};
  void accept(WorkspaceReceiver *workspaceReceiver){};
  void setModel(std::string name){};

  std::map<std::string, QVariant> m_options;
};

//=====================================================================================
// Functional tests
//=====================================================================================
class ReflCommandsTest : public CxxTest::TestSuite {

private:
public:
  // This pair of boilerplate methods prevent the suite being created statically
  // This means the constructor isn't called when running other tests
  static ReflCommandsTest *createSuite() { return new ReflCommandsTest(); }
  static void destroySuite(ReflCommandsTest *suite) { delete suite; }

  ReflCommandsTest() { FrameworkManager::Instance(); }

  void test_open_table_command() {
    NiceMock<MockIReflTablePresenter> mockPresenter;
    ReflOpenTableCommand command(&mockPresenter);

    // The presenter should be notified with the OpenTableFlag
    EXPECT_CALL(mockPresenter, notify(IReflTablePresenter::OpenTableFlag))
        .Times(Exactly(1));
    // Execute the command
    command.execute();
    // Verify expectations
    TS_ASSERT(Mock::VerifyAndClearExpectations(&mockPresenter));
  }

  void test_new_table_command() {
    NiceMock<MockIReflTablePresenter> mockPresenter;
    ReflNewTableCommand command(&mockPresenter);

    // The presenter should be notified with the NewTableFlag
    EXPECT_CALL(mockPresenter, notify(IReflTablePresenter::NewTableFlag))
        .Times(Exactly(1));
    // Execute the command
    command.execute();
    // Verify expectations
    TS_ASSERT(Mock::VerifyAndClearExpectations(&mockPresenter));
  }

  void test_save_table_command() {
    NiceMock<MockIReflTablePresenter> mockPresenter;
    ReflSaveTableCommand command(&mockPresenter);

    // The presenter should be notified with the SaveFlag
    EXPECT_CALL(mockPresenter, notify(IReflTablePresenter::SaveFlag))
        .Times(Exactly(1));
    // Execute the command
    command.execute();
    // Verify expectations
    TS_ASSERT(Mock::VerifyAndClearExpectations(&mockPresenter));
  }

  void test_save_table_as_command() {
    NiceMock<MockIReflTablePresenter> mockPresenter;
    ReflSaveTableAsCommand command(&mockPresenter);

    // The presenter should be notified with the SaveAsFlag
    EXPECT_CALL(mockPresenter, notify(IReflTablePresenter::SaveAsFlag))
        .Times(Exactly(1));
    // Execute the command
    command.execute();
    // Verify expectations
    TS_ASSERT(Mock::VerifyAndClearExpectations(&mockPresenter));
  }

  void test_import_table_command() {
    NiceMock<MockIReflTablePresenter> mockPresenter;
    ReflImportTableCommand command(&mockPresenter);

    // The presenter should be notified with the ImportTableFlag
    EXPECT_CALL(mockPresenter, notify(IReflTablePresenter::ImportTableFlag))
        .Times(Exactly(1));
    // Execute the command
    command.execute();
    // Verify expectations
    TS_ASSERT(Mock::VerifyAndClearExpectations(&mockPresenter));
  }

  void test_export_table_command() {
    NiceMock<MockIReflTablePresenter> mockPresenter;
    ReflExportTableCommand command(&mockPresenter);

    // The presenter should be notified with the ExportTableFlag
    EXPECT_CALL(mockPresenter, notify(IReflTablePresenter::ExportTableFlag))
        .Times(Exactly(1));
    // Execute the command
    command.execute();
    // Verify expectations
    TS_ASSERT(Mock::VerifyAndClearExpectations(&mockPresenter));
  }

  void test_options_command() {
    NiceMock<MockIReflTablePresenter> mockPresenter;
    ReflOptionsCommand command(&mockPresenter);

    // The presenter should be notified with the OptionsDialogFlag
    EXPECT_CALL(mockPresenter, notify(IReflTablePresenter::OptionsDialogFlag))
        .Times(Exactly(1));
    // Execute the command
    command.execute();
    // Verify expectations
    TS_ASSERT(Mock::VerifyAndClearExpectations(&mockPresenter));
  }

  void test_process_command() {
    NiceMock<MockIReflTablePresenter> mockPresenter;
    ReflProcessCommand command(&mockPresenter);

    // The presenter should be notified with the ProcessFlag
    EXPECT_CALL(mockPresenter, notify(IReflTablePresenter::ProcessFlag))
        .Times(Exactly(1));
    // Execute the command
    command.execute();
    // Verify expectations
    TS_ASSERT(Mock::VerifyAndClearExpectations(&mockPresenter));
  }

  void test_expand_command() {
    NiceMock<MockIReflTablePresenter> mockPresenter;
    ReflExpandCommand command(&mockPresenter);

    // The presenter should be notified with the ExpandSelectionFlag
    EXPECT_CALL(mockPresenter, notify(IReflTablePresenter::ExpandSelectionFlag))
        .Times(Exactly(1));
    // Execute the command
    command.execute();
    // Verify expectations
    TS_ASSERT(Mock::VerifyAndClearExpectations(&mockPresenter));
  }

  void test_plot_row_command() {
    NiceMock<MockIReflTablePresenter> mockPresenter;
    ReflPlotRowCommand command(&mockPresenter);

    // The presenter should be notified with the PlotRowFlag
    EXPECT_CALL(mockPresenter, notify(IReflTablePresenter::PlotRowFlag))
        .Times(Exactly(1));
    // Execute the command
    command.execute();
    // Verify expectations
    TS_ASSERT(Mock::VerifyAndClearExpectations(&mockPresenter));
  }

  void test_plot_group_command() {
    NiceMock<MockIReflTablePresenter> mockPresenter;
    ReflPlotGroupCommand command(&mockPresenter);

    // The presenter should be notified with the PlotGroupFlag
    EXPECT_CALL(mockPresenter, notify(IReflTablePresenter::PlotGroupFlag))
        .Times(Exactly(1));
    // Execute the command
    command.execute();
    // Verify expectations
    TS_ASSERT(Mock::VerifyAndClearExpectations(&mockPresenter));
  }

  void test_append_row_command() {
    NiceMock<MockIReflTablePresenter> mockPresenter;
    ReflAppendRowCommand command(&mockPresenter);

    // The presenter should be notified with the AppendRowFlag
    EXPECT_CALL(mockPresenter, notify(IReflTablePresenter::AppendRowFlag))
        .Times(Exactly(1));
    // Execute the command
    command.execute();
    // Verify expectations
    TS_ASSERT(Mock::VerifyAndClearExpectations(&mockPresenter));
  }

  void test_prepend_row_command() {
    NiceMock<MockIReflTablePresenter> mockPresenter;
    ReflPrependRowCommand command(&mockPresenter);

    // The presenter should be notified with the PrependRowFlag
    EXPECT_CALL(mockPresenter, notify(IReflTablePresenter::PrependRowFlag))
        .Times(Exactly(1));
    // Execute the command
    command.execute();
    // Verify expectations
    TS_ASSERT(Mock::VerifyAndClearExpectations(&mockPresenter));
  }

  void test_group_rows_command() {
    NiceMock<MockIReflTablePresenter> mockPresenter;
    ReflGroupRowsCommand command(&mockPresenter);

    // The presenter should be notified with the GroupRowsFlag
    EXPECT_CALL(mockPresenter, notify(IReflTablePresenter::GroupRowsFlag))
        .Times(Exactly(1));
    // Execute the command
    command.execute();
    // Verify expectations
    TS_ASSERT(Mock::VerifyAndClearExpectations(&mockPresenter));
  }

  void test_copy_selected_command() {
    NiceMock<MockIReflTablePresenter> mockPresenter;
    ReflCopySelectedCommand command(&mockPresenter);

    // The presenter should be notified with the CopySelectedFlag
    EXPECT_CALL(mockPresenter, notify(IReflTablePresenter::CopySelectedFlag))
        .Times(Exactly(1));
    // Execute the command
    command.execute();
    // Verify expectations
    TS_ASSERT(Mock::VerifyAndClearExpectations(&mockPresenter));
  }

  void test_cut_selected_command() {
    NiceMock<MockIReflTablePresenter> mockPresenter;
    ReflCutSelectedCommand command(&mockPresenter);

    // The presenter should be notified with the CutSelectedFlag
    EXPECT_CALL(mockPresenter, notify(IReflTablePresenter::CutSelectedFlag))
        .Times(Exactly(1));
    // Execute the command
    command.execute();
    // Verify expectations
    TS_ASSERT(Mock::VerifyAndClearExpectations(&mockPresenter));
  }

  void test_paste_selected_command() {
    NiceMock<MockIReflTablePresenter> mockPresenter;
    ReflPasteSelectedCommand command(&mockPresenter);

    // The presenter should be notified with the PasteSelectedFlag
    EXPECT_CALL(mockPresenter, notify(IReflTablePresenter::PasteSelectedFlag))
        .Times(Exactly(1));
    // Execute the command
    command.execute();
    // Verify expectations
    TS_ASSERT(Mock::VerifyAndClearExpectations(&mockPresenter));
  }

  void test_clear_selected_command() {
    NiceMock<MockIReflTablePresenter> mockPresenter;
    ReflClearSelectedCommand command(&mockPresenter);

    // The presenter should be notified with the ClearSelectedFlag
    EXPECT_CALL(mockPresenter, notify(IReflTablePresenter::ClearSelectedFlag))
        .Times(Exactly(1));
    // Execute the command
    command.execute();
    // Verify expectations
    TS_ASSERT(Mock::VerifyAndClearExpectations(&mockPresenter));
  }

  void test_delete_row_command() {
    NiceMock<MockIReflTablePresenter> mockPresenter;
    ReflDeleteRowCommand command(&mockPresenter);

    // The presenter should be notified with the DeleteRowFlag
    EXPECT_CALL(mockPresenter, notify(IReflTablePresenter::DeleteRowFlag))
        .Times(Exactly(1));
    // Execute the command
    command.execute();
    // Verify expectations
    TS_ASSERT(Mock::VerifyAndClearExpectations(&mockPresenter));
  }
};
#endif /* MANTID_CUSTOMINTERFACES_REFLCOMMANDSTEST_H */
