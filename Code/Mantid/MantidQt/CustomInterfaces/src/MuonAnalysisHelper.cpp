//----------------------------------------------------------------------
// Includes
//----------------------------------------------------------------------
#include "MantidQtCustomInterfaces/MuonAnalysis.h"
#include "MantidQtCustomInterfaces/MuonAnalysisHelper.h"
#include "MantidKernel/ConfigService.h"

#include "MantidQtAPI/UserSubWindow.h"

#include <boost/shared_ptr.hpp>
#include <fstream>  
//-----------------------------------------------------------------------------


namespace MantidQt
{
namespace CustomInterfaces
{
namespace Muon
{
  using namespace MantidQt::API;
  using namespace Mantid::Kernel;


/**
 * Add Greek letter to label from code 
 *
 */
void createMicroSecondsLabels(Ui::MuonAnalysis& m_uiForm)
{
  //the unicode code for the mu symbol is 956, doing the below keeps this file ASCII compatible
  static const QChar MU_SYM(956);
  m_uiForm.Time_Zero_label->setText(m_uiForm.Time_Zero_label->text() + QString(" (%1s)").arg(MU_SYM));
  m_uiForm.First_Good_Data_label->setText(m_uiForm.First_Good_Data_label->text() + QString(" (%1s)").arg(MU_SYM));
  m_uiForm.timeAxisStartAtLabel->setText(m_uiForm.timeAxisStartAtLabel->text() + QString(" (%1s)").arg(MU_SYM));
  m_uiForm.timeAxisFinishAtLabel->setText(m_uiForm.timeAxisFinishAtLabel->text() + QString(" (%1s)").arg(MU_SYM));
}


/**
 * Auto save various GUI selected items
 */
void autoSave(Ui::MuonAnalysis& m_uiForm)
{

  QSettings prevValues;
  prevValues.beginGroup("CustomInterfaces/MuonAnalysis/SaveOutput");

  // Get value for "dir". If the setting doesn't exist then use
  // the the path in "defaultsave.directory"
  QString prevPath = prevValues.value("dir", QString::fromStdString(
    ConfigService::Instance().getString("defaultsave.directory"))).toString();

  /*QString filter;
  filter.append("Files (*.XML *.xml)");
  filter += ";;AllFiles (*.*)";
  QString groupingFile = API::FileDialogHandler::getSaveFileName(this,
                                   "Save Grouping file as", prevPath, filter);

  if( ! groupingFile.isEmpty() )
  {
    saveGroupingTabletoXML(m_uiForm, groupingFile.toStdString());
    
    QString directory = QFileInfo(groupingFile).path();
    prevValues.setValue("dir", directory);
  }
*/
}

}
}
}
