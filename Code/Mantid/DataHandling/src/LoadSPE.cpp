//---------------------------------------------------
// Includes
//---------------------------------------------------
#include "MantidDataHandling/LoadSPE.h"
#include "MantidKernel/FileProperty.h"
#include "MantidKernel/UnitFactory.h"
#include "MantidDataObjects/Histogram1D.h"
#include <cstdio>

/// @cond
// Don't document this very long winded way of getting "degrees" to print on the axis.
namespace
{
  class Degrees : public Mantid::Kernel::Unit
  {
    const std::string unitID() const { return ""; }
    const std::string caption() const { return "Phi"; }
    const std::string label() const { return "degrees"; }
    void toTOF(std::vector<double>& xdata, std::vector<double>& ydata, const double& l1, const double& l2,
      const double& twoTheta, const int& emode, const double& efixed, const double& delta) const {}
    void fromTOF(std::vector<double>& xdata, std::vector<double>& ydata, const double& l1, const double& l2,
      const double& twoTheta, const int& emode, const double& efixed, const double& delta) const {}
  };
} // end anonynmous namespace
/// @endcond

namespace Mantid
{
namespace DataHandling
{

using namespace Kernel;
using namespace API;

// Register the algorithm into the AlgorithmFactory
DECLARE_ALGORITHM(LoadSPE)

//---------------------------------------------------
// Private member functions
//---------------------------------------------------

/**
 * Initialise the algorithm
 */
void LoadSPE::init()
{
  declareProperty(new FileProperty("Filename","", FileProperty::Load,
    std::vector<std::string>(1,"spe")),"Name of the SPE file to load" );
  declareProperty(new WorkspaceProperty<>("Workspace", "", Direction::Output),
    "The name to use for the output workspace" );
}

/**
 * Execute the algorithm
 */
void LoadSPE::exec()
{
  // Retrieve filename and try to open the file
  m_filename = getPropertyValue("Filename");

  FILE * speFile;
  speFile = fopen(m_filename.c_str(),"r");
  if ( !speFile )
  {
    g_log.error("Failed to open file: " + m_filename);
    throw Exception::FileError("Failed to open file: " , m_filename);
  }

  // The first two numbers are the number of histograms and the number of bins
  unsigned int nhist=0, nbins=0;
  int retval = fscanf(speFile,"%8u%8u\n",&nhist,&nbins);
  if ( retval != 2 ) reportFormatError("Header line");

  // Next line should be comment line: "### Phi Grid" or "### Q Grid"
  char comment[100];
  fgets(comment,100,speFile);
  if ( comment[0] != '#' ) reportFormatError(std::string(comment));

  // Create the axis that will hold the phi values
  Axis* phiAxis = new Axis(AxisType::Numeric,nhist+1);
  // Look at previously read comment field to see what unit vertical axis should have
  if ( comment[4] == 'Q' || comment[4] == 'q') 
  {
    phiAxis->unit() = UnitFactory::Instance().create("MomentumTransfer");
  }
  else 
  {
    phiAxis->unit() = boost::shared_ptr<Unit>(new Degrees);
  }

  // Read in phi grid
  for (unsigned int i = 0; i <= nhist; ++i)
  {
    double phi;
    retval = fscanf(speFile,"%10le",&phi);
    if ( retval != 1 ) 
    {
      std::stringstream ss;
      ss << "Reading phi value" << i;
      reportFormatError(ss.str());
    }
    phiAxis->setValue(i,phi);
  }
  // Read to EOL
  fgets(comment,100,speFile);

  // Next line should be comment line: "### Energy Grid"
  fgets(comment,100,speFile);
  if ( comment[0] != '#' ) reportFormatError(std::string(comment));

  // Now the X bin boundaries
  DataObjects::Histogram1D::RCtype XValues;
  MantidVec& X = XValues.access();
  X.resize(nbins+1);

  for (unsigned int i = 0; i <= nbins; ++i)
  {
    retval = fscanf(speFile,"%10le",&X[i]);
    if ( retval != 1 ) 
    {
      std::stringstream ss;
      ss << "Reading energy value" << i;
      reportFormatError(ss.str());
    }
  }
  // Read to EOL
  fgets(comment,100,speFile);

  // Now create the output workspace
  MatrixWorkspace_sptr workspace = WorkspaceFactory::Instance().create("Workspace2D",nhist,nbins+1,nbins);
  workspace->getAxis(0)->unit() = UnitFactory::Instance().create("DeltaE");
  workspace->isDistribution(true); // It should be a distribution
  workspace->setYUnitLabel("S(Phi,Energy)");
  // Replace the default spectrum axis with the phi values one
  workspace->replaceAxis(1,phiAxis);

  // Now read in the data spectrum-by-spectrum
  Progress progress(this,0,1,nhist);
  for (unsigned int j = 0; j < nhist; ++j) 
  {
    // Set the common X vector
    workspace->setX(j,XValues);
    // Read in the Y & E data
    readHistogram(speFile,workspace,j);

    progress.report();
  }

  // Close the file
  fclose(speFile);

  // Set the output workspace property
  setProperty("Workspace", workspace);
}

/** Reads in the data corresponding to a single spectrum
 *  @param speFile   The file handle
 *  @param workspace The output workspace
 *  @param index     The index of the current spectrum
 */
void LoadSPE::readHistogram(FILE* speFile, API::MatrixWorkspace_sptr workspace, int index)
{
  // First, there should be a comment line
  char comment[100];
  fgets(comment,100,speFile);
  if ( comment[0] != '#' ) reportFormatError(std::string(comment));

  // Then it's the Y values
  MantidVec& Y = workspace->dataY(index);
  const int nbins = workspace->blocksize();
  int retval;
  for (int i = 0; i < nbins; ++i)
  {
    retval = fscanf(speFile,"%10le",&Y[i]);
    //g_log.error() << Y[i] << std::endl;
    if ( retval != 1 ) 
    {
      std::stringstream ss;
      ss << "Reading data value" << i << " of histogram " << index;
      reportFormatError(ss.str());
    }
  }
  // Read to EOL
  fgets(comment,100,speFile);

  // Another comment line
  fgets(comment,100,speFile);
  if ( comment[0] != '#' ) reportFormatError(std::string(comment));

  // And then the error values
  MantidVec& E = workspace->dataE(index);
  for (int i = 0; i < nbins; ++i)
  {
    retval = fscanf(speFile,"%10le",&E[i]);
    if ( retval != 1 ) 
    {
      std::stringstream ss;
      ss << "Reading error value" << i << " of histogram " << index;
      reportFormatError(ss.str());
    }
  }
  // Read to EOL
  fgets(comment,100,speFile);

  return;
}

/** Called if the file is not formatted as expected
 *  @param what A string describing where the problem occurred
 *  @throws Mantid::Kernel::Exception::FileError terminating the algorithm
 */
void LoadSPE::reportFormatError(const std::string& what)
{
  g_log.error("Unexpected formatting in file " + m_filename + " : " + what);
  throw Exception::FileError("Unexpected formatting in file: " , m_filename);
}

} // namespace DataHandling
} // namespace Mantid
