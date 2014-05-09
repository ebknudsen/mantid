#include "MantidCrystal/PeakClusterProjection.h"

#include "MantidKernel/V3D.h"
#include "MantidAPI/IPeak.h"
#include "MantidAPI/IMDHistoWorkspace.h"
#include "MantidAPI/PeakTransformHKL.h"
#include "MantidAPI/PeakTransformQLab.h"
#include "MantidAPI/PeakTransformQSample.h"

#include <stdexcept>

using namespace Mantid::Kernel;
using namespace Mantid::API;

namespace
{
  /**
   * Free function to generate a peak transform from a MDWS.
   *
   * @param mdWS
   * @return A peak transform matching the workspace configuration.
   */
  PeakTransform_sptr makePeakTransform(IMDHistoWorkspace const * const mdWS)
  {
    const SpecialCoordinateSystem mdCoordinates = mdWS->getSpecialCoordinateSystem();
    PeakTransformFactory_sptr peakTransformFactory;
    if (mdCoordinates == QLab)
    {
      peakTransformFactory = boost::make_shared<PeakTransformQLabFactory>();
    }
    else if (mdCoordinates == QSample)
    {
      peakTransformFactory = boost::make_shared<PeakTransformQSampleFactory>();
    }
    else if (mdCoordinates == Mantid::API::HKL)
    {
      peakTransformFactory = boost::make_shared<PeakTransformHKLFactory>();
    }
    const std::string xDim = mdWS->getDimension(0)->getName();
    const std::string yDim = mdWS->getDimension(1)->getName();
    PeakTransform_sptr peakTransform = peakTransformFactory->createTransform(xDim, yDim);
    return peakTransform;
  }
}



namespace Mantid
{
namespace Crystal
{


  //----------------------------------------------------------------------------------------------
  /** Constructor
   */
  PeakClusterProjection::PeakClusterProjection(IMDHistoWorkspace_sptr& mdWS) : m_mdWS(mdWS)
  {
    auto coordinateSystem = mdWS->getSpecialCoordinateSystem();
    if(coordinateSystem == None)
    {
      throw std::invalid_argument("Input image IMDHistoWorkspace must have a known special coordinate system.");
    }

    if(mdWS->getNumDims() < 3)
    {
      throw std::invalid_argument("Need to have 3 or more dimension in the workspace.");
    }

    // Make a peak transform so that we can understand a peak in the context of the mdworkspace coordinate setup.
    m_peakTransform = makePeakTransform(mdWS.get());
  }

  /**
   * Get the signal value at the peak center.
   * @param peak
   * @param normalization : Optional normalization. Off by default.
   * @return signal value at peak center. NAN if the peak is not centered on the image.
   */
  Mantid::signal_t PeakClusterProjection::signalAtPeakCenter(IPeak& peak, MDNormalization normalization) const
  {
    const Mantid::Kernel::V3D& center = m_peakTransform->transformPeak(peak);
    return m_mdWS->getSignalAtVMD(center, normalization);
  }
    
  //----------------------------------------------------------------------------------------------
  /** Destructor
   */
  PeakClusterProjection::~PeakClusterProjection()
  {
  }
  


} // namespace Crystal
} // namespace Mantid
