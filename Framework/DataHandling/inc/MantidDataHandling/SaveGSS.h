#ifndef DATAHANDING_SAVEGSS_H_
#define DATAHANDING_SAVEGSS_H_

//---------------------------------------------------
// Includes
//---------------------------------------------------
#include "MantidAPI/Algorithm.h"
#include "MantidAPI/SpectrumInfo.h"
#include "MantidAPI/Run.h"
#include "MantidKernel/System.h"
#include "MantidKernel/cow_ptr.h"

#include <iosfwd>
#include <vector>

// Forward declare
namespace Mantid {
namespace HistogramData {
class Histogram;
}

namespace DataHandling {
/**
     Saves a focused data set
     into a three column GSAS format containing X_i, Y_i*step, and E_i*step.
   Exclusively for
     the crystallography package GSAS and data needs to be in time-of-flight
     For data where the focusing routine has generated several spectra (for
   example, multi-bank instruments),
     the option is provided for saving all spectra into a single file, separated
   by headers, or into
     several files that will be named "workspaceName_"+spectra_number.

     Required properties:
     <UL>
     <LI> InputWorkspace - The workspace name to save. </LI>
     <LI> Filename       - The filename for output </LI>
     </UL>

     Optional properties:
     <UL>
     <LI> SplitFiles - Option for splitting into N files for workspace with
   N-spectra</LI>
     <LI> Append     - Append to Filename, if it already exists (default:
   true).</LI>
     <LI> Bank       - The bank number of the first spectrum (default: 1)</LI>
     </UL>

     @author Laurent Chapon, ISIS Facility, Rutherford Appleton Laboratory
     @date 04/03/2009

     Copyright &copy; 2009-2010 ISIS Rutherford Appleton Laboratory, NScD Oak
   Ridge National Laboratory & European Spallation Source

     This file is part of Mantid.

     Mantid is free software; you can redistribute it and/or modify
     it under the terms of the GNU General Public License as published by
     the Free Software Foundation; either version 3 of the License, or
     (at your option) any later version.

     Mantid is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU General Public License for more details.

     You should have received a copy of the GNU General Public License
     along with this program.  If not, see <http://www.gnu.org/licenses/>.

     File change history is stored at: <https://github.com/mantidproject/mantid>
     Code Documentation is available at: <http://doxygen.mantidproject.org>
  */
class DLLExport SaveGSS : public Mantid::API::Algorithm {
public:
  /// Constructor
  SaveGSS();
  /// Algorithm's name
  const std::string name() const override { return "SaveGSS"; }
  /// Summary of algorithms purpose
  const std::string summary() const override {
    return "Saves a focused data set into a three column GSAS format.";
  }

  /// Algorithm's version
  int version() const override { return (1); }
  /// Algorithm's category for identification
  const std::string category() const override {
    return "Diffraction\\DataHandling;DataHandling\\Text";
  }

private:
  /// Initialisation code
  void init() override;
  /// Execution code
  void exec() override;

  /// Determines if all spectra have detectors
  bool areAllDetectorsValid() const;

  /// Turns the data associated with this spectra into a string stream
  void generateBankData(std::stringstream &outBuf, size_t specIndex) const;

  /// Generates the bank header and returns this as a string stream
  void generateBankHeader(std::stringstream &out,
                          const API::SpectrumInfo &spectrumInfo,
                          size_t specIndex) const;

  /// Generates the output which will be written to the GSAS file
  void generateGSASBuffer(size_t numOutFiles, size_t numOutSpectra);

  /// Generates the instrument header and returns this as a string stream
  void generateInstrumentHeader(std::stringstream &out, double l1) const;

  /// Generates the filename(s) and paths to write to and stores in member var
  void generateOutFileNames(size_t numberOfOutFiles);

  /// Returns the log value in a GSAS format as a string stream
  void getLogValue(std::stringstream &out, const API::Run &runInfo,
                   const std::string &name,
                   const std::string &failsafeValue = "UNKNOWN") const;

  /// Returns if the input workspace instrument is valid
  bool isInstrumentValid() const;

  /// Opens a new file stream at the path specified.
  std::ofstream openFileStream(const std::string &outFilePath);

  /// sets non workspace properties for the algorithm
  void setOtherProperties(IAlgorithm *alg, const std::string &propertyName,
                          const std::string &propertyValue,
                          int periodNum) override;

  /// Validates the user input and warns / throws on bad conditions
  void validateUserInput() const;

  /// Writes the current buffer to the user specified file path
  void writeBufferToFile(size_t numOutFiles, size_t numSpectra);

  // Writes the header for RALF data format to the buffer
  void writeRALFHeader(std::stringstream &out, int bank,
                       const HistogramData::Histogram &histo) const;

  /// Write out the data in RALF - ALT format
  void writeRALF_ALTdata(std::stringstream &out, const int bank,
                         const HistogramData::Histogram &histo) const;

  /// Write out the data in RALF - FXYE format
  void writeRALF_XYEdata(const int bank, const bool MultiplyByBinWidth,
                         std::stringstream &out,
                         const HistogramData::Histogram &histo) const;

  /// Write out the data in SLOG format
  void writeSLOGdata(const int bank, const bool MultiplyByBinWidth,
                     std::stringstream &out,
                     const HistogramData::Histogram &histo) const;

  /// Workspace
  API::MatrixWorkspace_const_sptr m_inputWS;
  /// The output buffer. This is either n spectra in one file,
  /// or n files with 1 spectra
  std::vector<std::stringstream> m_outputBuffer{};
  /// The output filename(s)
  std::vector<std::string> m_outFileNames{};
  /// Indicates whether all spectra have valid detectors
  bool m_allDetectorsValid{false};
  /// Holds pointer to progress bar
  std::unique_ptr<API::Progress> m_progress{nullptr};
};
}
}
#endif // DATAHANDING_SAVEGSS_H_
