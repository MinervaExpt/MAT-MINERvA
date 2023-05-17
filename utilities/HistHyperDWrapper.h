//=============================================================================
/*!
  @brief A wrapper around MnvH1Ds which allows one to interface more easily,
         intuitively, and directly with error bands and their constituent
         universes.

  This class allows one to more easily and intuitively interface with an
  MnvH1D's systematic error bands and their constituent universes.

  It was built to eliminate the distinction between lateral and vertical error
  bands, eliminate FillVert/LatErrorBand functions, and standardize as many
  MINERvA errors as possible, and linearize the process of looping events and
  filling error bands.

  e.g. To create an MnvH1D and fill it with the full set of standard genie
  systematics:

      std::map<std::string, std::vector<CVUniverse*>> error_bands = GetGenieSystematics<CVUniverse>(my_event_chain);
      PlotUtils::HistHyperDWrapper<CVUniverse> my_hw("E_{#nu}", nbins, xmin, xmax, error_bands);
      loop my_event_chain
        loop band in error_bands
          loop universe in band
            if(PassesCuts(universe))
              my_hw.univHist(universe)->Fill(universe->GetEnu(), universe->GetWeight);

  Pass a container of universe objects, and the constructor will make an
  MnvH1D and populate it with MnvVertErrorBands.
  Then, in your loop over events and universes, HistHyperDWrapper::univHist accesses
  the correct TH1 so you can fill it directly.

  The user must write her own CVUniverse class, which can be very simple,
  containing only these members:

  - virtual std::string ShortName() const = 0;
  - virtual std::string LatexName() const = 0;
  - virtual double GetWeight() const;
  - PlotUtils::ChainWrapper& m_chw;
  - double m_nsigma;
  - Long64_t m_entry;

  See a complete example [here](http://cdcvs0.fnal.gov/cgi-bin/public-cvs/cvsweb-public.cgi/Personal/bmesserl/SystematicsFramework/?cvsroot=mnvsoft)

  @author Ben Messerly
*/
//=============================================================================

#ifndef HISTHYPERDWRAPPER_H
#define HISTHYPERDWRAPPER_H

#include <map>
#include <vector>

// #include "PlotUtils/MnvH1D.h"
// #include "PlotUtils/MnvH2D.h"
#include "PlotUtils/Hist2DWrapper.h"
#include "PlotUtils/HistWrapper.h"
#include "PlotUtils/HyperDimLinearizer.h"

namespace PlotUtils {
template <typename T>
struct HistHyperDWrapper {
    // Default Constructor
    HistHyperDWrapper();

    // 1D constructors
    //! 1D Constructor from a vector of universes
    HistHyperDWrapper(const char* hist_name, const char* title,
                      int nBins, double xmin, double xmax,
                      std::vector<T*>& univs,
                      EAnalysisType type = k1D);

    //! 1D Constructor from a map of errors/universes, uniform bins
    HistHyperDWrapper(const char* hist_name, const char* title,
                      int nBins, double xmin, double xmax,
                      std::map<std::string, std::vector<T*> >& bands,
                      EAnalysisType type = k1D);

    //! 1D Constructor from a pre-existing MnvH1D and a vector of universes
    HistHyperDWrapper(MnvH1D* h1d, std::vector<T*>& univs,
                      bool clear_error_bands = false);

    //! 1D Constructor from a pre-existing MnvH1D and a map of errors/universes
    HistHyperDWrapper(MnvH1D* h1d, std::map<std::string, std::vector<T*> >& bands,
                      bool clear_error_bands = false);

// The next 2 constructors are hidden from our python bindings because they use a c++11 feature: delegating constructors
#ifndef __GCCXML__
    //! 1D Constructor from a vector of universes and variable bin widths
    HistHyperDWrapper(const char* name, const char* title,
                      const std::vector<double> bins,
                      std::vector<T*>& univs,
                      EAnalysisType type = k1D);

    //! 1D Constructor from a map of universes and variable bin widths
    HistHyperDWrapper(const char* name, const char* title,
                      const std::vector<double> bins,
                      std::map<std::string, std::vector<T*> >& bands,
                      EAnalysisType type = k1D);

#endif  //__GCCXML__

    //! 1D Constructor from a histogram to add universes later, uniform bins. HMS
    HistHyperDWrapper(const char* hist_name, const char* title,
                      int nBins, double xmin, double xmax,
                      EAnalysisType type = k1D);

    //! 1D Constructor from a histogram to add universes later, variable bins. HMS
    HistHyperDWrapper(const char* hist_name, const char* title,
                      int nBins, std::vector<double> bins,
                      EAnalysisType type = k1D);

    // 2D Constructors

    //! 2D Constructor from a vector of universes
    HistHyperDWrapper(const char* hist_name, const char* title,
                      int nBinsX, double xmin, double xmax,
                      int nBinsY, double ymin, double ymax,
                      std::vector<T*>& univs,
                      EAnalysisType type = k2D);

    //! 2D Constructor from a map of errors/universes
    HistHyperDWrapper(const char* hist_name, const char* title,
                      int nBins, double xmin, double xmax,
                      int nBinsY, double ymin, double ymax,
                      std::map<std::string, std::vector<T*> >& bands,
                      EAnalysisType type = k2D);

    //! 2D Constructor from a pre-existing MnvH2D and a vector of universes
    HistHyperDWrapper(MnvH2D* h2d, std::vector<T*>& univs, bool clear_error_bands = false);

    //! 2D Constructor from a pre-existing MnvH2D and a map of errors/universes
    HistHyperDWrapper(MnvH2D* h2d, std::map<std::string, std::vector<T*> >& bands, bool clear_error_bands = false);

#ifndef __GCCXML__
    //! 2D Constructor from a vector of universes and variable bin widths
    HistHyperDWrapper(const char* name, const char* title,
                      const std::vector<double> xBins, const std::vector<double> yBins,
                      std::vector<T*>& univs,
                      EAnalysisType type = k2D);

    //! 2D Constructor from a map of universes and variable bin widths
    HistHyperDWrapper(const char* name, const char* title,
                      const std::vector<double> xBins, const std::vector<double> yBins,
                      std::map<std::string, std::vector<T*> >& bands,
                      EAnalysisType type = k2D);

#endif  // __GCCXML__

    // Data members

    PlotUtils::HistWrapper<T>* hist;
    PlotUtils::Hist2DWrapper<T>* hist2D;

    // PlotUtils::MnvH1D* hist;                        /*!< The MnvH1D that is created and filled */
    // PlotUtils::MnvH2D* hist2D;                        /*!< The MnvH2D that is created and filled */

    // std::map<const T*, TH1D*> univToHistMap;        /*!< The map between univs and TH2 */
    // std::map<const T*, TH2D*> univToHist2DMap;        /*!< The map between univs and TH2 */

    std::map<std::string, int> nhistsAssignedSoFar; /*!< Counter for nunivs in each band.*/

    EAnalysisType analysisType;

    //! Synchronize the MnvH2D's CV with each of its error band's CVs.
    void SyncCVHistos();

    void Write(TFile& file);

    //! Access universe TH21 given a universe object
    TH1D* univHist(const T* univ) const;

    //! Access universe TH2 given a universe object
    TH2D* univHist2D(const T* univ) const;

    //! Fill universe TH1
    void FillUniverse(const T& univ, double valueX, double weight);

    //! Fill universe TH1
    void FillUniverse(const T* univ, double valueX, double weight);

    //! Fill universe TH2
    void FillUniverse(const T& univ, double valueX, double valueY, double weight);

    //! Fill universe TH2
    void FillUniverse(const T* univ, double valueX, double valueY, double weight);

    void AddUniverses(const std::string name, const T* univ, const int nhists, const int histno);

   private:
    // Add vertical error bands to the MnvH2D and associate them with universes
    void FillErrorBandsWithSysUni(std::vector<T*> univs, int nhists, std::string name);
};
}  // namespace PlotUtils

// Template classes must have function definitions available in header.
#include "HistHyperDWrapper.cxx"

#endif  // HISTHYPERDWRAPPER_H
