#ifndef HISTHYPERDWRAPPER_CXX
#define HISTHYPERDWRAPPER_CXX

#include "PlotUtils/FluxSystematics.cxx"  // PlotUtils::flux_reweighter
#include "utilities/HistHyperDWrapper.h"

using namespace PlotUtils;


// Default Constructor
template <typename T>
HistHyperDWrapper<T>::HistHyperDWrapper()
    : hist(nullptr), hist2D(nullptr), nhistsAssignedSoFar() {}

// ============================================================================
// 1D Constructors
// ============================================================================
// TODO: Make sure member variables are set correctly
//! 1D Constructor from a vector of universes
template <typename T>
HistHyperDWrapper<T>::HistHyperDWrapper(const char* hist_name, const char* title,
                                        int nBins, double xmin, double xmax,
                                        std::vector<T*>& univs)
    : analysisType(k1D),
      hist(new HistWrapper<T>(hist_name, title, nBins, xmin, xmax, univs)),
      hist2D(nullptr) {}

//! 1D Constructor from a map of universes
template <typename T>
HistHyperDWrapper<T>::HistHyperDWrapper(const char* hist_name, const char* title,
                                        int nBins, double xmin, double xmax,
                                        std::map<std::string, std::vector<T*> >& bands)
    : analysisType(k1D),
      hist(new HistWrapper<T>(hist_name, title, nBins, xmin, xmax, bands)),
      hist2D(nullptr) {}

//! Constructor from a pre-existing MnvH1D and a vector of universes
template <typename T>
HistHyperDWrapper<T>::HistHyperDWrapper(MnvH1D* h1d,
                                        std::vector<T*>& univs,
                                        bool clear_error_bands)
    : analysisType(k1D),
      hist(new HistWrapper<T>(h1d, univs, clear_error_bands)),
      hist2D(nullptr) {}

// Constructor from a pre-existing MnvH1D & map< string, vector<universes> >
template <typename T>
HistHyperDWrapper<T>::HistHyperDWrapper(MnvH1D* h1d,
                                        std::map<std::string, std::vector<T*> >& bands,
                                        bool clear_error_bands)
    : analysisType(k1D),
      hist(new HistWrapper<T>(h1d, bands, clear_error_bands)),
      hist2D(nullptr) {}

#ifndef __GCCXML__
template <class T>
HistHyperDWrapper<T>::HistHyperDWrapper(const char* name, const char* title,
                                        const std::vector<double> bins,
                                        std::map<std::string, std::vector<T*> >& bands)
    : HistHyperDWrapper(
          new MnvH1D(name, title, bins.size() - 1, bins.data()),
          bands) {}

// Construct from a vector of universes and variable-sized bins.
template <class T>
HistHyperDWrapper<T>::HistHyperDWrapper(const char* name, const char* title,
                                        const std::vector<double> bins,
                                        std::vector<T*>& univs)
    : HistHyperDWrapper(
          new MnvH1D(name, title, bins.size() - 1, bins.data()),
          univs) {}
#endif  //__GCCXML__

template <class T>
HistHyperDWrapper<T>::HistHyperDWrapper(const char* name, const char* title,
                                        int nBins, double xmin, double xmax)
    : analysisType(k1D),
      hist(new HistWrapper<T>(name, title, nBins, xmin, xmax)),
      hist2D(nullptr) {}

template <class T>
HistHyperDWrapper<T>::HistHyperDWrapper(const char* name, const char* title,
                                        int nBins, std::vector<double> bins)
    : analysisType(k1D),
      hist(new HistWrapper<T>(name, title, nBins, bins)),
      hist2D(nullptr) {}

// ============================================================================
// 2D Constructors
// ============================================================================
// TODO: Check member variables assigned correctly

// Constructor from a vector of universes
template <typename T>
HistHyperDWrapper<T>::HistHyperDWrapper(const char* hist_name, const char* title,
                                        int nBinsX, double xmin, double xmax,
                                        int nBinsY, double ymin, double ymax,
                                        std::vector<T*>& univs)
    : analysisType(k2D),
      hist2D(new Hist2DWrapper<T>(hist_name, title, nBinsX, xmin, xmax, nBinsY, ymin, ymax, univs)),
      hist(nullptr) {}

// Constructor from map< string, vector<universe> >
template <typename T>
HistHyperDWrapper<T>::HistHyperDWrapper(const char* hist_name, const char* title,
                                        int nBinsX, double xmin, double xmax,
                                        int nBinsY, double ymin, double ymax,
                                        std::map<std::string, std::vector<T*> >& bands)
    : analysisType(k2D),
      hist2D(new Hist2DWrapper<T>(hist_name, title, nBinsX, xmin, xmax, nBinsY, ymin, ymax, bands)),
      hist(nullptr) {}

// Constructor from template MnvH2D and vector of universes
template <typename T>
HistHyperDWrapper<T>::HistHyperDWrapper(MnvH2D* h2d,
                                        std::vector<T*>& univs,
                                        bool clear_error_bands) 
    : analysisType(k2D),
      hist2D(new Hist2DWrapper<T>(h2d, univs, clear_error_bands)),
      hist(nullptr) {}

// Constructor from template MnvH2D and map of universes
template <typename T>
HistHyperDWrapper<T>::HistHyperDWrapper(MnvH2D* h2d,
                                        std::map<std::string, std::vector<T*> >& bands,
                                        bool clear_error_bands)
    : analysisType(k2D),
      hist2D(new Hist2DWrapper<T>(h2d, bands, clear_error_bands)),
      hist(nullptr) {}

#ifndef __GCCXML__
// Constructor from a vector of universes and variable bin widths
template <typename T>
HistHyperDWrapper<T>::HistHyperDWrapper(const char* name, const char* title,
                                        const std::vector<double> xBins, 
                                        const std::vector<double> yBins,
                                        std::vector<T*>& univs)
    : HistHyperDWrapper(
          new MnvH2D(name, title, xBins.size() - 1, xBins.data(),
                     yBins.size() - 1, yBins.data()),
          univs) {}

// Constructor from a map of universes and variable bin widths
template <typename T>
HistHyperDWrapper<T>::HistHyperDWrapper(const char* name, const char* title,
                                        const std::vector<double> xBins,
                                        const std::vector<double> yBins,
                                        std::map<std::string, std::vector<T*> >& bands)
    : HistHyperDWrapper(
          new MnvH2D(name, title, xBins.size() - 1, xBins.data(),
                     yBins.size() - 1, yBins.data()),
          bands) {}
#endif //__GCCXML__

// ============================================================================
// Methods
// ============================================================================

// Copy an MnvH1D's CV histo to each of its vertical error band CV histos
template <typename T>
void HistHyperDWrapper<T>::SyncCVHistos() {
    hist->SyncCVHistos();
    hist2D->SyncCVHistos();
}

// TODO: add check to see if 1D or 2D?
template <typename T>
void HistHyperDWrapper<T>::FillUniverse(const T& univ,
                                        const double value, 
                                        const double weight) {
    hist->FillUniverse(univ, value, weight);
}

template <typename T>
void HistHyperDWrapper<T>::FillUniverse(const T* univ,
                                        const double value,
                                        const double weight) {
    hist->FillUniverse(*univ, value, weight);
}

// TODO: add check to see if 1D or 2D? Maybe if not 2D, then fill 1D with x value
template <typename T>
void HistHyperDWrapper<T>::FillUniverse(const T& univ,
                                        const double valueX,
                                        const double valueY,
                                        const double weight) {
    if (analysisType == k1D) {
        hist->FillUniverse(univ, valueX, weight);
    } else {  // if (analysisType == k2D)
        hist2D->FillUniverse(univ, valueX, valueY, weight);
    }
}

template <typename T>
void HistHyperDWrapper<T>::FillUniverse(const T* univ,
                                        const double valueX,
                                        const double valueY,
                                        const double weight) {
    if (analysisType == k1D) {
        hist->FillUniverse(*univ, valueX, weight);
    } else {  // if (analysisType == k2D)
        hist2D->FillUniverse(*univ, valueX, valueY, weight);
    }
}

template <typename T>
TH1D* HistHyperDWrapper<T>::univHist(const T* univ) const {
    return hist->univHist(univ);
}

template <typename T>
TH2D* HistHyperDWrapper<T>::univHist2D(const T* univ) const {
    return hist2D->univHist(univ);
}

template <typename T>
void HistHyperDWrapper<T>::AddUniverses(const std::string name, const T* univ,
                                        const int nhists, const int histno) {
    if (analysisType == k1D) { // 1D already has this method, 2D doesn't
        hist->AddUniverses(name, univ, nhists, histno);
    } else { // if (analysisType == k2D)
        assert(histno < nhists);
        if (name == "cv") {
            if (nhists != 1) {
                std::cerr << "HistHyperDWrapper constructor ERROR: CV should not have "
                             "more than one universe!"
                          << std::endl;
                std::exit(2);
            }

            hist2D->univToHistMap[univ] = hist2D->hist; // TODO: this look right?
            return;
        }
        std::cout << "try to add a " << univ->ShortName() << " as " << name
                  << std::endl;

        // add error band if not already there
        if (!hist2D->hist->HasVertErrorBand(name)) {
            std::cerr << "adding vert error band with " << name << " " << nhists
                      << std::endl;
            hist2D->hist->AddVertErrorBand(name, nhists);
        }
        // Connect each band's universe to the corresponding MnvH1D's TH1
        hist2D->univToHistMap[univ] = hist2D->hist->GetVertErrorBand(name)->GetHist(histno);

        assert(hist2D->hist);
    }
}

#endif  // #ifndef HISTHYPERDWRAPPER_CXX