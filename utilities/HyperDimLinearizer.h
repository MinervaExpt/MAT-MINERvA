#ifndef MNV_HYPERDIMLINEARIZER_h
#define MNV_HYPERDIMLINEARIZER_h 1
#include <iostream>
#include <utility>
#include <vector>

#include "PlotUtils/MnvH2D.h"
namespace PlotUtils {
enum EAnalysisType { k2D,
                     k1D,
                     k2D_lite,
                     k1D_lite };  
class HyperDimLinearizer {
   public:
    HyperDimLinearizer(std::vector<std::vector<double> > input, int type);                       // constructor
    HyperDimLinearizer(std::vector<std::vector<double> > input, EAnalysisType type);             // constructor

    std::pair<int, int> GetBin(std::vector<double> values);                                      // Template get bin for 2,3,4D cases
    std::vector<int> GetValues(int x);                                                           //
    std::vector<TH2D*> Get2DHistos(PlotUtils::MnvH2D* result, bool IncludeSys);                  // This is for type==0
    std::vector<PlotUtils::MnvH2D*> Get2DMnvHistos(PlotUtils::MnvH2D* result, bool IncludeSys);  // This is for type==0
    TH2D* Get2DHisto(PlotUtils::MnvH1D* result, bool IncludeSys);                                // This is for type==1, 2D result only!!
    PlotUtils::MnvH2D* Get2DMnvHisto(PlotUtils::MnvH1D* result, bool IncludeSys);                // This is for type==1, 2D result only!!
    void TestFunctionality();                                                                    // a bunch of prints.
    bool IsUnderflow(int lin_bin, int axis);                                                     // Check if a bin in linearized is an underflow bin in phase space, axis defaulted
    bool IsOverflow(int lin_bin, int axis);                                                      // Check if a bin in linearized is an overflow bin in phase space, axis defaulted
   private:
    int Get1DBin(double value, int el);  // Get the bin number for one of the dimensions

    // Type 0 is 2D (leaving y axis alone), Type 1 is 1D; these bin under/overflow through out
    // Type 2 is like Type 0 (2D), Type 3 is like Type 1 (1D); these use a global under/overflow bins
    // int m_analysis_type;
    EAnalysisType m_analysis_type;
    std::vector<int> m_el_size;                 // how many bins are we talking about
    std::vector<int> m_cell_size;               // For each axis, how big are the cells of a single bin of that axis in bin space
    std::vector<std::vector<double> > m_invec;  // internal vector of boundaries
    int m_n_global_x_bins;
};
}  // namespace PlotUtils
#endif
