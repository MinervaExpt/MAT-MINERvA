#ifndef PLOTUTILSDICT_H 
#define PLOTUTILSDICT_H 1

// Include files for PlotUtils dictionary.

/** @file PlotUtilsDict.h
 *  
 *
 *  @author Jeremy Wolcott <jwolcott@fnal.gov>
 *  @date   2012-11-25
 */
// ============================================================================
// PlotUtils
// ============================================================================

// here we need to include all the header files
// for the classes we want to make dictionaries of

#include <vector>

/*#include "NSFDefaults.h"
#include "AnaBinning.h"
#include "ArachneUtils.h"
#include "Exceptions.h"
#include "HistogramUtils.h"
#include "MnvFluxConstraint.h" 
#include "MnvAnaTuple.h"
#include "MnvPlotter.h"
#include "MnvRecoShifter.h"
#include "TargetUtils.h"
#include "MnvNormalization.h"
#include "POTCounter.h"
#include "FluxReweighter.h"
#include "FluxReweighterWithWiggleFit.h"
#include "HyperDimLinearizer.h"
#include "PhysicsVariables.h"
#include "MnvColors.h"
#include "GridCanvas.h"*/

// PlotUtils weight classes
#include "weighters/weightRPA.h"
#include "weighters/weight_2p2h.h"
#include "weighters/weightLowQ2Pi.h"
#include "weighters/weightDIS.h"
#include "weighters/weightZExp.h"

//PlotUtils systematic universes classes (new sys framework)
#include "universes/MinervaUniverse.h"
#include "universes/GenieSystematics.h"
#include "universes/GeantHadronSystematics.h"
#include "universes/MnvTuneSystematics.h"
#include "universes/MinosEfficiencySystematics.h"
#include "universes/MuonSystematics.h"
#include "universes/MuonResolutionSystematics.h"
#include "universes/MichelSystematics.h"
#include "universes/AngleSystematics.h"
#include "universes/ResponseSystematics.h"


//TODO: Do I need this?
//#include "ErrorHandler.h"

// this garbage is necessary so that gccxml is able to create dictionaries for these custom containers
// (since it otherwise doesn't know which specific version of these templated classes to instantiate)
// see: http://root.cern.ch/root/roottalk/roottalk10/0035.html
// somehow std::map<>s seem to be instantiated somewhere else, so explicit instantiation is not necessary?
#ifdef __GCCXML__
/*template class std::vector<PlotUtils::MnvEVD::Event>;                                       // the 'Events' typedef
template class std::pair<std::string, std::vector<PlotUtils::MnvEVD::Event> >;              // the 'EventGroup' typedef

template class std::map<std::string, std::vector<std::string> >;
template class std::pair<std::string, std::vector<std::string> >;*/

// The std::pair<>s for those std::map<>s don't seem to be generated though.
/*template class std::pair< std::string, PlotUtils::MnvLatErrorBand* >;
template class std::pair< std::string, PlotUtils::MnvVertErrorBand* >;
template class std::pair< std::string, TMatrixT<double>* >;
template class std::pair< std::string, PlotUtils::MnvVertErrorBand2D* >;
template class std::pair< std::string, PlotUtils::MnvLatErrorBand2D* >;
template class std::pair< std::string, PlotUtils::MnvVertErrorBand3D* >;
template class std::pair< std::string, PlotUtils::MnvLatErrorBand3D* >;*/

// Use extern keyword because these functions are instantiated in FluxReweighter.cxx already
/*extern template void PlotUtils::FluxReweighter::AddFluxErrorBand<PlotUtils::MnvH1D>(PlotUtils::MnvH1D*);
extern template void PlotUtils::FluxReweighter::AddFluxErrorBand<PlotUtils::MnvH2D>(PlotUtils::MnvH2D*);
extern  template MnvH1D* FluxReweighter::GetIntegratedFluxReweighted<MnvH1D>( int nuPDG,
                                                                 MnvH1D* template_hist,
                                                                 double min_energy,
                                                                 double max_energy,
                                                                 bool use_muon_correlations);
extern template MnvH2D* FluxReweighter::GetIntegratedFluxReweighted<MnvH2D>( int nuPDG,
                                                                 MnvH2D* template_hist,
                                                                 double min_energy,
                                                                 double max_energy,
                                                                 bool use_muon_correlations);*/

#endif
#endif // PLOTUTILSDICT_H

