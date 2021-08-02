#ifndef TARGETMASSSYSTEMATICS_CXX
#define TARGETMASSSYSTEMATICS_CXX

#include "NSFDefaults.h"
#include "TargetMassSystematics.h"
#include "TargetUtils.h"
#include <iostream>

// Helper functions -- get Weighters, containers of systematics universes
namespace PlotUtils {

  template <class T>
  std::vector<T*> GetTargetMassSystematics(typename T::config_t chain ) {
    std::vector<T*> ret;
    
    ret.push_back(new PlotUtils::TargetMassUniverse<T>(chain, -1.));
    ret.push_back(new PlotUtils::TargetMassUniverse<T>(chain, 1.));

    return ret;
  }


  template <class T>
  std::map< std::string, std::vector<T*> > GetTargetMassSystematicsMap(
      typename T::config_t chain ) {
    std::map< std::string, std::vector<T*> > ret;
    
    ret["TargetMass"].push_back(new PlotUtils::TargetMassUniverse<T>(chain, -1.));
    ret["TargetMass"].push_back(new PlotUtils::TargetMassUniverse<T>(chain, 1.));

    return ret;
  }

}


// Class Definitions
// Constructor
template<typename T>
PlotUtils::TargetMassUniverse<T>::TargetMassUniverse(
    typename T::config_t chw, double nsigma) : T(chw, nsigma) {}


template<typename T>
double PlotUtils::TargetMassUniverse<T>::GetTargetMassWeight() const {

  double cv_wgt = T::GetTargetMassWeight();

  int targetZ = T::GetTargetZTrue(); // atomic number of struck nucleus
  int targetVtxZ = T::GetVertexZTrue(); // z-coordinate of interaction vertex

  double wgt_shift;

  // if interaction is not in Nuclear Target region, use CH error
  if(targetVtxZ>PlotUtils::TargetProp::Tracker::Face){
    wgt_shift = NSFDefaults::ch_err;
  }
  // if interaction is in Nuclear Target region, error depends on
  // struck nucleus (one of the nuclear targets, else default to CH)
  else if( PlotUtils::TargetUtils::Get().InWaterTargetMC( T::GetVecElem("mc_vtx",0), T::GetVecElem("mc_vtx",1),
                                                         T::GetVecElem("mc_vtx",2), targetZ ) ){
    wgt_shift = NSFDefaults::h2o_err;
  }
  else if(targetZ==6){
    wgt_shift = NSFDefaults::c_err; 
  }
  else if(targetZ==26){
    wgt_shift = NSFDefaults::fe_err; 
  }
  else if(targetZ==82){
    wgt_shift = NSFDefaults::pb_err; 
  }
  // assume CH if not H2O, C, Fe, or Pb
  else{
    wgt_shift = NSFDefaults::ch_err;
  }
  return  cv_wgt + T::m_nsigma*wgt_shift;
}

//This systematic is special: its weight function is always 1 in the CV
//as of April 2021.  So, it doesn't even have a Reweighter.  Just
//including it in the list of systematics and using PlotUtils::Model applies it.
template <typename T>
double PlotUtils::TargetMassUniverse<T>::GetWeightRatioToCV() const {
  double cv_wgt = 1; //TODO: Watch very carefully for this to change

  int targetZ = T::GetTargetZTrue(); // atomic number of struck nucleus
  int targetVtxZ = T::GetVertexZTrue(); // z-coordinate of interaction vertex

  double wgt_shift;

  // if interaction is not in Nuclear Target regio, use CH error
  if(targetVtxZ>PlotUtils::TargetProp::Tracker::Face){
    wgt_shift = NSFDefaults::ch_err;
  }
  // if interaction is in Nuclear Target region, error depends on
  // struck nucleus (one of the nuclear targets, else default to CH)
  else if(targetZ==6){
    wgt_shift = NSFDefaults::c_err;
  }
  else if(targetZ==26){
    wgt_shift = NSFDefaults::fe_err;
  }
  else if(targetZ==82){
    wgt_shift = NSFDefaults::pb_err;
  }
  // assume CH if not C, Fe, or Pb
  else{
    wgt_shift = NSFDefaults::ch_err;
  }
  return  1 + T::m_nsigma*wgt_shift/cv_wgt;
}

// This method returns the MnvHnD that users will include in the denominator of
// their cross section calculation to normalize by the number of targets.
// It returns an MnvHnD, in the binning specified by the user, where the content of
// each bin is the number of targets (also specified by the user) in every systematic
// universe except the "TargetMass" universes which are appropriately shifted
template<class MnvHistoType>
MnvHistoType* PlotUtils::GetNTargetsHist(double nTargets,
                                         int targetZ,
                                         MnvHistoType* template_hist)
{

  // This is the hist that will be returned to the user. Start from a clean slate
  MnvHistoType* h_number_of_targets = (MnvHistoType*)template_hist->Clone("number_of_targets");
  h_number_of_targets->ClearAllErrorBands();
  h_number_of_targets->Reset();

  // Fetch correct number of targets
  double nTargets = 2.0;

  //CV first
  for(int i=0;i<h_flux_number_of_targets->GetSize();i++)
    h_flux_number_of_targets->SetBinContent(i,nTargets);

  // Systematically vary number of targets in "TargetMass" universes
  MnvVertErrorBand *errBand = h_flux_ppfx->GetVertErrorBand("TargetMass");
  const int universes = errBand->GetNHists();
  auto flux_sys_hists = GetVector(h_flux_integrated);
  for(int u=0;u<universes;++u) {
    TH1D* tmp_flux = new TH1D(*errBand->GetHist( u ));
    auto tmp_template = h_flux_integrated->GetCVHistoWithStatError();
    tmp_template.SetName(Form("Flux_integrated_universe_%d",u));
    double flux_uni = tmp_flux->Integral(ppfx_b_min,ppfx_b_max,"width");
    for(int i=0;i<h_flux_integrated->GetSize();i++)
      tmp_template.SetBinContent(i,flux_uni);
    flux_sys_hists.push_back(NewHist(tmp_template));
  }

  // Push the constructed error band into the return hist
  h_number_of_targets->AddVertErrorBand("TargetMass_{0}".format(targetZ),flux_sys_hists);
  
  // Fill CV for any systematic universes which don't have a non-CV nTargets
  h_number_of_targets->AddMissingErrorBandsAndFillWithCV(*template_hist);
  return h_number_of_targets;

}

template<typename T>
std::string PlotUtils::TargetMassUniverse<T>::ShortName() const { return "Target_Mass"; }


template<typename T>
std::string PlotUtils::TargetMassUniverse<T>::LatexName() const { return "Target Mass"; }


#endif // TARGETMASSSYSTEMATICS_CXX
