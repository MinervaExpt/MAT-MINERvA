#ifndef TARGETMASSSYSTEMATICS_CXX
#define TARGETMASSSYSTEMATICS_CXX

#include "utilities/NSFDefaults.h"
#include "universes/TargetMassSystematics.h"
#include "utilities/TargetUtils.h"
#include "MnvH1D.h"
#include "MnvH2D.h"
#include <iostream>

// Helper functions -- get Weighters, containers of systematics universes
namespace PlotUtils {

  template <class T>
  std::vector<T*> GetTargetMassSystematics(typename T::config_t chain ) {
    std::vector<T*> ret;
    
    ret.push_back(new PlotUtils::TargetMassScintillatorUniverse<T>(chain, -1.));
    ret.push_back(new PlotUtils::TargetMassScintillatorUniverse<T>(chain, 1.));
    ret.push_back(new PlotUtils::TargetMassCarbonUniverse<T>(chain, -1.));
    ret.push_back(new PlotUtils::TargetMassCarbonUniverse<T>(chain, 1.));
    ret.push_back(new PlotUtils::TargetMassWaterUniverse<T>(chain, -1.));
    ret.push_back(new PlotUtils::TargetMassWaterUniverse<T>(chain, 1.));
    ret.push_back(new PlotUtils::TargetMassIronUniverse<T>(chain, -1.));
    ret.push_back(new PlotUtils::TargetMassIronUniverse<T>(chain, 1.));
    ret.push_back(new PlotUtils::TargetMassLeadUniverse<T>(chain, -1.));
    ret.push_back(new PlotUtils::TargetMassLeadUniverse<T>(chain, 1.));

    return ret;
  }


  template <class T>
  std::map< std::string, std::vector<T*> > GetTargetMassSystematicsMap(
      typename T::config_t chain ) {
    std::map< std::string, std::vector<T*> > ret;
    
    ret["Target_Mass_CH"].push_back(new PlotUtils::TargetMassScintillatorUniverse<T>(chain, -1.));
    ret["Target_Mass_CH"].push_back(new PlotUtils::TargetMassScintillatorUniverse<T>(chain, 1.));
    ret["Target_Mass_C"].push_back(new PlotUtils::TargetMassCarbonUniverse<T>(chain, -1.));
    ret["Target_Mass_C"].push_back(new PlotUtils::TargetMassCarbonUniverse<T>(chain, 1.));
    ret["Target_Mass_H2O"].push_back(new PlotUtils::TargetMassWaterUniverse<T>(chain, -1.));
    ret["Target_Mass_H2O"].push_back(new PlotUtils::TargetMassWaterUniverse<T>(chain, 1.));
    ret["Target_Mass_Fe"].push_back(new PlotUtils::TargetMassIronUniverse<T>(chain, -1.));
    ret["Target_Mass_Fe"].push_back(new PlotUtils::TargetMassIronUniverse<T>(chain, 1.));
    ret["Target_Mass_Pb"].push_back(new PlotUtils::TargetMassLeadUniverse<T>(chain, -1.));
    ret["Target_Mass_Pb"].push_back(new PlotUtils::TargetMassLeadUniverse<T>(chain, 1.));

    return ret;
  }

  // Each of these methods returns the MnvHnD that users will include in the denominator of
  // their cross section calculation to normalize by the number of targets (based on which target)
  // Each returns an MnvHnD, in the binning specified by the user, where the content of
  // each bin is the number of targets (also specified by the user) in every systematic
  // universe except the "Target_Mass_X" universes which are appropriately shifted

  // Scintillator
  template<class MnvHistoType>
  MnvHistoType* GetNTargetsScintillatorHist(double nTargets, MnvHistoType* template_hist){
  
    // This is the hist that will be returned to the user. Start from a clean slate
    MnvHistoType* h_number_of_targets = (MnvHistoType*)template_hist->Clone("number_of_targets");
    h_number_of_targets->ClearAllErrorBands();
    h_number_of_targets->Reset();
  
    // Populate the CV hist of the MnvH1D first
    // Note that the CV number of targets is provided by the user
    for(int i=0;i<h_number_of_targets->GetSize();i++){
      h_number_of_targets->SetBinContent(i,nTargets);
    }
  
    // Create the only error band affected by this variation
    h_number_of_targets->AddVertErrorBand("Target_Mass_CH",2);
    for(int i=0;i<h_number_of_targets->GetSize();i++){
      h_number_of_targets->GetVertErrorBand("Target_Mass_CH")->GetHist(0)->SetBinContent(i,nTargets*(1-NSFDefaults::ch_err));
      h_number_of_targets->GetVertErrorBand("Target_Mass_CH")->GetHist(1)->SetBinContent(i,nTargets*(1+NSFDefaults::ch_err));
    }
  
    // Fill CV for any systematic universes which don't have a non-CV nTargets
    h_number_of_targets->AddMissingErrorBandsAndFillWithCV(*template_hist);
    return h_number_of_targets;
  
  }

  template MnvH1D* GetNTargetsScintillatorHist<MnvH1D>(double nTargets,
                                                      MnvH1D* template_hist);
  
  template MnvH2D* GetNTargetsScintillatorHist<MnvH2D>(double nTargets,
                                                      MnvH2D* template_hist);
 
  // Carbon 
  template<class MnvHistoType>
  MnvHistoType* GetNTargetsCarbonHist(double nTargets, MnvHistoType* template_hist){
  
    // This is the hist that will be returned to the user. Start from a clean slate
    MnvHistoType* h_number_of_targets = (MnvHistoType*)template_hist->Clone("number_of_targets");
    h_number_of_targets->ClearAllErrorBands();
    h_number_of_targets->Reset();
  
    // Populate the CV hist of the MnvH1D first
    // Note that the CV number of targets is provided by the user
    for(int i=0;i<h_number_of_targets->GetSize();i++){
      h_number_of_targets->SetBinContent(i,nTargets);
    }
  
    // Create the only error band affected by this variation
    h_number_of_targets->AddVertErrorBand("Target_Mass_C",2);
    for(int i=0;i<h_number_of_targets->GetSize();i++){
      h_number_of_targets->GetVertErrorBand("Target_Mass_C")->GetHist(0)->SetBinContent(i,nTargets*(1-NSFDefaults::c_err));
      h_number_of_targets->GetVertErrorBand("Target_Mass_C")->GetHist(1)->SetBinContent(i,nTargets*(1+NSFDefaults::c_err));
    }
  
    // Fill CV for any systematic universes which don't have a non-CV nTargets
    h_number_of_targets->AddMissingErrorBandsAndFillWithCV(*template_hist);
    return h_number_of_targets;
  
  }

  template MnvH1D* GetNTargetsCarbonHist<MnvH1D>(double nTargets,
                                                 MnvH1D* template_hist);
  
  template MnvH2D* GetNTargetsCarbonHist<MnvH2D>(double nTargets,
                                                 MnvH2D* template_hist);
 
  // Water 
  template<class MnvHistoType>
  MnvHistoType* GetNTargetsWaterHist(double nTargets, MnvHistoType* template_hist){
  
    // This is the hist that will be returned to the user. Start from a clean slate
    MnvHistoType* h_number_of_targets = (MnvHistoType*)template_hist->Clone("number_of_targets");
    h_number_of_targets->ClearAllErrorBands();
    h_number_of_targets->Reset();
  
    // Populate the CV hist of the MnvH1D first
    // Note that the CV number of targets is provided by the user
    for(int i=0;i<h_number_of_targets->GetSize();i++){
      h_number_of_targets->SetBinContent(i,nTargets);
    }
  
    // Create the only error band affected by this variation
    h_number_of_targets->AddVertErrorBand("Target_Mass_H2O",2);
    for(int i=0;i<h_number_of_targets->GetSize();i++){
      h_number_of_targets->GetVertErrorBand("Target_Mass_H2O")->GetHist(0)->SetBinContent(i,nTargets*(1-NSFDefaults::h2o_err));
      h_number_of_targets->GetVertErrorBand("Target_Mass_H2O")->GetHist(1)->SetBinContent(i,nTargets*(1+NSFDefaults::h2o_err));
    }
  
    // Fill CV for any systematic universes which don't have a non-CV nTargets
    h_number_of_targets->AddMissingErrorBandsAndFillWithCV(*template_hist);
    return h_number_of_targets;
  
  }

  template MnvH1D* GetNTargetsWaterHist<MnvH1D>(double nTargets,
                                                MnvH1D* template_hist);
  
  template MnvH2D* GetNTargetsWaterHist<MnvH2D>(double nTargets,
                                                MnvH2D* template_hist);
 
  // Iron 
  template<class MnvHistoType>
  MnvHistoType* GetNTargetsIronHist(double nTargets, MnvHistoType* template_hist){
  
    // This is the hist that will be returned to the user. Start from a clean slate
    MnvHistoType* h_number_of_targets = (MnvHistoType*)template_hist->Clone("number_of_targets");
    h_number_of_targets->ClearAllErrorBands();
    h_number_of_targets->Reset();
  
    // Populate the CV hist of the MnvH1D first
    // Note that the CV number of targets is provided by the user
    for(int i=0;i<h_number_of_targets->GetSize();i++){
      h_number_of_targets->SetBinContent(i,nTargets);
    }
  
    // Create the only error band affected by this variation
    h_number_of_targets->AddVertErrorBand("Target_Mass_Fe",2);
    for(int i=0;i<h_number_of_targets->GetSize();i++){
      h_number_of_targets->GetVertErrorBand("Target_Mass_Fe")->GetHist(0)->SetBinContent(i,nTargets*(1-NSFDefaults::fe_err));
      h_number_of_targets->GetVertErrorBand("Target_Mass_Fe")->GetHist(1)->SetBinContent(i,nTargets*(1+NSFDefaults::fe_err));
    }
  
    // Fill CV for any systematic universes which don't have a non-CV nTargets
    h_number_of_targets->AddMissingErrorBandsAndFillWithCV(*template_hist);
    return h_number_of_targets;
  
  }

  template MnvH1D* GetNTargetsIronHist<MnvH1D>(double nTargets,
                                                      MnvH1D* template_hist);
  
  template MnvH2D* GetNTargetsIronHist<MnvH2D>(double nTargets,
                                                      MnvH2D* template_hist);
 
  // Lead 
  template<class MnvHistoType>
  MnvHistoType* GetNTargetsLeadHist(double nTargets, MnvHistoType* template_hist){
  
    // This is the hist that will be returned to the user. Start from a clean slate
    MnvHistoType* h_number_of_targets = (MnvHistoType*)template_hist->Clone("number_of_targets");
    h_number_of_targets->ClearAllErrorBands();
    h_number_of_targets->Reset();
  
    // Populate the CV hist of the MnvH1D first
    // Note that the CV number of targets is provided by the user
    for(int i=0;i<h_number_of_targets->GetSize();i++){
      h_number_of_targets->SetBinContent(i,nTargets);
    }
  
    // Create the only error band affected by this variation
    h_number_of_targets->AddVertErrorBand("Target_Mass_Pb",2);
    for(int i=0;i<h_number_of_targets->GetSize();i++){
      h_number_of_targets->GetVertErrorBand("Target_Mass_Pb")->GetHist(0)->SetBinContent(i,nTargets*(1-NSFDefaults::pb_err));
      h_number_of_targets->GetVertErrorBand("Target_Mass_Pb")->GetHist(1)->SetBinContent(i,nTargets*(1+NSFDefaults::pb_err));
    }
  
    // Fill CV for any systematic universes which don't have a non-CV nTargets
    h_number_of_targets->AddMissingErrorBandsAndFillWithCV(*template_hist);
    return h_number_of_targets;
  
  }

  template MnvH1D* GetNTargetsLeadHist<MnvH1D>(double nTargets,
                                                      MnvH1D* template_hist);
  
  template MnvH2D* GetNTargetsLeadHist<MnvH2D>(double nTargets,
                                                      MnvH2D* template_hist);
  
}


// Class Definitions Scintillator
// Constructor
template<typename T>
PlotUtils::TargetMassScintillatorUniverse<T>::TargetMassScintillatorUniverse(
    typename T::config_t chw, double nsigma) : T(chw, nsigma) {}

template<typename T>
double PlotUtils::TargetMassScintillatorUniverse<T>::GetTargetMassWeight() const {

  double cv_wgt = T::GetTargetMassWeight();

  int targetVtxZ = T::GetVertexZTrue(); // z-coordinate of interaction vertex

  double wgt_shift;

  // If interaction is not in Nuclear Target region, use CH error
  if(targetVtxZ>PlotUtils::TargetProp::Tracker::Face){
    wgt_shift = NSFDefaults::ch_err;
  }
  // If interaction is in Nuclear Target region, error depends on
  // struck nucleus (one of the nuclear targets, else default to CH)
  else if( PlotUtils::TargetUtils::Get().InPassiveTargetVolMC( T::GetVecElem("mc_vtx",0), T::GetVecElem("mc_vtx",1),
                                                               T::GetVecElem("mc_vtx",2)) ){
    wgt_shift = 0.;
  } 
  // assume CH if not H2O, C, Fe, or Pb
  else{
    wgt_shift = NSFDefaults::ch_err;
  }
  return  cv_wgt + T::m_nsigma*wgt_shift;
}

template<typename T>
double PlotUtils::TargetMassScintillatorUniverse<T>::GetWeightRatioToCV() const {

  double cv_wgt = T::GetTargetMassWeight();

  int targetVtxZ = T::GetVertexZTrue(); // z-coordinate of interaction vertex

  double wgt_shift;

  // If interaction is not in Nuclear Target region, use CH error
  if(targetVtxZ>PlotUtils::TargetProp::Tracker::Face){
    wgt_shift = NSFDefaults::ch_err;
  }
  // If interaction is in Nuclear Target region, error depends on
  // struck nucleus (one of the nuclear targets, else default to CH)
  else if( PlotUtils::TargetUtils::Get().InPassiveTargetVolMC( T::GetVecElem("mc_vtx",0), T::GetVecElem("mc_vtx",1),
                                                               T::GetVecElem("mc_vtx",2)) ){
    wgt_shift = 0.;
  } 
  // assume CH if not H2O, C, Fe, or Pb
  else{
    wgt_shift = NSFDefaults::ch_err;
  }
  return  1. + T::m_nsigma*wgt_shift/cv_wgt;
}

template<typename T>
std::string PlotUtils::TargetMassScintillatorUniverse<T>::ShortName() const { return "Target_Mass_CH"; }

template<typename T>
std::string PlotUtils::TargetMassScintillatorUniverse<T>::LatexName() const { return "Target Mass CH"; }

// Class Definitions Carbon
// Constructor
template<typename T>
PlotUtils::TargetMassCarbonUniverse<T>::TargetMassCarbonUniverse(
    typename T::config_t chw, double nsigma) : T(chw, nsigma) {}

template<typename T>
double PlotUtils::TargetMassCarbonUniverse<T>::GetTargetMassWeight() const {

  double cv_wgt = T::GetTargetMassWeight();

  double wgt_shift;

  if( PlotUtils::TargetUtils::Get().InCarbonTargetVolMC( T::GetVecElem("mc_vtx",0), T::GetVecElem("mc_vtx",1),
                                                         T::GetVecElem("mc_vtx",2)) ){
    wgt_shift = NSFDefaults::c_err; 
  }
  // Only shift interactions that truly occurred in carbon (targets)
  else{
    wgt_shift = 0.;
  }
  return  cv_wgt + T::m_nsigma*wgt_shift;
}

template<typename T>
double PlotUtils::TargetMassCarbonUniverse<T>::GetWeightRatioToCV() const {

  double cv_wgt = T::GetTargetMassWeight();

  double wgt_shift;

  if( PlotUtils::TargetUtils::Get().InCarbonTargetVolMC( T::GetVecElem("mc_vtx",0), T::GetVecElem("mc_vtx",1),
                                                       T::GetVecElem("mc_vtx",2)) ){
    wgt_shift = NSFDefaults::c_err; 
  }
  // Only shift interactions that truly occurred in lead
  else{
    wgt_shift = 0.;
  }
  return  1. + T::m_nsigma*wgt_shift/cv_wgt;
}

template<typename T>
std::string PlotUtils::TargetMassCarbonUniverse<T>::ShortName() const { return "Target_Mass_C"; }

template<typename T>
std::string PlotUtils::TargetMassCarbonUniverse<T>::LatexName() const { return "Target Mass C"; }

// Class Definitions Water
// Constructor
template<typename T>
PlotUtils::TargetMassWaterUniverse<T>::TargetMassWaterUniverse(
    typename T::config_t chw, double nsigma) : T(chw, nsigma) {}

template<typename T>
double PlotUtils::TargetMassWaterUniverse<T>::GetTargetMassWeight() const {

  double cv_wgt = T::GetTargetMassWeight();

  int targetZ = T::GetTargetZTrue(); // atomic number of struck nucleus

  double wgt_shift;

  if( PlotUtils::TargetUtils::Get().InWaterTargetMC( T::GetVecElem("mc_vtx",0), T::GetVecElem("mc_vtx",1),
                                                     T::GetVecElem("mc_vtx",2), targetZ ) ){
    wgt_shift = NSFDefaults::h2o_err;
  }
  // Only shift interactions that truly occurred in water
  else{
    wgt_shift = 0.;
  }
  return  cv_wgt + T::m_nsigma*wgt_shift;
}

template<typename T>
double PlotUtils::TargetMassWaterUniverse<T>::GetWeightRatioToCV() const {

  double cv_wgt = T::GetTargetMassWeight();

  int targetZ = T::GetTargetZTrue(); // atomic number of struck nucleus

  double wgt_shift;

  if( PlotUtils::TargetUtils::Get().InWaterTargetMC( T::GetVecElem("mc_vtx",0), T::GetVecElem("mc_vtx",1),
                                                     T::GetVecElem("mc_vtx",2), targetZ ) ){
    wgt_shift = NSFDefaults::h2o_err;
  }
  // Only shift interactions that truly occurred in water
  else{
    wgt_shift = 0.;
  }
  return  1. + T::m_nsigma*wgt_shift/cv_wgt;
}

template<typename T>
std::string PlotUtils::TargetMassWaterUniverse<T>::ShortName() const { return "Target_Mass_H2O"; }

template<typename T>
std::string PlotUtils::TargetMassWaterUniverse<T>::LatexName() const { return "Target Mass H2O"; }

// Class Definitions Iron 
// Constructor
template<typename T>
PlotUtils::TargetMassIronUniverse<T>::TargetMassIronUniverse(
    typename T::config_t chw, double nsigma) : T(chw, nsigma) {}

template<typename T>
double PlotUtils::TargetMassIronUniverse<T>::GetTargetMassWeight() const {

  double cv_wgt = T::GetTargetMassWeight();

  double wgt_shift;

  if( PlotUtils::TargetUtils::Get().InIronTargetVolMC( T::GetVecElem("mc_vtx",0), T::GetVecElem("mc_vtx",1),
                                                       T::GetVecElem("mc_vtx",2)) ){
    wgt_shift = NSFDefaults::fe_err; 
  }
  // Only shift interactions that truly occurred in iron 
  else{
    wgt_shift = 0.;
  }
  return  cv_wgt + T::m_nsigma*wgt_shift;
}

template<typename T>
double PlotUtils::TargetMassIronUniverse<T>::GetWeightRatioToCV() const {

  double cv_wgt = T::GetTargetMassWeight();

  double wgt_shift;

  if( PlotUtils::TargetUtils::Get().InIronTargetVolMC( T::GetVecElem("mc_vtx",0), T::GetVecElem("mc_vtx",1),
                                                       T::GetVecElem("mc_vtx",2)) ){
    wgt_shift = NSFDefaults::fe_err; 
  }
  // Only shift interactions that truly occurred in lead
  else{
    wgt_shift = 0.;
  }
  return  1. + T::m_nsigma*wgt_shift/cv_wgt;
}

template<typename T>
std::string PlotUtils::TargetMassIronUniverse<T>::ShortName() const { return "Target_Mass_Fe"; }

template<typename T>
std::string PlotUtils::TargetMassIronUniverse<T>::LatexName() const { return "Target Mass Fe"; }

// Class Definitions Lead
// Constructor
template<typename T>
PlotUtils::TargetMassLeadUniverse<T>::TargetMassLeadUniverse(
    typename T::config_t chw, double nsigma) : T(chw, nsigma) {}

template<typename T>
double PlotUtils::TargetMassLeadUniverse<T>::GetTargetMassWeight() const {

  double cv_wgt = T::GetTargetMassWeight();

  double wgt_shift;

  if( PlotUtils::TargetUtils::Get().InLeadTargetVolMC( T::GetVecElem("mc_vtx",0), T::GetVecElem("mc_vtx",1),
                                                       T::GetVecElem("mc_vtx",2)) ){
    wgt_shift = NSFDefaults::pb_err; 
  }
  // Only shift interactions that truly occurred in lead
  else{
    wgt_shift = 0.;
  }
  return  cv_wgt + T::m_nsigma*wgt_shift;
}

template<typename T>
double PlotUtils::TargetMassLeadUniverse<T>::GetWeightRatioToCV() const {

  double cv_wgt = T::GetTargetMassWeight();

  double wgt_shift;

  if( PlotUtils::TargetUtils::Get().InLeadTargetVolMC( T::GetVecElem("mc_vtx",0), T::GetVecElem("mc_vtx",1),
                                                       T::GetVecElem("mc_vtx",2)) ){
    wgt_shift = NSFDefaults::pb_err; 
  }
  // Only shift interactions that truly occurred in lead
  else{
    wgt_shift = 0.;
  }
  return  1. + T::m_nsigma*wgt_shift/cv_wgt;
}

template<typename T>
std::string PlotUtils::TargetMassLeadUniverse<T>::ShortName() const { return "Target_Mass_Pb"; }

template<typename T>
std::string PlotUtils::TargetMassLeadUniverse<T>::LatexName() const { return "Target Mass Pb"; }


#endif // TARGETMASSSYSTEMATICS_CXX
