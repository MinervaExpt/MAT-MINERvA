#ifndef DiffractiveUncSystematics_CXX
#define DiffractiveUncSystematics_CXX

#include "DiffractiveUncSystematics.h"
#include "TargetUtils.h"
#include <iostream>

// Helper functions -- get Weighters, containers of systematics universes
namespace Minerva {

  template <class T>
  std::vector<T*> GetDiffractiveUncSystematics(typename T::config_t chain, 
                                double fracTrkUnc = fracDiffractiveTrackerUnc, double fracTarUnc = fracDiffractiveTargetUnc ) {
    std::vector<T*> ret;
    ret.push_back(new Minerva::DiffractiveUncUniverse<T>(chain, -1., fracTrkUnc, fracTarUnc));
    ret.push_back(new Minerva::DiffractiveUncUniverse<T>(chain, 1. , fracTrkUnc, fracTarUnc));

    return ret;
  }

  template <class T>
  std::map< std::string, std::vector<T*> > GetDiffractiveUncSystematicsMap(typename T::config_t chain, 
                                double fracTrkUnc = fracDiffractiveTrackerUnc, double fracTarUnc = fracDiffractiveTargetUnc ) {
    std::map< std::string, std::vector<T*> > ret;
    
    ret["DiffractiveUnc"].push_back(new Minerva::DiffractiveUncUniverse<T>(chain, -1., fracTrkUnc, fracTarUnc));
    ret["DiffractiveUnc"].push_back(new Minerva::DiffractiveUncUniverse<T>(chain, 1. , fracTrkUnc, fracTarUnc));

    return ret;
  }

}

// Class Definitions
// Constructor
template<typename T>
Minerva::DiffractiveUncUniverse<T>::DiffractiveUncUniverse(
    typename T::config_t chw, double nsigma, double fracTrkUnc, double fracTarUnc ) : 
      T(chw, nsigma), m_fracTrkUnc(fracTrkUnc), m_fracTarUnc(fracTarUnc) {}

template<typename T>
double Minerva::DiffractiveUncUniverse<T>::GetDiffractiveUncWeight() const {
  //Is this a coherent event?
  if( T::GetInt("mc_intType") != 4 ) return 1;

  //Diffractive should really come off of hydrogen only, but this estimating from measurements we have
  double fracDiffUnc = 0;
  //Are you in NTR?
  if( Minerva::TargetUtils::Get().InNukeRegion( T::GetVecElem("mc_vtx",0), 
                            T::GetVecElem("mc_vtx",1), T::GetVecElem("mc_vtx",2) ) ) fracDiffUnc = m_fracTarUnc; 
  //Tracker?
  if( Minerva::TargetUtils::Get().InTracker( T::GetVecElem("mc_vtx",0), 
                            T::GetVecElem("mc_vtx",1), T::GetVecElem("mc_vtx",2) ) ) fracDiffUnc = m_fracTrkUnc; 

  double shift_val = 1 + T::m_nsigma * fracDiffUnc;
  return shift_val;
}

template<typename T>
std::string Minerva::DiffractiveUncUniverse<T>::ShortName() const { return "DiffractiveUnc"; }

template<typename T>
std::string Minerva::DiffractiveUncUniverse<T>::LatexName() const { return "Diffractive Pion Uncertainty"; }

#endif // DiffractiveUncSystematics_CXX
