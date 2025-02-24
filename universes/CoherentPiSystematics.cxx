#ifndef CoherentPiSystematics_CXX
#define CoherentPiSystematics_CXX

#include "universes/CoherentPiSystematics.h"
#include "utilities/TargetUtils.h"
#include <iostream>

// Helper functions -- get Weighters, containers of systematics universes
namespace PlotUtils {

  template <class T>
  std::vector<T*> GetCoherentPiSystematics(typename T::config_t chain,
                              double fracTrkUnc = fracCoherentPiUncTracker, double fracTarUnc = fracCoherentPiUncTarget ) {
    std::vector<T*> ret;
    ret.push_back(new PlotUtils::CoherentPiUniverse<T>( chain, -1., fracTrkUnc, fracTarUnc ));
    ret.push_back(new PlotUtils::CoherentPiUniverse<T>( chain, 1. , fracTrkUnc, fracTarUnc ));

    return ret;
  }

  template <class T>
  std::map< std::string, std::vector<T*> > GetCoherentPiSystematicsMap(typename T::config_t chain,
                              double fracTrkUnc = fracCoherentPiUncTracker, double fracTarUnc = fracCoherentPiUncTarget ) {
    std::map< std::string, std::vector<T*> > ret;
    
    ret["CoherentPiUnc"].push_back(new PlotUtils::CoherentPiUniverse<T>( chain, -1., fracTrkUnc, fracTarUnc ));
    ret["CoherentPiUnc"].push_back(new PlotUtils::CoherentPiUniverse<T>( chain, 1. , fracTrkUnc, fracTarUnc ));

    return ret;
  }

}

// Class Definitions
// Constructor

template<typename T>
PlotUtils::CoherentPiUniverse<T>::CoherentPiUniverse(
    typename T::config_t chw, double nsigma, double fracTrkUnc, double fracTarUnc ) : 
      T(chw, nsigma), m_fracTrkUnc(fracTrkUnc), m_fracTarUnc(fracTarUnc) {}

template<typename T>
double PlotUtils::CoherentPiUniverse<T>::GetCoherentPiWeight(double thpi_true, double tpi_true) const {
  //Is this a coherent event?
  if( T::GetInt("mc_intType") != 4 ) return 1;

    //Coherent uncertainty different in the nuke region vs tracker region
    double fracCohPiUnc = 0;
    //Are you in NTR?
    if( PlotUtils::TargetUtils::Get().InNukeRegion( T::GetVecElem("mc_vtx",0), 
                              T::GetVecElem("mc_vtx",1), T::GetVecElem("mc_vtx",2) ) ) fracCohPiUnc = m_fracTarUnc; 
    //Tracker?
    if( PlotUtils::TargetUtils::Get().InTracker( T::GetVecElem("mc_vtx",0), 
                              T::GetVecElem("mc_vtx",1), T::GetVecElem("mc_vtx",2) ) ) fracCohPiUnc = m_fracTrkUnc; 


  double shift_val = 1 + T::m_nsigma * fracCohPiUnc;
  return shift_val; // do I need to multiply by the weight or not, we shall see
}

template<typename T>
std::string PlotUtils::CoherentPiUniverse<T>::ShortName() const { return "CoherentPiUnc"; }

template<typename T>
std::string PlotUtils::CoherentPiUniverse<T>::LatexName() const { return "Coherent Pion Uncertainty"; }

#endif // CoherentSystematics_CXX
