#ifndef MLVertexSystematics_CXX
#define MLVertexSystematics_CXX

// Created by Anezka Klustova: a.klustova20@imperial.ac.uk/klustova.a@gmail.com
// September 2023
// Neutrino and antineutrino flat normalization uncertainty for ML vertexing

#include "universes/MLVertexSystematics.h"
#include <iostream>

// Helper functions -- get Weighters, containers of systematics universes
namespace PlotUtils {

  template <class T>
  std::vector<T*> GetMLVertexSystematics(typename T::config_t chain) {
    std::vector<T*> ret;
    ret.push_back(new PlotUtils::MLVertexUniverse<T>(chain, -1.));
    ret.push_back(new PlotUtils::MLVertexUniverse<T>(chain, 1. ));

    return ret;
  }

  template <class T>
  std::map< std::string, std::vector<T*> > GetMLVertexSystematicsMap(typename T::config_t chain ) {
    std::map< std::string, std::vector<T*> > ret;
    
    ret["ML_Vertex"].push_back(new PlotUtils::MLVertexUniverse<T>(chain, -1.));
    ret["ML_Vertex"].push_back(new PlotUtils::MLVertexUniverse<T>(chain, 1.));
    std::cout << "ML Vertex Systematics created. " << std::endl;
    return ret;
  }

}

// Class Definitions
// Constructor
template<typename T>
PlotUtils::MLVertexUniverse<T>::MLVertexUniverse(
    typename T::config_t chw, double nsigma) : T(chw, nsigma) {}

template<typename T>
double PlotUtils::MLVertexUniverse<T>::GetMLVertexWeight() const {
  if( T::IsTruth() ) return 1; //No reweighting for truth events 
  
  double cv_wgt = T:: GetMLVertexWeight();

  double wgt_shift = 0;
  // is neutrino? or antineutino
  if( T::GetInt("mc_incoming") == 14 ){
    wgt_shift = NSFDefaults::NuVertexMLUnc; 
  }
  if( T::GetInt("mc_incoming") == -14 ){
    wgt_shift = NSFDefaults::AntinuVertexMLUnc; 
  }

  return  cv_wgt + T::m_nsigma * wgt_shift;
}

template<typename T>
double PlotUtils::MLVertexUniverse<T>::GetWeightRatioToCV() const {

  double cv_wgt = T:: GetMLVertexWeight();

  double wgt_shift = 0;
  // is neutrino? or antineutino
  if( T::GetInt("mc_incoming") == 14 ){
    wgt_shift = NSFDefaults::NuVertexMLUnc; 
  }
  if( T::GetInt("mc_incoming") == -14 ){
    wgt_shift = NSFDefaults::AntinuVertexMLUnc; 
  }

  return  1. + T::m_nsigma*wgt_shift/cv_wgt;
}

template<typename T>
std::string PlotUtils::MLVertexUniverse<T>::ShortName() const { return "ML_Vertex"; }

template<typename T>
std::string PlotUtils::MLVertexUniverse<T>::LatexName() const { return "ML_Vertex"; }

#endif // MLVertexSystematics_CXX
