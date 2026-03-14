#ifndef COHERENTPISYSTEMATICS_H
#define COHERENTPISYSTEMATICS_H

#include "PlotUtils/TreeWrapper.h"
#include "utilities/NSFDefaults.h"

namespace PlotUtils{
  //Right now, just hardcoding this.  
  //To estimate uncertainty on diffractive pion events, putting an uncertainty on coherent events
  //  to simulate our lack of knowledge about diffractive events 
  //  50% for tracker, 20% for nuclear targets
  
  //Looking at Alex's coherent xsec in plastic for pion T<500 MeV
  // The data w/ error is close to MC 
  // The precision on the data is about 20%
  // (yes, plastic, but we're assuming its roughly the same for targets)
  // 24/02/2025: Anezka lower unc looking at data err in Alex's paper, about 10% in scintillator, nuclear target about 30%
  static double fracCoherentPiUncTracker = 0.1; // orig 0.2
  static double fracCoherentPiUncTarget = 0.3; // orig 0.4

  template<class T>
  class CoherentPiUniverse: public T
  {
    public:
      CoherentPiUniverse(typename T::config_t chw, double nsigma, double fracTrkUnc = fracCoherentPiUncTracker , double fracTarUnc = fracCoherentPiUncTarget );
      
      virtual double GetCoherentPiWeight( double thpi_true /*deg*/, double tpi_true /*GeV*/ ) const /*override*/;
      
      virtual std::string ShortName() const /*override*/;
      virtual std::string LatexName() const /*override*/;
      virtual bool IsVerticalOnly()  const  { return true; }/*override*/
    private:
      double m_fracTrkUnc;
      double m_fracTarUnc;
  };
  
  }

// Explicit Instantiation, needed for loading into python
//! Make sure you put this into Ana/PlotUtils/dict/PlotUtilsDict.h
//template class

// For template classes, the header needs to know about the definitions
#include "CoherentPiSystematics.cxx"

#endif // COHERENTPISYSTEMATICS_H

