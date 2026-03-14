#ifndef MLVertexSystematics_H
#define MLVertexSystematics_H

#include "utilities/NSFDefaults.h" // define ML uncertainty in NSd defaults
#include "PlotUtils/TreeWrapper.h"

namespace PlotUtils{
  // To estimate uncertainty on ML vertex for neutrinos and antineutrinos
  // 0.1 % uncertainty assigned based on antineutrino study in MINERvA Document 29111-v1

  template<class T>
  class MLVertexUniverse: public T
  {
    public:
      MLVertexUniverse(typename T::config_t chw, double nsigma );
      
      virtual double GetMLVertexWeight() const /*override*/;
      double GetWeightRatioToCV() const;
      
      virtual std::string ShortName() const /*override*/;
      virtual std::string LatexName() const /*override*/;
      virtual bool IsVerticalOnly()  const  { return true; }/*override*/
  };
  
}

// Explicit Instantiation, needed for loading into python
//! Make sure you put this into Ana/PlotUtils/dict/PlotUtilsDict.h
//template class

// For template classes, the header needs to know about the definitions
#include "MLVertexSystematics.cxx"

#endif // MLVertexSystematics_H
