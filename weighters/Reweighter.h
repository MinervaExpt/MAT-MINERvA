//File: Reweighter.h
//Brief: A Reweighter changes the CV model into a different model using just a multiplicative
//       constant.  All vertical systematics are implemented by taking ratios to such weights.
//       Some Reweighters are mutually exclusive, and others are only needed for specific systematics.
//Author: Andrew Olivier aolivier@ur.rochester.edu

#ifndef PLOTUTILS_REWEIGHTER_H
#define PLOTUTILS_REWEIGHTER_H


#if __cplusplus < 201103L
  #define override
#endif

//c++ includes
#include <string>
#include <vector>

namespace PlotUtils
{
  //Ugly hack to allow us to include both Cut.h and Reweighter.h (from MAT-MINERvA) in the same file or either one individually.  The right solution is usually to create ano
  //ther file, but I don't want to confuse people with a file that just creates an empty structure :(
  #ifndef PLOTUTILS_DETAIL_EMPTY
  #define PLOTUTILS_DETAIL_EMPTY
  namespace detail
  {
    struct empty {};
  }
  #endif //PLOTUTILS_DETAIL_EMPTY

  template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
  class Reweighter
  {
    public:
      Reweighter() = default;
      virtual ~Reweighter() = default;

      virtual double GetWeight(const UNIVERSE& univ, const EVENT& event) const = 0;
      virtual std::string GetName() const = 0;

      virtual bool DependsReco() const = 0;
      /*virtual bool DependsTruth() const;*/ //Not needed as of time of writing.

      virtual bool IsCompatible(const Reweighter& /*other*/) const { return true; }
      virtual std::vector<UNIVERSE*> GetRequiredUniverses() const { return std::vector<UNIVERSE*>{}; }
  };
}

#endif //PLOTUTILS_REWEIGHTER_H
