//File: GeantNeutronCVReweighter.h
//Brief: A Reweighter that changes the inelastic interaction cross section for FS neutrons.
//Author: Andrew Olivier aolivier@ur.rochester.edu

#ifndef PLOTUTILS_GeantNeutronCVREWEIGHTER_H
#define PLOTUTILS_GeantNeutronCVREWEIGHTER_H

//PlotUtils includes
#include "utilities/NSFDefaults.h"
#include "universes/GeantHadronSystematics.h"

//Reweighter includes
#include "weighters/Reweighter.h"

namespace PlotUtils
{
  template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
  class GeantNeutronCVReweighter: public Reweighter<UNIVERSE, EVENT>
  {
    public:
      GeantNeutronCVReweighter(): Reweighter<UNIVERSE, EVENT>()
      {
      }

      virtual ~GeantNeutronCVReweighter() = default;

      double GetWeight(const UNIVERSE& univ, const EVENT& /*event*/) const override
      {
	//Removing the following line after confirmation that it should not be the standard practice for this weight, and this should be applied to the efficiency denominator to combat biasing from the reweight. -David L. 10/09/2024
        //if(univ.IsTruth()) return 1; // No efficiency reweighting for truth events

        univ.SetupMHRWeighter(); //TODO: do this once in the constructor
        return PlotUtils::weight_hadron<PlotUtils::TreeWrapper*>().reweightNeutronCV(univ);
      }

      std::string GetName() const override { return "GeantNeutronCV"; }
      bool DependsReco() const override { return false; }
  };
}

#endif //PLOTUTILS_GeantNeutronCVREWEIGHTER_H
