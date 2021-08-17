//File: MINOSEfficiencyReweighter.h
//Brief: A Reweighter that changes the simulated FS particle model.
//Author: Andrew Olivier aolivier@ur.rochester.edu

#ifndef PLOTUTILS_MINOSEfficiencyREWEIGHTER_H
#define PLOTUTILS_MINOSEfficiencyREWEIGHTER_H

//PlotUtils includes
#include "utilities/NSFDefaults.h"
#include "universes/GenieSystematics.cxx" //IsNonResPi()
#include "universes/MnvTuneSystematics.cxx" //IsCCRes()
#include "weighters/MinosMuonEfficiencyCorrection.h"
#include "utilities/MnvNormalization.h"

//Reweighter includes
#include "weighters/Reweighter.h"

namespace PlotUtils
{
  template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
  class MINOSEfficiencyReweighter: public Reweighter<UNIVERSE, EVENT>
  {
    public:
      MINOSEfficiencyReweighter(): Reweighter<UNIVERSE, EVENT>()
      {
      }

      virtual ~MINOSEfficiencyReweighter() = default;

      double GetWeight(const UNIVERSE& univ, const EVENT& /*event*/) const override
      {
        if(univ.IsTruth()) return 1.; // No efficiency reweigting for truth events
        if(!univ.isFHC() && !univ.isRHC()) return 1.; // LE or nonstandard playlist.  Don't have away of dealing with this yet

        if (univ.IsPlaylistME(univ.GetPlaylist())) {  // The correction factors are different
                                      // between ME and LE
        double pmu = univ.GetPmuMinos() / 1000;  // GetCorrection expects GeV
        return PlotUtils::MinosMuonEfficiencyCorrection::Get(univ.isFHC()).GetCorrection(pmu, univ.GetBatchPOT(), univ.isFHC());
        } else {                       // Assume if not ME, then it's LE
          double pmu = univ.GetPmuMinos();  // MnVnormalizer GetCorrection expects MeV
        #ifndef __CINT__
          static PlotUtils::MnvNormalizer mnvNormalizer =
              PlotUtils::MnvNormalizer("Eroica", univ.GetPlaylist());
        #endif // __CINT__
          return mnvNormalizer.GetCorrection(pmu);
        }
      }

      std::string GetName() const override { return "MINOSEfficiency"; }
      bool DependsReco() const override { return true; }
  };
}

#endif //PLOTUTILS_MINOSEfficiencyREWEIGHTER_H
