//File: LowRecoil2p2hReweighter.h
//Brief: A Reweighter that changes the simulated FS particle model.
//Author: Andrew Olivier aolivier@ur.rochester.edu

#ifndef PLOTUTILS_LowRecoil2p2hREWEIGHTER_H
#define PLOTUTILS_LowRecoil2p2hREWEIGHTER_H

//PlotUtils includes
#include "utilities/NSFDefaults.h"
#include "universes/MnvTuneSystematics.h"

//Reweighter includes
#include "weighters/Reweighter.h"

namespace PlotUtils
{
  template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
  class LowRecoil2p2hReweighter: public Reweighter<UNIVERSE, EVENT>
  {
    public:
      LowRecoil2p2hReweighter(const int mode = 0): Reweighter<UNIVERSE, EVENT>(), fMode(mode)
      {
      }

      virtual ~LowRecoil2p2hReweighter() = default;

      double GetWeight(const UNIVERSE& univ, const EVENT& /*event*/) const override
      {
        //variation 0 is the CV
        return PlotUtils::GetLowRecoil2p2hWeight(univ, univ.Getq0True() / 1000 /* GeV */,
                                           univ.Getq3True() / 1000 /* GeV */,
                                           fMode);
      }

      std::string GetName() const override { return "LowRecoil2p2hTune"; }
      bool DependsReco() const override { return false; }

    private:
      int fMode; //Turns on different 2p2h variations.  0 is our CV for MnvTunev1.  1 and 2 weight only np and nn/pp events.  3 weights up QE events instead of 2p2h.
  };
}

#endif //PLOTUTILS_LowRecoil2p2hREWEIGHTER_H
