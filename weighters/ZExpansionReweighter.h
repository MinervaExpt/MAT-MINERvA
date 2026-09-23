// File: ZExpansionReweighter.h
// Brief: A Reweighter to change from dipole formalism to z-expansion axial form factor from the 2016 paper fits to bubble chamber data (used in MnvTune vX.4.Z)
// Author: Noah Harvey Vaughan, vaughann@oregonstate.edu, nhvaughan on github

#ifndef PLOTUTILS_ZExpansionREWEIGHTER_H
#define PLOTUTILS_ZExpansionREWEIGHTER_H

// Reweighter includes
#include "weighters/Reweighter.h"
#include "weighters/weightZExp.h"

namespace PlotUtils {
template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
class ZExpansionReweighter : public Reweighter<UNIVERSE, EVENT> {
   public:
    ZExpansionReweighter() : Reweighter<UNIVERSE, EVENT>(), fCalculator(std::string(std::getenv("MPARAMFILESROOT")) + "/data/Reweight/Z_Expansion_Reweight_v2126.root") {
    }

    virtual ~ZExpansionReweighter() = default;

    double GetWeight(const UNIVERSE& univ, const EVENT& /*event*/) const override {
        // Needs to be CCQE
        if (univ.GetInt("mc_intType") != 1) return 1.0;

        // Not on hydrogren? Based off Dan's implementation in old framework, not sure if this makes sense...
        // if (univ.GetInt("mc_targetZ")<6) return 1.0;
        
        // variation 0 is the CV
        return fCalculator.getWeight(univ.GetQ2True() * 1e-6 /*GeV^2*/);
        // Systematics are added when you set UseZExpansionReweight on MinervaUniverse, and it will make 100 universes when calling GetGenieSystematics or GetGenieSystematicsMap
    }

    std::string GetName() const override { return "ZExpansion"; }
    bool DependsReco() const override { return false; }

    // std::vector<UNIVERSE*> GetRequiredUniverses() const override {
    //     return std::vector<UNIVERSE*>{};  // TODO: Return the zexp universes here?
    // }
   private:
    mutable PlotUtils::weightZExp fCalculator;
};
}  // namespace PlotUtils

#endif  // PLOTUTILS_ZExpansionREWEIGHTER_H
