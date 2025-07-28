//File: GenericVerticalSystematic.h
//Brief: Wrapper over a PlotUtils::Reweighter<> that turns a CV reweight into a systematic Universe.
//Author: Andrew Olivier aolivier@ur.rochester.edu

//MAT-MINERvA includes
#include "weighters/Reweighter.h"

namespace PlotUtils
{
  template <class BASE_UNIVERSE, class EVENT>
  class GenericVerticalUniverse: public BASE_UNIVERSE
  {
    public:
      GenericVerticalUniverse(typename BASE_UNIVERSE::config_t& config, std::unique_ptr<PlotUtils::Reweighter<BASE_UNIVERSE, EVENT>>&& weighter, double scale = 1.0, std::string name = "") : BASE_UNIVERSE(config), fReweight(std::move(weighter)), fScaleAboutCV(scale), fNameExt(name)
      {
      }

      ~GenericVerticalUniverse() override = default;

      double GetWeightRatioToCV() const override
      {
        const EVENT empty; //TODO: Update GetWeightRatioToCV() to work with EVENT?  For now, this harness just won't work with Reweighters that rely on a specific EVENT format.
        // The following modifications were added to allow for the difference between this reweighted universe and the CV to be treated more easily as a systematic shift. Arbitrary scaling magnitude allows one to determine how many sigma the difference should be consisdered as. In the later function which grabs the systematic map for these reweight universes. -David L. July 29, 2022
        if (fScaleAboutCV == 0)
          return 1.0;
        else if (fScaleAboutCV > 0)
          return fScaleAboutCV * fReweight->GetWeight(*this, empty);
        else 
          return (2.0 + fScaleAboutCV * fReweight->GetWeight(*this, empty));
      }

      std::string ShortName() const override { return ((fNameExt != "") ? fReweight->GetName() + "_" + fNameExt : fReweight->GetName()); }
      bool IsVerticalOnly() const override { return true; } //By definition

    private:
      std::unique_ptr<PlotUtils::Reweighter<BASE_UNIVERSE, EVENT>> fReweight;
      double fScaleAboutCV;
      std::string fNameExt;
  };
}
