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
      GenericVerticalUniverse(typename BASE_UNIVERSE::config_t& config, std::unique_ptr<PlotUtils::Reweighter<BASE_UNIVERSE, EVENT>>&& weighter): BASE_UNIVERSE(config), fReweight(std::move(weighter))
      {
      }

      ~GenericVerticalUniverse() override = default;

      double GetWeightRatioToCV() const override
      {
        const EVENT empty; //TODO: Update GetWeightRatioToCV() to work with EVENT?  For now, this harness just won't work with Reweighters that rely on a specific EVENT format.
        return fReweight->GetWeight(*this, empty);
      }

      std::string ShortName() const override { return fReweight->GetName(); }
      bool IsVerticalOnly() const override { return true; } //By definition

    private:
      std::unique_ptr<PlotUtils::Reweighter<BASE_UNIVERSE, EVENT>> fReweight;
  };
}
