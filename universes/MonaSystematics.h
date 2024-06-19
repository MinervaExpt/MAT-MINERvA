#ifndef MONASystematics_h
#define MONASystematics_h

//==============================================================================
// Get Several standard MINERvA systematics
//==============================================================================

#include "PlotUtils/GenericVerticalSystematic.h"
#include "PlotUtils/NeutronInelasticReweighter.h"

namespace PlotUtils
{
  std::map<std::string,std::vector<int>> MonaMapDefault = {{"nGamma",{1000060120, 2112}},
                {"threeAlpha",{1000020040, 100002040, 100002040, 2112}},
                {"Bnp",{1000050110, 2112, 2212}}};

  template <class T>
  std::map< std::string, std::vector<T*> > GetMonaSystematicMap(typename T::config_t chain)
  {
    // return map
    std::map< std::string, std::vector<T*> > ret;

    ret["NeutronInelasticsReweight"].push_back(new PlotUtils::GenericVerticalUniverse<T, PlotUtils::detail::empty>(chain, std::unique_ptr<PlotUtils::Reweighter<T, PlotUtils::detail::empty>>(new NeutronInelasticReweighter<T>(MonaMapDefault)), 1.0));
    ret["NeutronInelasticsReweight"].push_back(new PlotUtils::GenericVerticalUniverse<T, PlotUtils::detail::empty>(chain, std::unique_ptr<PlotUtils::Reweighter<T, PlotUtils::detail::empty>>(new NeutronInelasticReweighter<T>(MonaMapDefault)), -1.0));

    std::cout << "MoNA systematics created." << std::endl;

    return ret;
  }
}
#endif  // MONASystematics_h
