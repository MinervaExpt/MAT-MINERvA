//File: NeutronInelasticReweighter.h
//Brief: A Reweighter that changes MINERvA's neutron inelastic cross sections for several channels into the inelastic cross sections from low energy neutron data.
//Author: Andrew Olivier aolivier@ur.rochester.edu

//MAT-MINERvA includes
#include "utilities/TargetUtils.h"

//ROOT includes
#include "TH1D.h"
#include "TGraph.h"
#include "TSpline3.h"
#include "TF1.h"
#include "TFile.h"

//c++ includes
#include <vector>
#include <map>
#include <memory>
#include <numeric>

template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
class NeutronInelasticReweighter: public PlotUtils::Reweighter<UNIVERSE, EVENT>
{
  public:
    NeutronInelasticReweighter(const std::string& oldTotalInel, const std::string& newTotalInel, const std::map<std::string, std::vector<int>>& fileNameToFS): fTotalInelastic(oldTotalInel, newTotalInel, {}), fGeometry()
    {
      //Load fKinENormalization from a file.  Do this first because it can fail.
      const std::string kinEFileName = "", 
                        kinENormHistName = "";
      auto oldPwd = gDirectory;
      std::unique_ptr<TFile> kinEFile(TFile::Open(kinEFileName.c_str()));
      fKinENormalization = dynamic_cast<TH1D*>(kinEFile.Get(kinENormHistName.c_str())->Clone()); //Make a Clone() so I don't have to keep kinEFile open while the job runs.
      if(!fKinENormalization) throw std::runtime_error("Failed to load a histogram named " + kinENormHistName + " from a file named " + kinEFileName + " for neutron inelastic reweight normalization.");
      fKinENormalization.SetDirectory(nullptr); //Make sure fKineENormalization is no longer tied to its parent object's file because that file will eventually be closed.
      gDirectory = oldPwd;

      //Load interaction channels from files
      for(const auto& channel: fileNameToFS) fChannels.emplace_back(channel.first, channel.first, channel.second);

      //Create an "Other" Channel that preserves the total cross section
      fOther.fOldSigmaRatio = TF1([this](const double* x, const double* /*p*/)
                                  {
                                    double result = 1;
                                    for(const auto& channel: this->fChannels) result -= channel->fOldSigmaRatio->Eval(x[0]);
                                    return result;
                                  });
      fOther.fNewSigmaRatio = TF1([this](const double* x, const double* /*p*/)
                                  {
                                    double result = 1;
                                    for(const auto& channel: this->fChannels) result -= channel->fNewSigmaRatio->Eval(x[0]);
                                    return result;
                                  });;
      fOther.fMin = std::max_element(fChannels.begin(), fChannels.end(), [](const auto& channel) { return channel.fMin; })->fMin;
      fOther.fMax = std::min_element(fChannels.begin(), fChannels.end(), [](const auto& channel) { return channel.fMax; })->fMax;
    }

    ~NeutronInelasticReweighter() = default;

    double GetWeight(const UNIVERSE& univ, const EVENT& /*event*/) const;
    std::string GetName() const { return "NeutronInelasticExclusives"; }

    bool DependsReco() const { return false; }

  private:
    struct Channel
    {
      std::multiset<int> fInelasticChildren;

      TF1 fOldSigmaRatio; //N.B.: Using TF1s instead of TSpline3s or TGraphs so I can use Integral() function like Jeffrey did
      TF1 fNewSigmaRatio;

      double fMin;
      double fMax;

      Channel(const std::string& oldFile, const std::string& newFile, const std::vector<int> inelChildren): fInelasticChildren(inelChildren.begin(), inelChildren.end())
      {
        //TODO: oldRatioGraph might actually come from a TFile the way things are written right now
        TGraph oldRatioGraph(oldFile.c_str());
        fOldSigmaRatioSpline = TSpline3(oldRatioGraph, oldFile.substr(oldFile.rfind("/"), oldFile.find(".")-1 - oldFile.rfind("/")));
        fOldSigmaRatio = TF1([this](const double* x, const double* /*p*/){ return this->fOldSigmaRatioSpline.Eval(x[0]); });
                                                                                                                                                                         
        TGraph newRatioGraph(newFile.c_str());
        fNewSigmaRatioSpline = TSpline3(newRatioGraph, newFile.substr(newFile.rfind("/"), newFile.find(".")-1 - newFile.rfind("/")));
        fNewSigmaRatio = TF1([this](const double* x, const double* /*p*/){ return this->fNewSigmaRatioSpline.Eval(x[0]); });
                                                                                                                                                                         
        fMin = std::max(fOldSigmaRatioSpline.GetMin(), fNewSigmaRatioSpline.GetMin());
        fMax = std::min(fOldSigmaRatioSpline.GetMax(), fNewSigmaRatioSpline.GetMax());
      }

      private:
        //I think I have to keep these TSpline3 objects around because they're referenced by the TF1s :(
        TSpline3 fOldSigmaRatioSpline;
        TSpline3 fNewSigmaRatioSpline;
    };

    std::vector<Channel> fChannels; //channels that will be reweighted
    Channel fOther; //All other channels that aren't reweighted are lumped into one.  This keeps the total inelastic cross section the same.
    //Channel fTotalInelastic; //TODO: Not needed as long as I can get away with just weighting each neutron by exclusive cross section ratio

    TH1D* fKinENormalization; //Normalization to keep the overall neutrino cross section the same in kinetic energy and angle

    PlotUtils::TargetUtils fGeometry;

    double getNonInteractingWeight(const Channel& channel, const double density, const double Ti, const double Tf) const;
    double getInteractingWeight(const Channel& channel, const double /*density*/, const double Ti, const double Tf) const;

    double evalSigmaRatio(const TF1& ratioFunc, double Ti, double Tf, const double min, const double max) const;
};

template <class UNIVERSE, class EVENT>
double NeutronInelasticReweighter::GetWeight(const UNIVERSE& univ, const EVENT& /*event*/) const
{
  double weight = 1;

  constexpr double neutronMass = 939.6; //MeV/c^2
  const std::string prefix = "neutronInelasticReweight"; //Beginning of branch names for inelastic reweighting

  const int nNeutrons = univ.GetInt(prefix + "NPaths");
  const auto startEnergyPerPoint = univ.GetVecDouble(prefix + "InitialE"),
             endEnergyPerPoint = univ.GetVecDouble(prefix + "FinalE"),
             densityPerPoint = univ.GetVecDouble(prefix + "ColumnarDensity"),
             xPerPoint = univ.GetVecDouble(prefix + "PosX"),
             yPerPoint = univ.GetVecDouble(prefix + "PosY"),
             zPerPoint = univ.GetVecDouble(prefix + "PosZ");
  const auto nPointsPerNeutron = univ.GetVecInt(prefix + "NTrajPoints"),
             nInelasticChildren = univ.GetVecInt(prefix + "NInelasticChildren"),
             allInelChildren = univ.GetVecInt(prefix + "InelasticChildPDGs"),
             materialPerPoint = univ.GetVecInt(prefix + "Nuke"),
             intCodePerPoint = univ.GetVecInt(prefix + "IntCodePerSegment");

  if(nPointsPerNeutron.empty()) return 1.;

  int endPoint = -1;
  auto endInelasticChild = allInelChildren.begin();
  for(int whichNeutron = 0; whichNeutron < nNeutrons; ++whichNeutron)
  {
    const int startPoint = endPoint;
    const auto startInelasticChild = endInelasticChild;
    endPoint += nPointsPerNeutron[whichNeutron];
    endInelasticChild += nInelasticChildren[whichNeutron];

    //TODO: Only turn this on if I'm changing the total cross section.  Otherwise, it's always 1.
    //If I read MnvHadronReweight literally and keep the total cross section constant, getNonInteractingWeight()
    //is always 1.  There's also an inelastic-only reweight that does apply a non-interacting weight using the
    //inelastic cross section.  Is it correct to just leave that out?  I could make my tuples a lot simpler if so.

    //Possibly-elastic points where inelastic interaction did not happen
    //Stop before the last point because it may have ended with an inelastic interaction
    /*for(int whichPoint = startPoint; whichPoint < endPoint; ++whichPoint)
    {
      //N.B.: material of -6 seems to be a special flag Jeffrey added to denote CH scintillator as opposed to pure carbon from target 3.
      if(materialPerPoint[whichPoint] != -6) continue; //We only have data for CH scintillator

      if(fGeometry.InTracker(xPerPoint[endPoint], yPerPoint[lastPoint], zPerPoint[lastPoint]))
      {
        const double Ti = startEnergyPerPoint[whichPoint] - neutronMass,
                     Tf = endEnergyPerPoint[whichPoint] - neutronMass;
        const double density = densityPerPoint[whichPoint];

        for(const auto& channel: fChannels) weight *= getNonInteractingWeight(channel, density, Ti, Tf);
        weight *= getNonInteractingWeight(fOther, density, Ti, Tf);
      }
    }*/

    if(fGeometry.InTracker(xPerPoint[endPoint - 1], yPerPoint[lastPoint - 1], zPerPoint[lastPoint - 1]) && materialPerPoint[lastPoint - 1] == -6)
    {
      //A multi-set is a collection of numbers with a count of how many times each number came up.
      std::multiset<int> inelasticChildren(allInelChildren.begin() + startInelasticChild, allInelChildren.begin() + endInelasticChild);
      inelasticChildren.erase(22); //Ignore photons because GEANT tends to emit extra low energy photons to distribute binding energy

      const double Ti = startEnergyPerPoint[endPoint - 1],
                   Tf = endEnergyPerPoint[endPoint - 1],
                   density = densityPerPoint[endPoint - 1];
      const int intCode = intCodePerPoint[endPoint - 1];

      //Inelastic interactions end any TG4Trajectory.  Figure out whether this
      //trajectory ended with an inelastic interaction.  If so, is it one of
      //the channels I'm reweighting?

      if(intCode == 1 || intCode == 4) //If there was an inelastic interaction
      {
        const auto foundChannel = std::find(fChannels.begin(), fChannels.end(),
                                            [&inelasticChildren](const auto& channel)
                                            {
                                              return channel.fInelasticChildren == inelasticChildren;
                                            });
        if(foundChannel != fChannels.end()) weight *= getInteractingWeight(*foundChannel, density, Ti, Tf);
        else getInteractingWeight(fOther, density, Ti, Tf); //Other channel
      }
      /*else //If this trajectory ended by some process other than an inelastic interaction
      {
        //TODO: Only turn this on if I'm changing the total cross section.  Otherwise, it's always 1.
        //TODO: I think this is equivalent to reweighting based on the total cross section because I'm multiplying exponentials with the same
        //      coefficients.  But, I may pick up a larger roundoff error this way.  Do I need the total inelastic cross section to use
        //      getNonInteractingWeight() anyway?
        for(const auto& channel: fChannels) weight *= getNonInteractingWeight(channel, density, Ti, Tf);
        weight *= getNonInteractingWeight(fOther, density, Ti, Tf);
      }*/
    } //If last point is in the tracker and CH scintillator

    //Divide by a kinematics-dependent normalization factor to keep the total neutrino cross section constant.
    //TODO: Aaron only does this for the leading particle in the original MnvHadronReweight
    //TODO: Aaron does this using FS particle branches because he only cares about FS particles.  Do I have momentum components for all neutrons?
    //      If not, I'm tempted to try just reweighting in neutron KE first.
    //      Nope, I don't have neutron direction.  Trying neutron KE until I see that it's a problem.
    weight /= fKinENormalization->GetBinContent(fKinENormalization->FindBin(startEnergyPerPoint[startPoint] - neutronMass));
  } //For each neutron

  return weight;
}

template <class UNIVERSE, class EVENT>
double NeutronInelasticReweighter<UNIVERSE, EVENT>::getNonInteractingWeight(const Channel& channel, const double density, const double Ti, const double Tf) const
{
  //TODO: Need distance and density correction.  Multiply density by 4.626e22 for MINERvA's CH scintillator
  //TODO: Multiply by total inelastic cross section.
  return exp(-1.0 * density * distance * (evalSigmaRatio(channel.fOldSigmaRatio, Ti, Tf, channel.fMin, channel.fMax) - evalSigmaRatio(channel.fNewSigmaRatio, Ti, Tf, channel.fMin, channel.fMax)));
}

template <class UNIVERSE, class EVENT>
double NeutronInelasticReweighter<UNIVERSE, EVENT>::getInteractingWeight(const Channel& channel, const double /*density*/, const double Ti, const double Tf) const
{
  //TODO: Need distance and density correction
  //I don't need to reweight based on the total cross section because I'm implicitly keeping it the same.
  /*const double denom = 1. - exp(-1. * density * distance * evalSigmaRatio(fTotalInelastic.fOldSigmaRatio, Ti, Tf, fTotalInelastic.fMin, fTotalInelastic.fMax));
  if(denom <= 0) return 0;
  const double num = 1. - exp(-1. * density * distance * evalSigmaRatio(fTotalInelastic.fNewSigmaRatio, Ti, Tf, fTotalInelastic.fMin, fTotalInelastic.fMax));*/

  const double a = evalSigmaRatio(channel.fNewSigmaRatio, Ti, Tf, channel.fMin, channel.fMax);
  const double b = evalSigmaRatio(channel.fOldSigmaRatio, Ti, Tf, channel.fMin, channel.fMax);
  //return num / denom * a / b;
  return a / b;
}

//Adapt to graph evaluation pitfalls
template <class UNIVERSE, class EVENT>
double NeutronInelasticReweighter::evalSigmaRatio(const TF1& ratioFunc, double Ti, double Tf, const double min, const double max) const
{
  //Prefer rounding into the range where we have data over interpolating off the end of a spline
  //"clamp" Ti and Tf to min/max of ratioFunc
  Ti = std::min(Ti, max);
  Tf = std::min(Tf, max);

  //Some strange "linear interpolation towards 0" that Jeffrey does.  He also comments that this never happens in MnvHadronReweight because
  //the "HD neutron cross section" goes down to 1 MeV.
  //N.B.: ratioFunc wraps over a cubic spline to data
  if(Ti < min) Ti *= ratioFunc.Eval(Ti)/min;
  if(Tf < min) Tf *= ratioFunc.Eval(Tf)/min;

  if(Ti == Tf) return ratioFunc.Eval(Ti);
  return ratioFunc.Integral(Ti, Tf)/(Tf - Ti); //TF1::Integral() is supposedly a Gaussian quadrature algorithm in some cases
}
