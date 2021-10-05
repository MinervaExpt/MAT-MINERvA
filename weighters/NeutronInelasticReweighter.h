//File: NeutronInelasticReweighter.h
//Brief: A Reweighter that changes MINERvA's neutron inelastic cross sections for several channels into the inelastic cross sections from low energy neutron data.
//Author: Andrew Olivier aolivier@ur.rochester.edu

//MAT-MINERvA includes
#include "utilities/TargetUtils.h"

//ROOT includes
#include "TH1D.h"
#include "TGraph.h"
#include "TSpline.h"
#include "TF1.h"
#include "TFile.h"

//c++ includes
#include <vector>
#include <map>
#include <memory>
#include <numeric>

namespace
{
  constexpr double scintDensityToNucleons = 4.626e22;
}

template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
class NeutronInelasticReweighter: public PlotUtils::Reweighter<UNIVERSE, EVENT>
{
  public:
    NeutronInelasticReweighter(const std::map<std::string, std::vector<int>>& fileNameToFS): fTotalInelastic("inelastic", {}), fGeometry()
    {
      fChannels.reserve(fileNameToFS.size()); //If I don't use this, the program will often crash.  std::vector::emplace_back() will have to
                                              //reallocate memory many times.  When it does that, it copies the old Channels is made and then
                                              //deletes the originals.  But the copied TF1s hold lambda functions that still point at the
                                              //original (now deleted) Channels.

      //Load fKinENormalization from a file.  Do this first because it can fail.
      const std::string kinEFileName = "", 
                        kinENormHistName = "";
      auto oldPwd = gDirectory;
      /*std::unique_ptr<TFile> kinEFile(TFile::Open(kinEFileName.c_str()));
      fKinENormalization = dynamic_cast<TH1D*>(kinEFile->Get(kinENormHistName.c_str())->Clone()); //Make a Clone() so I don't have to keep kinEFile open while the job runs.
      if(!fKinENormalization) throw std::runtime_error("Failed to load a histogram named " + kinENormHistName + " from a file named " + kinEFileName + " for neutron inelastic reweight normalization.");
      fKinENormalization->SetDirectory(nullptr); //Make sure fKineENormalization is no longer tied to its parent object's file because that file will eventually be closed.*/

      //Load total elastic cross section from a file
      {
        std::unique_ptr<TFile> totalElasticFile(TFile::Open("cross_section.root"));
        assert(totalElasticFile);
        auto totalElasticGraph = dynamic_cast<TGraph*>(totalElasticFile->Get("elastic"));
        assert(totalElasticGraph);
        fTotalElasticSpline = TSpline3("elastic", totalElasticGraph);
      }

      //Load interaction channels from files
      for(const auto& channel: fileNameToFS) fChannels.emplace_back(channel.first, channel.second);

      //Create an "Other" Channel that preserves the total cross section
      /*fOther.fMin = std::max_element(fChannels.begin(), fChannels.end(), [](const auto& lhs, const auto& rhs) { return lhs.fMin < rhs.fMin; })->fMin;
      fOther.fMax = std::min_element(fChannels.begin(), fChannels.end(), [](const auto& lhs, const auto& rhs) { return lhs.fMax < rhs.fMax; })->fMax;*/

      fLowestMinKE = std::min_element(fChannels.begin(), fChannels.end(), [](const auto& lhs, const auto& rhs) { return lhs.fMin < rhs.fMin; })->fMin;
      fHighestMaxKE = std::max_element(fChannels.begin(), fChannels.end(), [](const auto& lhs, const auto& rhs) { return lhs.fMax < rhs.fMax; })->fMax;

      fTotalElastic = TF1("total", [this](double* x, double* /*p*/) { return fTotalElasticSpline.Eval(x[0]); }, fLowestMinKE, fHighestMaxKE, 0);

      /*fOther.fOldSigmaRatio = TF1("old", [this](double* x, double* p)
                                         {
                                           double result = 1;
                                           for(const auto& channel: this->fChannels) result -= channel.fOldSigmaRatio.Eval(x[0]);
                                           return result;
                                         }, fOther.fMin, fOther.fMax, 0);
      fOther.fNewSigmaRatio = TF1("new", [this](double* x, double* p)
                                         {
                                           double result = 1;
                                           for(const auto& channel: this->fChannels) result -= channel.fNewSigmaRatio.Eval(x[0]);
                                           return result;
                                         }, fOther.fMin, fOther.fMax, 0);*/

      gDirectory = oldPwd;
    }

    ~NeutronInelasticReweighter() = default;

    double GetWeight(const UNIVERSE& univ, const EVENT& /*event*/) const;
    std::string GetName() const { return "NeutronInelasticExclusives"; }

    bool DependsReco() const { return false; }

  private:
    struct Channel
    {
      std::multiset<int> fInelasticChildren;

      mutable TF1 fOldSigmaRatio; //N.B.: Using TF1s instead of TSpline3s or TGraphs so I can use Integral() function like Jeffrey did
      mutable TF1 fNewSigmaRatio;

      double fMin;
      double fMax;

      Channel(): fInelasticChildren(), fOldSigmaRatio(), fNewSigmaRatio(), fMin(), fMax(), fOldSigmaRatioSpline(), fNewSigmaRatioSpline()
      {
      }

      Channel(const std::string& channelName, const std::vector<int> inelChildren): fInelasticChildren(inelChildren.begin(), inelChildren.end())
      {
        const std::string oldFileName = "cross_section.root"; //TODO: Use PlotUtilsROOT or something to find this file
        std::unique_ptr<TFile> oldGraphFile(TFile::Open(oldFileName.c_str()));
        if(!oldGraphFile) throw std::runtime_error("Failed to open a file named " + oldFileName + " for a GEANT cross section graph in InelasticNeutronReweighter::Channel.");
        //TGraph oldRatioGraph(oldFile.c_str());
        auto oldRatioGraph = dynamic_cast<TGraph*>(oldGraphFile->Get(channelName.c_str()));
        if(!oldRatioGraph) throw std::runtime_error("Failed to load a TGraph named " + channelName + " from a file named " + oldFileName + " for NeutronInelasticReweighter::Channel");
        fOldSigmaRatioSpline = TSpline3(channelName.c_str(), oldRatioGraph);

        TGraph newRatioGraph((channelName + ".csv").c_str());
        fNewSigmaRatioSpline = TSpline3(channelName.c_str(), &newRatioGraph);

        fMin = std::max(fOldSigmaRatioSpline.GetXmin(), fNewSigmaRatioSpline.GetXmin());
        fMax = std::min(fOldSigmaRatioSpline.GetXmax(), fNewSigmaRatioSpline.GetXmax());

        //fOldSigmaRatio = TF1("old", [this](double* x, double* /*p*/){ return fOldSigmaRatioSpline.Eval(x[0]); }, fMin, fMax, 0);
        //fNewSigmaRatio = TF1("new", [this](double* x, double* /*p*/){ return fNewSigmaRatioSpline.Eval(x[0]); }, fMin, fMax, 0);

        fOldSigmaRatio = TF1("old", this, &Channel::evalOldSpline, fMin, fMax, 0);
        fNewSigmaRatio = TF1("new", this, &Channel::evalNewSpline, fMin, fMax, 0);
      }

      double evalOldSpline(double* x, double* /*p*/) const { return fOldSigmaRatioSpline.Eval(x[0]); }
      double evalNewSpline(double* x, double* /*p*/) const { return fNewSigmaRatioSpline.Eval(x[0]); }

      private:
        //I think I have to keep these TSpline3 objects around because they're referenced by the TF1s :(
        TSpline3 fOldSigmaRatioSpline;
        TSpline3 fNewSigmaRatioSpline;
    };

    std::vector<Channel> fChannels; //channels that will be reweighted
    //Channel fOther; //All other channels that aren't reweighted are lumped into one.  This keeps the total inelastic cross section the same.
    Channel fTotalInelastic; //Wouldn't be needed if I could get away with just weighting each neutron by exclusive cross section ratio
    mutable TF1 fTotalElastic; //I really need the total cross section.  But the inelastic changes for MoNA while I'm assuming that the elastic doesn't.
    TSpline3 fTotalElasticSpline; //This is referenced by a TF1, so I need to keep it alive as long as fTotalElastic is alive

    //KE range which at least some channel covers.  Any neutrons outside of this range just get a weight of 1.
    double fLowestMinKE;
    double fHighestMaxKE;

    TH1D* fKinENormalization; //Normalization to keep the overall neutrino cross section the same in kinetic energy and angle

    PlotUtils::TargetUtils fGeometry;

    double getNonInteractingWeight(const Channel& channel, const double density, const double Ti, const double Tf) const;
    double getInteractingWeight(const Channel& channel, const double /*density*/, const double Ti, const double Tf) const;

    double evalSigmaRatio(TF1& ratioFunc, double Ti, double Tf, const double min, const double max) const;
};

template <class UNIVERSE, class EVENT>
double NeutronInelasticReweighter<UNIVERSE, EVENT>::GetWeight(const UNIVERSE& univ, const EVENT& /*event*/) const
{
  double weight = 1;

  constexpr double neutronMass = 939.6; //MeV/c^2
  const std::string prefix = "truth_neutronInelasticReweight"; //Beginning of branch names for inelastic reweighting

  const int nNeutrons = univ.GetInt((prefix + "NPaths").c_str());
  const auto startEnergyPerPoint = univ.GetVecDouble((prefix + "InitialE").c_str()),
             endEnergyPerPoint = univ.GetVecDouble((prefix + "FinalE").c_str()),
             densityPerPoint = univ.GetVecDouble((prefix + "ColumnarDensity").c_str()),
             xPerPoint = univ.GetVecDouble((prefix + "PosX").c_str()),
             yPerPoint = univ.GetVecDouble((prefix + "PosY").c_str()),
             zPerPoint = univ.GetVecDouble((prefix + "PosZ").c_str());
  const auto nPointsPerNeutron = univ.GetVecInt((prefix + "NTrajPointsSaved").c_str()),
             nInelasticChildren = univ.GetVecInt((prefix + "NInelasticChildren").c_str()),
             allInelChildren = univ.GetVecInt((prefix + "InelasticChildPDGs").c_str()),
             materialPerPoint = univ.GetVecInt((prefix + "Nuke").c_str()),
             intCodePerPoint = univ.GetVecInt((prefix + "IntCodePerSegment").c_str());

  if(nPointsPerNeutron.empty()) return 1.;

  int endPoint = 0,
      endInelasticChild = 0;
  for(int whichNeutron = 0; whichNeutron < nNeutrons; ++whichNeutron)
  {
    const int startPoint = endPoint,
              startInelasticChild = endInelasticChild;
    endPoint += nPointsPerNeutron[whichNeutron];
    endInelasticChild += nInelasticChildren[whichNeutron];

    //TODO: Only turn this on if I'm changing the total cross section.  Otherwise, it's always 1.
    //If I read MnvHadronReweight literally and keep the total cross section constant, getNonInteractingWeight()
    //is always 1.  There's also an inelastic-only reweight that does apply a non-interacting weight using the
    //inelastic cross section.  Is it correct to just leave that out?  I could make my tuples a lot simpler if so.

    //Possibly-elastic points where inelastic interaction did not happen
    //Stop before the last point because it may have ended with an inelastic interaction
    for(int whichPoint = startPoint; whichPoint < endPoint; ++whichPoint)
    {
      //N.B.: material of -6 seems to be a special flag Jeffrey added to denote CH scintillator as opposed to pure carbon from target 3.
      if(materialPerPoint[whichPoint] != -6) continue; //We only have data for CH scintillator

      if(fGeometry.InTracker(xPerPoint[endPoint], yPerPoint[endPoint], zPerPoint[endPoint]))
      {
        const double Ti = startEnergyPerPoint[whichPoint] - neutronMass,
                     Tf = endEnergyPerPoint[whichPoint] - neutronMass;
        const double density = densityPerPoint[whichPoint];

        for(const auto& channel: fChannels) weight *= getNonInteractingWeight(channel, density, Ti, Tf);
        //weight *= getNonInteractingWeight(fOther, density, Ti, Tf);
      }
    }

    if(startPoint != endPoint && fGeometry.InTracker(/*xPerPoint[endPoint - 1]*/ xPerPoint.at(endPoint - 1), yPerPoint[endPoint - 1], zPerPoint[endPoint - 1]) && materialPerPoint[endPoint - 1] == -6)
    {
      //A multi-set is a collection of numbers with a count of how many times each number came up.
      std::multiset<int> inelasticChildren(allInelChildren.begin() + startInelasticChild, allInelChildren.begin() + endInelasticChild);
      inelasticChildren.erase(22); //Ignore photons because GEANT tends to emit extra low energy photons to distribute binding energy

      //Break up very short-lived nuclei
      if(inelasticChildren.count(1000040080))
      {
        inelasticChildren.insert(1000020040);
        inelasticChildren.insert(1000020040);
        inelasticChildren.erase(1000040080);
      }

      const double Ti = startEnergyPerPoint[endPoint - 1] - neutronMass,
                   Tf = endEnergyPerPoint[endPoint - 1] - neutronMass,
                   density = densityPerPoint[endPoint - 1];
      const int intCode = intCodePerPoint[endPoint - 1];

      //If a particle is outside the range covered by the MoNA paper's data, that's OK.
      //Don't apply a weight for that particle.
      if(Ti < fLowestMinKE || Ti > fHighestMaxKE) continue;

      //Inelastic interactions end any TG4Trajectory.  Figure out whether this
      //trajectory ended with an inelastic interaction.  If so, is it one of
      //the channels I'm reweighting?

      if(intCode == 1 || intCode == 4) //If there was an inelastic interaction
      {
        const auto foundChannel = std::find_if(fChannels.begin(), fChannels.end(),
                                               [&inelasticChildren](const auto& channel)
                                               {
                                                 return channel.fInelasticChildren == inelasticChildren;
                                               });
        if(foundChannel != fChannels.end()) weight *= getInteractingWeight(*foundChannel, density, Ti, Tf);
        //else getInteractingWeight(fOther, density, Ti, Tf); //Other channel
      }
      else //If this trajectory ended by some process other than an inelastic interaction
      {
        //TODO: Only turn this on if I'm changing the total cross section.  Otherwise, it's always 1.
        //TODO: I think this is equivalent to reweighting based on the total cross section because I'm multiplying exponentials with the same
        //      coefficients.  But, I may pick up a larger roundoff error this way.  Do I need the total inelastic cross section to use
        //      getNonInteractingWeight() anyway?
        for(const auto& channel: fChannels) weight *= getNonInteractingWeight(channel, density, Ti, Tf);
        //weight *= getNonInteractingWeight(fOther, density, Ti, Tf);
      }
    } //If last point is in the tracker and CH scintillator

    //Divide by a kinematics-dependent normalization factor to keep the total neutrino cross section constant.
    //TODO: Aaron only does this for the leading particle in the original MnvHadronReweight
    //TODO: Aaron does this using FS particle branches because he only cares about FS particles.  Do I have momentum components for all neutrons?
    //      If not, I'm tempted to try just reweighting in neutron KE first.
    //      Nope, I don't have neutron direction.  Trying neutron KE until I see that it's a problem.
    //weight /= fKinENormalization->GetBinContent(fKinENormalization->FindBin(startEnergyPerPoint[startPoint] - neutronMass));
  } //For each neutron

  return weight;
}

template <class UNIVERSE, class EVENT>
double NeutronInelasticReweighter<UNIVERSE, EVENT>::getNonInteractingWeight(const Channel& channel, const double density, const double Ti, const double Tf) const
{
  //TODO: Use total cross section, not just total inelastic.  I don't think I can get out of involving the elastic cross section here.
  //      Oh, wait, I just end up subtracting the total elastic part if it's the same for both.  Nevermind?
  return exp(-1.0 * density * scintDensityToNucleons * (evalSigmaRatio(channel.fOldSigmaRatio, Ti, Tf, channel.fMin, channel.fMax) - evalSigmaRatio(channel.fNewSigmaRatio, Ti, Tf, channel.fMin, channel.fMax)));
}

template <class UNIVERSE, class EVENT>
double NeutronInelasticReweighter<UNIVERSE, EVENT>::getInteractingWeight(const Channel& channel, const double density, const double Ti, const double Tf) const
{
  //I don't need to reweight based on the total cross section because I'm implicitly keeping it the same.
  const double totalElastic = evalSigmaRatio(fTotalElastic, Ti, Tf, fLowestMinKE, fHighestMaxKE);
  const double denom = 1. - exp(-1. * density * scintDensityToNucleons * (evalSigmaRatio(fTotalInelastic.fOldSigmaRatio, Ti, Tf, fTotalInelastic.fMin, fTotalInelastic.fMax) + totalElastic));
  if(denom <= 0) return 0;
  const double num = 1. - exp(-1. * density * scintDensityToNucleons * (evalSigmaRatio(fTotalInelastic.fNewSigmaRatio, Ti, Tf, fTotalInelastic.fMin, fTotalInelastic.fMax) + totalElastic));

  const double a = evalSigmaRatio(channel.fNewSigmaRatio, Ti, Tf, channel.fMin, channel.fMax);
  const double b = evalSigmaRatio(channel.fOldSigmaRatio, Ti, Tf, channel.fMin, channel.fMax);
  return num / denom * a / b;
  //return a / b; //Case for when not changing the total inelastic cross section
}

//Adapt to graph evaluation pitfalls
template <class UNIVERSE, class EVENT>
double NeutronInelasticReweighter<UNIVERSE, EVENT>::evalSigmaRatio(TF1& ratioFunc, double Ti, double Tf, const double min, const double max) const
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
