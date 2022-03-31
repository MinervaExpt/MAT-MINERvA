//File: NeutronInelasticReweighter.h
//Brief: A Reweighter that changes MINERvA's neutron inelastic cross sections for several channels into the inelastic cross sections from low energy neutron data.
//Author: Andrew Olivier aolivier@ur.rochester.edu

//MAT includes
#include "PlotUtils/ErrorHandler.h" //For ROOT::exception in case it's being used to react to non-existent files.

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
  constexpr double scintDensityToNucleons = 4.626e22 * 1e-27; //Nucleons per cubic millimeter times cm^2 per millibarn
  constexpr double neutronMass = 939.6; //MeV/c^2

  double findFirstNonZeroPointX(const TGraph& graph)
  {
    const int nPoints = graph.GetN();
    double x, y;
    for(int whichPoint = 0; whichPoint < nPoints; ++whichPoint)
    {
      graph.GetPoint(whichPoint, x, y);
      if(y > 0) return x;
    }

    return nPoints;
  }

  double findLastNonZeroPointX(const TGraph& graph)
  {
    const int nPoints = graph.GetN();
    double x, y;
    for(int whichPoint = nPoints-1; whichPoint >= 0; --whichPoint)
    {
      graph.GetPoint(whichPoint, x, y);
      if(y > 0) return x;
    }

    return 0;
  }

  //Analytical integral for a TSpline3.  This should be much faster and even more accurate than the TF1::Integral() that Jeffrey used in MnvHadronReweight.
  double antiderivative(TSpline3& spline, const double point, const int whichKnot)
  {
    double knotX, knotY, b, c, d;
    spline.GetCoeff(whichKnot, knotX, knotY, b, c, d);
    //const double dx = point - knotX;
  
    const double knotX2 = knotX*knotX, knotX3 = knotX*knotX*knotX;
    return point*(knotY - b*knotX + c*knotX2 - d*knotX3 + point*(b*0.5 - c*knotX + 3./2.*d*knotX2 + point*(c/3. - d*knotX + point*d*0.25)));
    //dx*dx*dx = (point - knotX)*(point - knotX)*(point - knotX)
    //         = (point*point - 2*point*knotX + knotX*knotX)*(point - knotX)
    //         = (point*point*point - 2*point*point*knotX + point*knotX*knotX - point*point*knotX + 2*point*knotX*knotX - knotX*knotX*knotX)
    //         = (point*point*point - 3*point*point*knotX + 3*knotX*knotX*point - knotX*knotX*knotX)
  }
  
  double integral(TSpline3& spline, const double start, const double end)
  {
    //The knot for a given x always has a data point at fX < x.
    const auto startKnot = spline.FindX(start),
               endKnot = spline.FindX(end);
  
    double knotX, knotY;
    spline.GetKnot(startKnot+1, knotX, knotY);
  
    double integral = antiderivative(spline, std::min(knotX, end), startKnot) - antiderivative(spline, start, startKnot); //First point.  Apparently, I was always getting this right?
  
    double secondKnotX = knotX;  //I never use knotY anyway, so just reuse it to save one variable.
    for(int whichKnot = startKnot+1; whichKnot < endKnot; ++whichKnot)
    {
      knotX = secondKnotX;
      spline.GetKnot(whichKnot+1, secondKnotX, knotY);
      integral += antiderivative(spline, secondKnotX, whichKnot) - antiderivative(spline, knotX, whichKnot);
    }
  
    if(startKnot != endKnot)
    {
      spline.GetKnot(endKnot, knotX, knotY);
      integral += antiderivative(spline, end, endKnot) - antiderivative(spline, knotX, endKnot);
    }
  
    return integral;
  }
}

template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
class NeutronInelasticReweighter: public PlotUtils::Reweighter<UNIVERSE, EVENT>
{
  public:
    NeutronInelasticReweighter(const std::map<std::string, std::vector<int>>& fileNameToFS): fTotalInelastic("inelastic", {}), fKinENormalization(nullptr), fGeometry()
    {
      fChannels.reserve(fileNameToFS.size()); //If I don't use this, the program will often crash.  std::vector::emplace_back() will have to
                                              //reallocate memory many times.  When it does that, it copies the old Channels is made and then
                                              //deletes the originals.  But the copied TF1s hold lambda functions that still point at the
                                              //original (now deleted) Channels.

      //Load fKinENormalization from a file.  Do this first because it can fail.
      const std::string kinEFileName = "MoNA_FS_normalizations.root", 
                        kinENormHistName = "Tracker_Signal_FSParticleKE_Truth_Neutron";
      auto oldPwd = gDirectory;
      try
      {
        std::unique_ptr<TFile> kinEFile(TFile::Open(kinEFileName.c_str()));
        if(kinEFile)
        {
          fKinENormalization = dynamic_cast<TH1D*>(kinEFile->Get(kinENormHistName.c_str())->Clone()); //Make a Clone() so I don't have to keep kinEFile open while the job runs.
          //if(!fKinENormalization) throw std::runtime_error("Failed to load a histogram named " + kinENormHistName + " from a file named " + kinEFileName + " for neutron inelastic reweight normalization.");
          if(fKinENormalization) fKinENormalization->SetDirectory(nullptr); //Make sure fKineENormalization is no longer tied to its parent object's file because that file will eventually be closed.
        }
      }
      catch(const ROOT::exception& /*e*/)
      {
        std::cerr << "Failed to load neutron inelastic reweight's renormalization file from " << kinEFileName << ".  Proceeding without renormalization...\n";
      }

      //Load total elastic cross section from a file
      {
        std::string weightFileDir = "";
        if(std::getenv("PLOTUTILSROOT")) weightFileDir = std::string(std::getenv("PLOTUTILSROOT")) + "/data/neutronInelasticReweight/";

        std::unique_ptr<TFile> totalElasticFile(TFile::Open((weightFileDir + "minerva_neutron_cross_sections.root").c_str()));
        assert(totalElasticFile);
        auto totalElasticGraph = dynamic_cast<TGraph*>(totalElasticFile->Get("elastic"));
        assert(totalElasticGraph);
        fTotalElasticSpline = TSpline3("elastic", totalElasticGraph);
      }

      //Load interaction channels from files
      for(const auto& channel: fileNameToFS) fChannels.emplace_back(channel.first, channel.second);

      fLowestMinKE = std::min_element(fChannels.begin(), fChannels.end(), [](const auto& lhs, const auto& rhs) { return lhs.fMin < rhs.fMin; })->fMin;
      fHighestMaxKE = std::max_element(fChannels.begin(), fChannels.end(), [](const auto& lhs, const auto& rhs) { return lhs.fMax < rhs.fMax; })->fMax;

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

      double fMin;
      double fMax;

      //Associate previous cross section integrals with each channel for when Ti and Tf don't change much.
      /*mutable double fOldSigmaCache;
      mutable double fNewSigmaCache;*/

      Channel(): fInelasticChildren(), fMin(), fMax(), fOldSigmaRatioSpline(), fNewSigmaRatioSpline()
      {
      }

      Channel(const std::string& channelName, const std::vector<int> inelChildren): fInelasticChildren(inelChildren.begin(), inelChildren.end())
      {
        std::string weightFileDir = "";
        if(std::getenv("PLOTUTILSROOT")) weightFileDir = std::string(std::getenv("PLOTUTILSROOT")) + "/data/neutronInelasticReweight/";

        const std::string oldFileName = weightFileDir + "minerva_neutron_cross_sections.root";
        std::unique_ptr<TFile> oldGraphFile(TFile::Open(oldFileName.c_str()));
        if(!oldGraphFile) throw std::runtime_error("Failed to open a file named " + oldFileName + " for a GEANT cross section graph in InelasticNeutronReweighter::Channel.");
        auto oldRatioGraph = dynamic_cast<TGraph*>(oldGraphFile->Get(channelName.c_str()));
        if(!oldRatioGraph) throw std::runtime_error("Failed to load a TGraph named " + channelName + " from a file named " + oldFileName + " for NeutronInelasticReweighter::Channel");
        fOldSigmaRatioSpline = TSpline3(channelName.c_str(), oldRatioGraph);

        TGraph newRatioGraph((weightFileDir + channelName + ".csv").c_str());
        fNewSigmaRatioSpline = TSpline3(channelName.c_str(), &newRatioGraph);

        fMin = std::max(findFirstNonZeroPointX(*oldRatioGraph), findFirstNonZeroPointX(newRatioGraph));
        fMax = std::min(findLastNonZeroPointX(*oldRatioGraph), findLastNonZeroPointX(newRatioGraph));
      }

      double evalOldSpline(double* x, double* /*p*/) const { return fOldSigmaRatioSpline.Eval(x[0]); }
      double evalNewSpline(double* x, double* /*p*/) const { return fNewSigmaRatioSpline.Eval(x[0]); }

      //I think I have to keep these TSpline3 objects around because they're referenced by the TF1s :(
      //TODO: Get the ROOT authors of TSpline3 to be const-correct!
      mutable TSpline3 fOldSigmaRatioSpline;
      mutable TSpline3 fNewSigmaRatioSpline;
    };

    std::vector<Channel> fChannels; //channels that will be reweighted
                                            //and because ROOT couldn't be bothered to get const-ness right on TSpline3::Eval().
    //Channel fOther; //All other channels that aren't reweighted are lumped into one.  This keeps the total inelastic cross section the same.
    Channel fTotalInelastic; //Wouldn't be needed if I could get away with just weighting each neutron by exclusive cross section ratio
    mutable TSpline3 fTotalElasticSpline;

    //KE range which at least some channel covers.  Any neutrons outside of this range just get a weight of 1.
    double fLowestMinKE;
    double fHighestMaxKE;

    TH1D* fKinENormalization; //Normalization to keep the overall neutrino cross section the same in kinetic energy and angle

    PlotUtils::TargetUtils fGeometry;

    double getInteractingWeight(const Channel& channel, const double density, const double Ti, const double Tf) const;
    double getOtherInelasticWeight(const double density, const double Ti, const double Tf) const;
    double getConstantChannelWeight(const double density, const double Ti, const double Tf) const;
    double getNoInteractionWeight(const double density, const double Ti, const double Tf) const;

    double evalSigmaRatio(TSpline3& sigmaSpline, double Ti, double Tf, const double min, const double max) const;
};

template <class UNIVERSE, class EVENT>
double NeutronInelasticReweighter<UNIVERSE, EVENT>::GetWeight(const UNIVERSE& univ, const EVENT& /*event*/) const
{
  double weight = 1;

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

    //Possibly-elastic points where inelastic interaction did not happen
    //Stop before the last point because it may have ended with an inelastic interaction
    for(int whichPoint = startPoint; whichPoint < endPoint-1; ++whichPoint)
    {
      const int intCode = intCodePerPoint[whichPoint];
      assert(intCode == 0 || intCode == 2 || intCode == 3);

      //N.B.: material of -6 seems to be a special flag Jeffrey added to denote CH scintillator as opposed to pure carbon from target 3.
      if(materialPerPoint[whichPoint] == 6 && fGeometry.InTracker(xPerPoint[whichPoint], yPerPoint[whichPoint], zPerPoint[whichPoint]))
      {
        const double Ti = startEnergyPerPoint[whichPoint] - ::neutronMass,
                     Tf = endEnergyPerPoint[whichPoint] - ::neutronMass;
        const double density = densityPerPoint[whichPoint];

        //Developer's note: oldInel is calculated on both if branches.  Right now, this is done by independent functions.  Calculating it outside the if block would at least make a GPU happier.  Not so sure about a CPU...
        if(intCode == 3) weight *= getConstantChannelWeight(density, Ti, Tf); //elastic interacting
        else weight *= getNoInteractionWeight(density, Ti, Tf); //if no interaction.  Signaled by intCode is 0 or 2, according to MnvHadronReweight comments
      } //If point is in tracker
    } //For each point in whichNeutron's trajectory
    assert(!isinf(weight));

    //Weight for final trajectory point.  This is a special case because it's the only time inelastic interactions can happen.
    if(startPoint != endPoint && fGeometry.InTracker(xPerPoint[endPoint - 1], yPerPoint[endPoint - 1], zPerPoint[endPoint - 1]) && materialPerPoint[endPoint - 1] == -6)
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

      const double Ti = startEnergyPerPoint[endPoint - 1] - ::neutronMass,
                   Tf = endEnergyPerPoint[endPoint - 1] - ::neutronMass,
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
        else weight *= getConstantChannelWeight(density, Ti, Tf); //getOtherInelasticWeight(density, Ti, Tf); //getConstantChannelWeight(density, Ti, Tf);
      }
      else if(intCode == 3) weight *= getConstantChannelWeight(density, Ti, Tf); //Trajectory ends with an elastic interaction
      else weight *= getNoInteractionWeight(density, Ti, Tf); //If this trajectory ended by some process other than an inelastic interaction

      if(isinf(weight)) std::cout << "weight = " << weight << " is now inf at Ti = " << Ti << ", Tf = " << Tf << ", and intCode = " << intCode << ".  getNoInteractionWeight returns " << getNoInteractionWeight(density, Ti, Tf) << "\n";
    } //If last point is in the tracker and CH scintillator

    //Divide by a kinematics-dependent normalization factor to keep the total neutrino cross section constant.
    //TODO: Aaron only does this for the leading particle in the original MnvHadronReweight
    //TODO: Aaron does this using FS particle branches because he only cares about FS particles.  Do I have momentum components for all neutrons?
    //      If not, I'm tempted to try just reweighting in neutron KE first.
    //      Nope, I don't have neutron direction.  Trying neutron KE until I see that it's a problem.
    //weight /= fKinENormalization->GetBinContent(fKinENormalization->FindBin(startEnergyPerPoint[startPoint] - ::neutronMass));

    assert(!isinf(weight));
  } //For each neutron

  //Now, normalize so that the number of FS neutrons in the entire playlist doesn't change.
  //This should keep the total neutrino cross section from changing too.
  if(fKinENormalization)
  {
    const auto fsPDGs = univ.GetVecInt("mc_FSPartPDG");
    const auto fsEnergies = univ.GetVecDouble("mc_FSPartE");
    const auto fsPx = univ.GetVecDouble("mc_FSPartPx"), fsPy = univ.GetVecDouble("mc_FSPartPy"), fsPz = univ.GetVecDouble("mc_FSPartPz");

    const size_t nFSPart = fsPDGs.size();
    for(size_t whichFS = 0; whichFS < nFSPart; ++whichFS)
    {
      if(fsPDGs[whichFS] == 2112) //for each FS neutron
      {
        const double KE = fsEnergies[whichFS] - ::neutronMass;
        const double cosTheta = fsPz[whichFS]/std::sqrt(fsPx[whichFS]*fsPx[whichFS] + fsPy[whichFS]*fsPy[whichFS]);
        weight /= fKinENormalization->GetBinContent(fKinENormalization->FindBin(KE, cosTheta));
      }
    }
  }

  assert(!isinf(weight));
  assert(!isnan(weight));
  return weight;
}

template <class UNIVERSE, class EVENT>
double NeutronInelasticReweighter<UNIVERSE, EVENT>::getInteractingWeight(const Channel& channel, const double density, const double Ti, const double Tf) const
{
  if(Tf < channel.fMin || Ti > channel.fMax) return 1.; //When given KE outside the range where I have splines to compare to, don't reweight.

  //I don't need to reweight based on the total cross section because I'm implicitly keeping it the same.
  const double totalElastic = evalSigmaRatio(fTotalElasticSpline, Ti, Tf, fLowestMinKE, fHighestMaxKE);
  const double denom = 1. - exp(-1. * density * scintDensityToNucleons * (evalSigmaRatio(fTotalInelastic.fOldSigmaRatioSpline, Ti, Tf, fTotalInelastic.fMin, fTotalInelastic.fMax) + totalElastic));
  if(denom <= 0) return 0;
  const double num = 1. - exp(-1. * density * scintDensityToNucleons * (evalSigmaRatio(fTotalInelastic.fNewSigmaRatioSpline, Ti, Tf, fTotalInelastic.fMin, fTotalInelastic.fMax) + totalElastic));

  const double a = evalSigmaRatio(channel.fNewSigmaRatioSpline, Ti, Tf, channel.fMin, channel.fMax);
  const double b = evalSigmaRatio(channel.fOldSigmaRatioSpline, Ti, Tf, channel.fMin, channel.fMax);
  assert(!isinf(num / denom * a / b));
  return num / denom * a / b;
  //return a / b; //Case for when not changing the total inelastic cross section
}

//Weight for a channel that I'm not reweighting while still keeping the total inelastic cross section at the predicted value.
template <class UNIVERSE, class EVENT>
double NeutronInelasticReweighter<UNIVERSE, EVENT>::getOtherInelasticWeight(const double density, const double Ti, const double Tf) const
{
  int nChannelsActive = 0;
  double oldKnownInelastic = 0, newKnownInelastic = 0;
  for(const auto& channel: fChannels)
  {
    if(Tf >= channel.fMin && Ti <= channel.fMax)
    {
      ++nChannelsActive;
      oldKnownInelastic += evalSigmaRatio(channel.fOldSigmaRatioSpline, Ti, Tf, channel.fMin, channel.fMax);
      newKnownInelastic += evalSigmaRatio(channel.fNewSigmaRatioSpline, Ti, Tf, channel.fMin, channel.fMax);
    }
  }
  if(nChannelsActive < 2) return 0;

  const double totalElastic = evalSigmaRatio(fTotalElasticSpline, Ti, Tf, fLowestMinKE, fHighestMaxKE),
               oldTotalInelastic = evalSigmaRatio(fTotalInelastic.fOldSigmaRatioSpline, Ti, Tf, fTotalInelastic.fMin, fTotalInelastic.fMax),
               newTotalInelastic = evalSigmaRatio(fTotalInelastic.fNewSigmaRatioSpline, Ti, Tf, fTotalInelastic.fMin, fTotalInelastic.fMax);
  const double denom = 1. - exp(-1. * density * scintDensityToNucleons * (oldTotalInelastic + totalElastic));
  if(denom <= 0) return 0;
  const double num = 1. - exp(-1. * density * scintDensityToNucleons * (newTotalInelastic + totalElastic));

  double a = newTotalInelastic - newKnownInelastic;
  double b = oldTotalInelastic - oldKnownInelastic;

  //std::cout << "Other channel ratio is " << a / b << "\n";

  //TODO: Remove the following debugging lines
  //if(a < 0 && a > -1) a = 0; //Small disagreement between splines where there's just the nGamma spline.

  /*if(a < 0)
  {
    std::cout << "Got a negative new cross section at Ti = " << Ti << " and Tf = " << Tf << " for \"Other\" channel: " << a << "\nTotal new cross section is " << newTotalInelastic << "\n";
    std::cout << "Channels in this region are:\n";
    for(const auto& channel: fChannels)
    {
      std::cout << channel.fNewSigmaRatioSpline.GetName() << ": ";
      if(Tf < channel.fMin || Ti > channel.fMax) std::cout << "0\n";
      else std::cout << evalSigmaRatio(channel.fNewSigmaRatioSpline, Ti, Tf, channel.fMin, channel.fMax) << "\n";
    }
  }
  if(b < 0)
  {
    std::cout << "Got a negative old cross section for \"Other\" channel: " << b << "\nTotal old cross section is " << oldTotalInelastic << "\n";
    std::cout << "Channels in this region are:\n";
    for(const auto& channel: fChannels)
    {
      std::cout << channel.fOldSigmaRatioSpline.GetName() << ": ";
      if(Tf < channel.fMin || Ti > channel.fMax) std::cout << "0\n";
      else std::cout << evalSigmaRatio(channel.fOldSigmaRatioSpline, Ti, Tf, channel.fMin, channel.fMax) << "\n";
    }
  }*/

  assert(!isinf(num / denom * a / b));
  assert(a >= 0);
  assert(b >= 0);
  return num / denom * a / b;
}

//Weight for a channel that I'm not actually reweighting.  It turns out not to be 1 if I work out the math for MnvHadronReweight.
//I use it in multiple places, so I'm making it a function to force myself to be consistent.
template <class UNIVERSE, class EVENT>
double NeutronInelasticReweighter<UNIVERSE, EVENT>::getConstantChannelWeight(const double density, const double Ti, const double Tf) const
{
  const double totalElastic = evalSigmaRatio(fTotalElasticSpline, Ti, Tf, fLowestMinKE, fHighestMaxKE);
  const double oldTotal = evalSigmaRatio(fTotalInelastic.fOldSigmaRatioSpline, Ti, Tf, fTotalInelastic.fMin, fTotalInelastic.fMax) + totalElastic;
  const double denom = 1. - exp(-1. * density * scintDensityToNucleons * oldTotal);
  if(denom > 0) //Otherwise, don't reweight at all for this step
  {
    const double newTotal = evalSigmaRatio(fTotalInelastic.fNewSigmaRatioSpline, Ti, Tf, fTotalInelastic.fMin, fTotalInelastic.fMax) + totalElastic;
    const double num = 1. - exp(-1. * density * scintDensityToNucleons * newTotal);
    //Ratio of elastic fractions before and after change reduces to ratio of total cross sections when elastic stays the same!
    assert(!isinf(num / denom * oldTotal / newTotal));
    return num / denom * oldTotal / newTotal;
  }
  return 1; //else
}

template <class UNIVERSE, class EVENT>
double NeutronInelasticReweighter<UNIVERSE, EVENT>::getNoInteractionWeight(const double density, const double Ti, const double Tf) const
{
  const double oldInel = evalSigmaRatio(fTotalInelastic.fOldSigmaRatioSpline, Ti, Tf, fTotalInelastic.fMin, fTotalInelastic.fMax);
  const double newInel = evalSigmaRatio(fTotalInelastic.fNewSigmaRatioSpline, Ti, Tf, fTotalInelastic.fMin, fTotalInelastic.fMax);
  assert(!isinf(exp(-1.0 * density * scintDensityToNucleons * (newInel - oldInel))));
  return exp(-1.0 * density * scintDensityToNucleons * (newInel - oldInel)); //Should be total cross section difference, but elastic cancels out when it stays the same
}

//Adapt to graph evaluation pitfalls
template <class UNIVERSE, class EVENT>
double NeutronInelasticReweighter<UNIVERSE, EVENT>::evalSigmaRatio(TSpline3& sigmaSpline, double Ti, double Tf, const double min, const double max) const
{
  //std::cout << "Starting with Ti = " << Ti << " and Tf = " << Tf << std::endl;

  //Prefer rounding into the range where we have data over interpolating off the end of a spline
  //"clamp" Ti and Tf to min/max of ratioFunc
  Ti = std::min(Ti, max);
  Tf = std::min(Tf, max);

  //Some strange "linear interpolation towards 0" that Jeffrey does.  He also comments that this never happens in MnvHadronReweight because
  //the "HD neutron cross section" goes down to 1 MeV.
  //N.B.: ratioFunc wraps over a cubic spline to data
  /*if(Ti < min) Ti *= ratioFunc.Eval(Ti)/min; //TODO: If Ti < min, then evaluating the spline at Ti could return crazy results!
  if(Tf < min) Tf *= ratioFunc.Eval(Tf)/min;*/

  //TODO: If Ti, Tf are outside the domain of ratioFunc, I'd rather just return a weight of 1 for this event.  I think this is guaranteed by other functions now?
  Ti = std::max(min, Ti);
  Tf = std::max(min, Tf);

  double result = 0.;
  if(fabs(Ti - Tf) < 1e-6 || Ti - Tf < 0) //If Ti - Tf < 0, then the difference is probably pretty small anyway.
  {
    //std::cout << "For a function named " << ratioFunc.GetName() << ", Ti = " << Ti << " is close to Tf = " << Tf << ".  min = " << min << ".  Returning " << ratioFunc.Eval(Ti) << std::flush << std::endl;
    //return ratioFunc.Eval(Ti);
    result = sigmaSpline.Eval(Ti);
  }
  else result = integral(sigmaSpline, Tf, Ti)/(Ti - Tf); //ratioFunc.Integral(Ti, Tf, 1e-6)/(Tf - Ti); //TF1::Integral() is supposedly a Gaussian quadrature algorithm in some cases

  if(result < 0) std::cout << "result = " << result << " < 0!  Ti = " << Ti << ", Tf = " << Tf << " for spline " << sigmaSpline.GetTitle() << ".  Ti - Tf = " << Ti - Tf << "\n";
  assert(result >= 0);
  return result;
}
