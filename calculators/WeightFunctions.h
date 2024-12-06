#ifndef WEIGHTFUNCTIONS_H
#define WEIGHTFUNCTIONS_H

#include "TruthFunctions.h"
#include "PlotUtils/NSFDefaults.h" 

// Get Weights
virtual double GetGenieWeight() const {
  double nonResPiWgt = UseNonResPiReweight() && PlotUtils::IsNonResPi(*this)
                     ? PlotUtils::kNonResPiWeight : 1.;
  double deutWgt = UseDeuteriumGeniePiTune() && PlotUtils::IsCCRes(*this) ? 
                   ( PlotUtils::GetGenieParReweight(*this,"truth_genie_wgt_MaRES", 
                                                  NSFDefaults::DEUTERIUM_MaRES,
                                                  NSFDefaults::GENIE_MaRES, 
                                                  NSFDefaults::GENIE_MaRES_1Sig ) * 
                    NSFDefaults::DEUTERIUM_RES_NORM ) : 1.;
  double zexpWgt = UseZExpansionFaReweight()? GetZExpWeight() : 1;
  return nonResPiWgt * deutWgt *zexpWgt;
}

virtual double GetRPAWeight() const {
  const int variation = 0;  // CV
  return PlotUtils::GetRPAWeight(*this, Getq0True() / 1000 /* GeV */,
                                 Getq3True() / 1000 /* GeV */, variation,
                                 IsProcessingNX(),
                                 GetTargetZTrue(),
				 GetAnalysisNuPDG(),
                                 GetRPAMaterials());
}

virtual double GetLowRecoil2p2hWeight() const {
  const int variation = 0;  // CV
  return PlotUtils::GetLowRecoil2p2hWeight(*this, Getq0True() / 1000 /* GeV */,
                                           Getq3True() / 1000 /* GeV */,
                                           variation);
}

virtual double GetLowQ2PiWeight(std::string channel) const {
  int variation = 0;  // CV
  if( channel == "MENU1PI" )
  {
    //if( PlotUtils::IsCCCoh(*this) ) return 1.;
    if( !PlotUtils::IsCCNucleonPion(*this) ) return 1.;
    else
      return PlotUtils::weight_lowq2pi().getWeight(GetQ2True() * 1e-6 /*GeV^2*/,
                                                   channel, variation, GetInt("mc_targetNucleus"));
  }

  if (!PlotUtils::IsCCRes(*this) )
    return 1.;
  else
    return PlotUtils::weight_lowq2pi().getWeight(GetQ2True() * 1e-6 /*GeV^2*/,
                                                 channel, variation, GetInt("mc_targetNucleus"));
}

virtual double GetCoherentPiWeight(double thpi_true /*deg*/,
                                   double tpi_true /*GeV*/) const {
  if (GetInt("mc_intType") != 4) return 1.;
  assert(tpi_true > 0. && "GetCoherentPiWeight failed with tpi < 0.");
  assert(thpi_true > 0. && "GetCoherentPiWeight failed with thpi < 0.");
  return PlotUtils::weight_coherent().get_combined_weight(thpi_true, tpi_true);
}

virtual double GetFluxAndCVWeight(double Enu = -99. /*GeV*/,
                                  int nu_pdg = -99.) const {
  if (Enu == -99.) Enu = GetDouble("mc_incomingE") * 1e-3;
  // For LE, electron-neutrino fluxes aren't available, so we should force
  // nu_pdg to have absolute value 14. The +/- doesn't matter, because FRW sets
  // up both, anyways. This may change in the future if electron-neutrino LE
  // fluxes become available.
  if (nu_pdg == -99 && !IsPlaylistME(GetPlaylist()))
    nu_pdg = 14;
  else if (nu_pdg == -99)
    nu_pdg = GetInt("mc_incoming");
  return PlotUtils::flux_reweighter(GetPlaylist(), nu_pdg, UseNuEConstraint(),
                                    GetNFluxUniverses())
      .GetFluxCVWeight(Enu, nu_pdg);
}

virtual double GetTargetMassWeight() const {
  return 1.0;  
}

virtual double GetFSIWeight(int iWeight) const {
  static PlotUtils::weight_fsi weight_FSI;
  weight_FSI.UseTrackingThreshold();
  weight_FSI.calcWeights(
      GetInt("mc_incoming"), GetInt("mc_primaryLepton"), GetInt("mc_charm"),
      GetInt("mc_intType"), GetInt("mc_targetA"), GetInt("mc_targetZ"),
      GetInt("mc_resID"), GetInt("mc_er_nPart"), GetVecInt("mc_er_ID"),
      GetVecInt("mc_er_status"), GetVecInt("mc_er_FD"), GetVecInt("mc_er_LD"),
      GetVecInt("mc_er_mother"), GetVecDouble("mc_er_Px"),
      GetVecDouble("mc_er_Py"), GetVecDouble("mc_er_Pz"),
      GetVecDouble("mc_er_E"), m_entry);
  if (iWeight == 0)
    return weight_FSI.GetElasticWeight(1) * weight_FSI.GetAbsorptionWeight();
  if (iWeight == 1) return weight_FSI.GetElasticWeight(1);
  if (iWeight == 2) return weight_FSI.GetAbsorptionWeight();
  return 1;
}

virtual double GetMKWeight() const { 
  return PlotUtils::weight_mk().getWeight(m_chw, m_entry);
}

virtual double GetGeantHadronWeight() const {
  if (m_is_truth) return 1.;  // No efficiency reweighting for truth events
  //Don't want to spend time creating thing weighter if you're not using this
  if( m_use_mhrw_neutronCV_reweight )
  {
    SetupMHRWeighter();
    return PlotUtils::weight_hadron<PlotUtils::TreeWrapper*>().reweightNeutronCV( *this );
  }
  else return 1.;
}

virtual double GetZExpWeight() const {
  if (GetInt("mc_intType")!=1) return 1; // removed condition (mc_targetZ < 6) should be applied to H,D,He   
  const double q2 = GetDouble("mc_Q2")/(1000*1000); // Convert to GeV
  static PlotUtils::weightZExp zExpWeighter = PlotUtils::weightZExp("$MPARAMFILESROOT/data/Reweight/Z_Expansion_Reweight_v2126.root");
  return zExpWeighter.getWeight(q2);
}


virtual double GetUntrackedPionWeight() const{
  double weight2 = 1.0;


  // The Following are weights for MnvTunev4 with LowQ3Pion tune bug fix to apply to only single pions Aug 10 2023
  //std::vector<double> tpiweights = {0.139855,0.15291,0.391857,0.632001,0.804802,0.858462,0.580174,0.928241,1.15861,0.989334,1.10499,1.18818,0.835503,0.950168,1.28937,1.44994,1.32275,0.994931,1.14384,1.14747,1.08596,1.12892,1.08526,0.427295,0.248945,0.555842,0.460325,0.459072,0.70149,0.988857,1.18221,1.01983,1.27742};
  //std::vector<double> tpilowbins = {1., 10., 15., 20., 25., 30., 36., 42., 48., 54.,60., 66., 72., 78.,  84., 90., 100., 110., 125., 140., 155., 175., 200., 225., 250., 275., 300., 325., 350., 400., 500., 700., 1000.};
  //These weights where obtained with the whole data set using the P4 tuplas
  std::vector<double> tpiweights = {0.267183,0.218322,0.372796,0.58721,0.767524,0.880305,0.669767,0.817111,1.09273,0.995627,0.916708,1.24354,1.21146,1.12187,1.25325,1.19151,1.03823,1.23792,1.19056,1.22908,0.988201,1.03294,0.901374,0.757748,0.755932,0.638574,0.493987,0.391947,0.323265,0.452765,0.594541,0.768459,0.658024,0.873622};
  std::vector<double> tpilowbins = {0.0, 10.0 ,15.0, 20.0, 25.0, 30.0,36.0,42.0,48.0,54.0,60.0,66.0,72.0,78.0,84.0,90.0,96.0,102.0,110.0,125.0,140.0,155.0,175.0,200.0,225.0,250.0,275.0,300.0,325.0, 350.0,400.0,500.0,700.0, 1000.0};

  int nPions = 0;
  double tpi = 999999.;
  int pdgsize = GetInt("mc_nFSPart");
  for (int i = 0; i< pdgsize; i++)
  {
    int pdg = GetVecElem("mc_FSPartPDG", i);
    if (pdg != 211) continue;
    nPions++;
    double energy = GetVecElem("mc_FSPartE", i);
    double momentumx = GetVecElem("mc_FSPartPx", i);
    double momentumy = GetVecElem("mc_FSPartPy", i);
    double momentumz = GetVecElem("mc_FSPartPz", i);
    double pionmomentum = TMath::Sqrt(pow(momentumx, 2) + pow(momentumy,2)+pow(momentumz,2));
    double pionmass = TMath::Sqrt(pow(energy, 2) - pow(pionmomentum, 2));  
    double KE = energy - pionmass;
    if (tpi > KE) tpi = KE;
  }
  double Eavail = 0.0;
  for (int i = 0; i< pdgsize; i++)
  {
    int pdg = GetVecElem("mc_FSPartPDG", i);
    double energy = GetVecElem("mc_FSPartE", i); // hopefully this is in MeV
    if (abs(pdg) > 1e9) continue; //ignore nuclear fragments
    else if (abs(pdg) == 11 || abs(pdg) == 13) continue; //ignore leptons   
    else if (abs(pdg) == 211) Eavail+= energy - 139.5701; // subtracting pion mass to get Kinetic energy
    else if (pdg == 2212) Eavail += energy - 938.27201; // proton
    else if (pdg == 2112) continue; //Skip neutrons
    else if (pdg == 111) Eavail += energy; // pi0
    else if (pdg == 22) Eavail += energy; // photons
    else if (pdg >= 2000) Eavail += energy - 938.27201;
    else if (pdg <= -2000) Eavail += energy + 938.27201;
    else Eavail += energy;
  }

  double ptmu = GetPlepTrue() * sin(GetThetalepTrue());

  if ( nPions == 0) return 1.0;
  else if (Eavail > 1200 || ptmu > 1800) return 1.0;
  else if(GetInt("mc_intType") == 4) return 1.0;
  else {
    for (int i = 0; i< (int)tpilowbins.size(); i++){
      if (i < (int)tpilowbins.size() and tpi >= tpilowbins[i] and tpi < tpilowbins[i+1]){
        weight2 = tpiweights[i];
        break;
      }
      else if (tpi >= 1000.){
        weight2 = 1.0; //abs(tpiweights[29]);
        break; 
      }
    }
    //if (angle > 0.10) weight2 = 0.90*weight2; //Correcting for Forward going Pions

    //std::cout << "tpi re weight: " << weight2 << std::endl;
    return weight2;
  }
};

virtual double GetChargedPionTuneWeight(double tpi_true, double q2_true) const {
  // TODO consider putting in both lowq2 and tpi weights here
  // As-is, you must call them in your CV separately, with no associated
  // systematic except this one.
  return PlotUtils::weight_CCOnePi().get_weight(tpi_true, q2_true);
//  return 1;
}


#endif  // WEIGHTFUNCTIONS
