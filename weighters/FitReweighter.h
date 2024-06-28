//File: FitReweighter.h
//Brief: Reweight after sideband fit
//Author:Andrew Olivier & Vladyslav Syrotenko

#ifndef PLOTUTILS_FITREWEIGHTER_H
#define PLOTUTILS_FITREWEIGHTER_H

//PlotUtils includes
#include "utilities/NSFDefaults.h"
#include "universes/MnvTuneSystematics.h"

//Reweighter includes
#include "weighters/Reweighter.h"

// #include ""


namespace
{

bool passTrueSingleChargedPion(int nfspart, int* fspartPDG){
          int num_pion = 0;
          int num_charged_pion = 0;
          for(int i = 0; i < nfspart; i++){
            if(  fspartPDG[i]  == 211 || fspartPDG[i] == -211 ){
              num_charged_pion++;
              num_pion++;
            }
            if(  fspartPDG[i]  == 111 ){
              num_pion++;
            }
          }
          if(num_charged_pion==1 && num_pion==1) return true; //only found 1 pion and it was charged
          return false;
        }

bool passTrueCCQELike( int* mc_FSPartPDG, double* mc_FSPartE, int mc_nFSPart )
        {
          int genie_n_muons         = 0;
          int genie_n_mesons        = 0;
          int genie_n_heavy_baryons_plus_pi0s = 0;
          int genie_n_photons       = 0;
          int genie_n_protons       = 0; //antinu 

          for(int i = 0; i < mc_nFSPart; ++i) {
            int pdg =  mc_FSPartPDG[i];
            double energy = mc_FSPartE[i]; 
            double proton_E = 1058.272;
            //removing the 1020 MeV proton KE cut as per Minerba's Suggestion.
            if(true)proton_E=938.28;
            //The photon energy cut is hard-coded at 10 MeV at present. We're happy to make it general, if the need arises ! 
            if( abs(pdg) == 13 ) genie_n_muons++;
            else if( pdg == 22 && energy >10 ) genie_n_photons++;
            else if( abs(pdg) == 211 || abs(pdg) == 321 || abs(pdg) == 323 || pdg == 111 || pdg == 130 || pdg == 310 || pdg == 311 || pdg == 313 ) genie_n_mesons++;
            else if( pdg == 3112 || pdg == 3122 || pdg == 3212 || pdg == 3222 || pdg == 4112 || pdg == 4122 || pdg == 4212 || pdg == 4222 || pdg == 411 || pdg == 421 || pdg == 111 ) genie_n_heavy_baryons_plus_pi0s++;
            else if( pdg == 2212 && energy > proton_E ) genie_n_protons++; //antinu
          }

          //Definition of CCQE-like: 1 muon (from neutrino) and no mesons/heavy baryons in final state
          //Any number of final state nucleons (protons or neutrons) allowed
          //Photons from nuclear de-excitation are kept. These tend to be < 10 MeV. Events with photons from other sources are excluded. 
          //GENIE simulates nuclear de-excitations only for Oxygen atoms at present.  
          if(true){
            if( genie_n_muons         == 1 && 
          genie_n_mesons        == 0 && 
          genie_n_heavy_baryons_plus_pi0s == 0 && 
          genie_n_photons       == 0 ) return true;
          }
          else{
            if( genie_n_muons         == 1 && 
          genie_n_mesons        == 0 && 
          genie_n_heavy_baryons_plus_pi0s == 0 && 
          genie_n_photons       == 0 &&
          genie_n_protons        == 0 ) return true;
          }
          return false;
        }
bool passTrueSingleNeutralPion(int nfspart, int* fspartPDG){
        int num_pion = 0;
        int num_neutral_pion = 0;
        for(int i = 0; i < nfspart; i++){
          if(  fspartPDG[i]  == 211 || fspartPDG[i] == -211 ){
            num_pion++;
          }
          if(  fspartPDG[i]  == 111 ){
            num_neutral_pion++;
            num_pion++;
          }
        }
        if(num_neutral_pion==1 && num_pion==1) return true; //only found 1 pion and it was neutral
        return false;
     }

bool passTrueMultiPion(int nfspart, int* fspartPDG, int nerpart, int* erpartID, int* erpartstatus){
        //if( passTrueSingleChargedPion( nfspart, fspartPDG)|| passTrueSingleNeutralPion( nfspart, fspartPDG) ) return false; //Already ID'd as a single pion event
        int num_pion = 0;
        bool found_eta = false;
        for(int i = 0; i < nfspart; i++){
          if(  fspartPDG[i]  == 211 || fspartPDG[i] == -211 || fspartPDG[i] == 111){
            num_pion++;
          }
        }
        
        for(int i = 0; i < nerpart; i++){
          if(true){
            if( erpartstatus[i] == 14 && erpartID[i] == 221){
        found_eta = true;
        break;
            }
          }
          else{
            if( erpartstatus[i] == -14 && erpartID[i] == 221){
        found_eta = true;
        break;
            }
          }
        }

        if(num_pion>1||found_eta) return true; //only found 1 pion and it was charged
        return false;
      }



}

namespace PlotUtils
{
  template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
  class FitReweighter: public Reweighter<UNIVERSE, EVENT>
  {
    public:
      FitReweighter(std::map<std::string, std::vector<UNIVERSE*>>& univs): Reweighter<UNIVERSE, EVENT>()
      {
        std::unique_ptr<TFile> fSF(TFile::Open((std::string(getenv("PLOTUTILSROOT")) + "/../etc/scale_factors.root").c_str(), "READ")); // NEED TO SPECIFY ADDRESS

        //std::string(getenv("PLOTUTILSROOT")) + "/../etc/scale_factors.root"; //You might want to check that getenv("PLOTUTILSROOT") isn't NULL eventually.  It returns NULL if the environment variable PLOTUTILSROOT isn't set at all.  source opt/bin/setup.sh should set it for you.

        MnvH1D* weightsFromFileBkg1 = (MnvH1D*) fSF->Get("PiP")->Clone();
        MnvH1D* weightsFromFileBkg2 = (MnvH1D*) fSF->Get("PiN")->Clone();
        MnvH1D* weightsFromFileBkg3 = (MnvH1D*) fSF->Get("NPi")->Clone();

        bkg1WeightsWrapper.hist = weightsFromFileBkg1;
        bkg2WeightsWrapper.hist = weightsFromFileBkg2;
        bkg3WeightsWrapper.hist = weightsFromFileBkg3;
  

        for (auto band:univs) {
          for (size_t whichUniv = 0; whichUniv < band.second.size(); whichUniv++) {
            bkg1WeightsWrapper.AddUniverses(band.first, band.second[whichUniv], band.second.size(), whichUniv);
            bkg2WeightsWrapper.AddUniverses(band.first, band.second[whichUniv], band.second.size(), whichUniv);
            bkg3WeightsWrapper.AddUniverses(band.first, band.second[whichUniv], band.second.size(), whichUniv);
          }
        }
      }

      virtual ~FitReweighter() = default;

      double GetWeight(const UNIVERSE& univ, const EVENT& /*event*/) const override
      {
        //variation 0 is the CV
        // return PlotUtils::GetRPAWeight(univ, univ.Getq0True() / 1000 /* GeV */,
        //                          univ.Getq3True() / 1000 /* GeV */, 0,
        //                          univ.IsProcessingNX());
        //return bkg1WeightsWrapper.univHist(&univ)->GetBinContent(1);
        //here instead of just simple return we will need to port an algorithm of identifying prope sideband component, and return an appropriate number

        
        int nfspart = univ.GetInt("mc_nFSPart");
        std::vector<int> fspartPDG = univ.GetVecInt("mc_FSPartPDG");

        int nerpart = univ.GetInt("mc_er_nPart");
        std::vector<int> erpartID = univ.GetVecInt("mc_er_ID");
        std::vector<int> erpartstatus = univ.GetVecInt("mc_er_status");



        bool trueSingleChargedPion = passTrueSingleChargedPion(nfspart, fspartPDG.data());
        bool trueSingleNeutralPion = passTrueSingleNeutralPion(nfspart, fspartPDG.data());
        bool trueMultiPion = passTrueMultiPion (nfspart, fspartPDG.data(), nerpart, erpartID.data(), erpartstatus.data());


        //TrueCCQELike
        

        int mc_incoming = univ.GetInt("mc_incoming");
        int mc_current = univ.GetInt("mc_current");
        std::vector<int> mc_FSPartPDG = univ.GetVecInt("mc_FSPartPDG");
        std::vector<double> mc_FSPartE = univ.GetVecDouble("mc_FSPartE");
        int mc_nFSPart = univ.GetInt("mc_nFSPart");

        bool trueCCQELike = mc_incoming && mc_current && passTrueCCQELike(mc_FSPartPDG.data(), mc_FSPartE.data(), mc_nFSPart);
        

        if (!trueCCQELike) {

          if (trueSingleChargedPion) {
          	//std::cout << ""
            std::cout << "!!!TEST" << bkg1WeightsWrapper.univHist(&univ)->GetBinContent(1) << std::endl;
            std::cout << "!!!TEST" << univ.ShortName() << std::endl;
            return bkg1WeightsWrapper.univHist(&univ)->GetBinContent(1);
          } else if (trueSingleNeutralPion) {
            return bkg2WeightsWrapper.univHist(&univ)->GetBinContent(1);
          } else if (trueMultiPion) {
            return bkg3WeightsWrapper.univHist(&univ)->GetBinContent(1);
          } 

        } 
        return 1;
        
      }
                 



        // //MICHEL:
        // bool improved_michel = false;
        // if (univ.GetInt("improved_michel_vertex_type_sz") > 0) improved_michel = true;

        // if (improved_michel) return false;
        // else return true;

        // //BLOBS:
        // int n_blobs = 0;
        // int nblobs_max = 1;
        // for(int k=0; k < univ.GetInt("nonvtx_iso_blobs_start_position_z_in_prong_sz"); ++k){
        //   if(univ.GetVecElem("nonvtx_iso_blobs_start_position_z_in_prong", k) > 4750) n_blobs++;
        // }
        // bool  result =  true;
        // if( n_blobs > nblobs_max ) {
        //   result = false;
        // }
        // return result;

        // MICHEL & BLOB



        ////


        //SingleCharged:
        

      //Neutral
//   bool CCQENuCuts::passTrueSingleNeutralPion(CCQENuEvent* event){
//   return passTrueSingleNeutralPion(event->mc_nFSPart, event->mc_FSPartPDG);
// }

  
  //npi
//   bool CCQENuCuts::passTrueMultiPion(CCQENuEvent* event){
//   return passTrueMultiPion(event->mc_nFSPart, event->mc_FSPartPDG, event->mc_er_nPart, event->mc_er_ID, event->mc_er_status);
// }
//   bool CCQENuCuts::passTrueMultiPion(int nfspart, int* fspartPDG, int nerpart, int* erpartID, int* erpartstatus){
//   if( passTrueSingleChargedPion( nfspart, fspartPDG)|| passTrueSingleNeutralPion( nfspart, fspartPDG) ) return false; //Already ID'd as a single pion event
//   int num_pion = 0;
//   bool found_eta = false;
//   for(int i = 0; i < nfspart; i++){
//     if(  fspartPDG[i]  == 211 || fspartPDG[i] == -211 || fspartPDG[i] == 111){
//       num_pion++;
//     }
//   }
  
//   for(int i = 0; i < nerpart; i++){
//     if(neutrinoMode){
//       if( erpartstatus[i] == 14 && erpartID[i] == 221){
//   found_eta = true;
//   break;
//       }
//     }
//     else{
//       if( erpartstatus[i] == -14 && erpartID[i] == 221){
//   found_eta = true;
//   break;
//       }
//     }
//   }

//   if(num_pion>1||found_eta) return true; //only found 1 pion and it was charged
//   return false;
// }
      
      //TrueCCQELike
// bool CCQENuCuts::passTrueCCQELike( CCQENuEvent* event ){
//   if(neutrinoMode) return ( event->mc_incoming==14 && event->mc_current==1 && passTrueCCQELike( event->mc_FSPartPDG, event->mc_FSPartE, event->mc_nFSPart ) );
//   else return ( event->mc_incoming==-14 && event->mc_current==1 && passTrueCCQELike( event->mc_FSPartPDG, event->mc_FSPartE, event->mc_nFSPart ) );
// }


    




      

      std::string GetName() const override { return "SideBandFit"; }
      bool DependsReco() const override { return true; }
      private: 
        PlotUtils::HistWrapper<UNIVERSE> bkg1WeightsWrapper;
        PlotUtils::HistWrapper<UNIVERSE> bkg2WeightsWrapper;
        PlotUtils::HistWrapper<UNIVERSE> bkg3WeightsWrapper;
  };
};

#endif //PLOTUTILS_FITREWEIGHTER_H
