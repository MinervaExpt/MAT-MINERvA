#ifndef LOWRECOILPIONCUTS_H
#define LOWRECOILPIONCUTS_H

#include "PlotUtils/Cut.h"

namespace LowRecoilPion {

const double kMAX3DDIST = 2600.;
const double kMAX2DDIST = 150.;

template <class UNIVERSE, class EVENT>
class hasMichel : public PlotUtils::Cut<UNIVERSE, EVENT> {
 public:
  hasMichel() : PlotUtils::Cut<UNIVERSE, EVENT>("Event Has Michel ") {}

  // Make michels and return whether there are any satisfying basic quality
  // cuts.
  static bool hasMichelCut(const UNIVERSE& univ, EVENT& evt) {
    //
    // This function does two things:
    //
    // (1) Fill (and return by reference) the EVENT with quality Michels, where
    // quality means (a) is fitted, (b) has a vtx time diff < 0.4, AND (c) lies
    // within the tracker or ecal. Note, this is THE place where Michel objects
    // are made.
    //
    // (2) Returns whether there are any such quality michels.
    //
    // Subsequent cut classes/functions will remove michels from this EVENT based
    // on further quality criteria and again return whether any such quality
    // michels remain.
    //
    // TODO This function sets some stuff in the MichelEvent, (e.g. muon PT, q3,
    // eavail) which does NOT need to be in here. These are event-wide values
    // that can be accessed from the univ object.

    // if (univ.ShortName() != "cv") return true; //evt.m_allmichels.clear();
    // int nclusters = univ.GetNClusters(); //univ.GetNonMuonClusters();
    // //univ.GetNClusters(); if (nclusters > 500.) return false;
    int nmichels = univ.GetNMichels();
    if (nmichels < 1) return bool(false);
    int nfittedmich = univ.GetFittedMichelsOnly();
    if (nfittedmich < 1) return false;
    int nclusters = univ.GetNClusters();
    if (univ.ShortName() != "cv") {
      if (!evt.m_allmichels.empty()) return bool(true);
      // else return bool(true);
    }
    // if (nclusters > 500.) return false;

    // std::cout << "Number of Clusters in Event are " << nclusters <<
    // std::endl;
    for (int i = 0; i < nmichels; ++i) {
      Michel current_michel(univ, i);

      // if (univ.ShortName() == "cv") Michel current_michel(univ, i);
      // else {
      //	if (!evt.m_allmichels.empty()) auto current_michel =
      //evt.m_allmichels[i];}

      // Michel* current_michel= new Michel(univ, i);
      if (current_michel.is_fitted != 1) continue;
      if (abs(current_michel.vtx_michel_timediff) < 0.400)
        continue;  // < 0.400 is to reject dead michels that happen during dead
                   // time. >-0.250 is to see what matches look like for michels
                   // that happen before neutrino event.

      double z1 = current_michel.m_z1;
      double z2 = current_michel.m_z2;
      // Michel is in Tracker and ECAL
      if (z1 < 5980. || z2 < 5980.)
        continue;
      else if (z1 > 9038. || z2 > 9038.)
        continue;
      evt.m_allmichels.push_back(current_michel);
    }
    // std::cout << "Has Michel Cut Done looping over Michels " << std::endl;
    // Filling Event Level Info needed for 2D selection
    // double lowtpiinevent = univ.GetTrueTpi();
    // evt.lowTpi = lowtpiinevent;
    evt.nclusters = nclusters;
    evt.pT_reco = univ.GetMuonPT();
    evt.q3_reco = univ.Getq3_lowrecoil();
    evt.eavail_reco = univ.NewEavail();
    if (evt.m_allmichels.empty())
      return bool(false);
    else
      return bool(true);
    // return !evt.m_allmichels.empty();
  }

 private:
  using Michel = typename std::remove_reference<decltype(
      std::declval<EVENT>().m_nmichels.front())>::type;

  bool checkCut(const UNIVERSE& univ, EVENT& evt) const {
    return hasMichelCut(univ, evt);
  }
};

template <class UNIVERSE, class EVENT>
class BestMichelDistance2D : public PlotUtils::Cut<UNIVERSE, EVENT> {
 public:
  BestMichelDistance2D(const double maxDistance = kMAX2DDIST)
      : PlotUtils::Cut<UNIVERSE, EVENT>(
            "Per Michel 2D Distance in at least two views is < " +
            std::to_string(maxDistance) + "mm"),
        m_maxDistance(maxDistance) {}

  // Update our michels and return whether there any satisfying 2D distance
  // cut.
  static bool BestMichelDistance2DCut(const UNIVERSE& univ, EVENT& evt,
                                  const double max_dist = kMAX2DDIST) {
      // This function will get an integer for the best match type of the Michel.
      // It compares distance between Michel and whatever its best match is to find
      // the Best type of Michel for a single Michel.

      //std::cout << "Implementing Michel Cut with 2D distance of " << max_dist << " mm" << std::endl;
      std::vector<Michel> nmichelspass;
      int nmichels = evt.m_allmichels.size();
      if (nmichels == 0) return bool(false);
      for (unsigned int i = 0; i < evt.m_allmichels.size(); i++)
      {
        int upvtxmatch = 0;
        int downvtxmatch = 0;
        int upclusmatch = 0;
        int downclusmatch = 0;



        // For Vertex Match Check to see if 2D distance cut will
        double upvtxXZ = evt.m_allmichels[i].up_to_vertex_XZ;
        double downvtxXZ = evt.m_allmichels[i].down_to_vertex_XZ;
        double upvtxUZ = evt.m_allmichels[i].up_to_vertex_UZ;
        double downvtxUZ = evt.m_allmichels[i].down_to_vertex_UZ;
        double upvtxVZ = evt.m_allmichels[i].up_to_vertex_VZ;
        double downvtxVZ = evt.m_allmichels[i].down_to_vertex_VZ;


        std::vector<double> upvtx = {upvtxXZ, upvtxUZ, upvtxVZ};
        std::vector<double> downvtx = {upvtxXZ, upvtxUZ, upvtxVZ};

        std::sort(upvtx.begin(), upvtx.end());
        std::sort(downvtx.begin(), downvtx.end());

        //if(upvtx[0] < max_dist && upvtx[1] < max_dist) evt.m_allmichels[i].passable_matchtype.at(0) = 1;
        //if(downvtx[0] < max_dist && downvtx[1] < max_dist) evt.m_allmichels[i].passable_matchtype.at(1) = 2;       


        if (upvtxXZ < max_dist && (upvtxUZ < max_dist || upvtxVZ < max_dist)) evt.m_allmichels[i].passable_matchtype.at(0) = 1;
        else if (upvtxUZ < max_dist && (upvtxXZ < max_dist || upvtxVZ < max_dist)) evt.m_allmichels[i].passable_matchtype.at(0) = 1;
        else if (upvtxVZ < max_dist && (upvtxXZ < max_dist || upvtxUZ < max_dist)) evt.m_allmichels[i].passable_matchtype.at(0) = 1;
        else evt.m_allmichels[i].passable_matchtype.at(0) = -1;

        if (downvtxXZ < max_dist && (downvtxUZ < max_dist || downvtxVZ < max_dist)) evt.m_allmichels[i].passable_matchtype.at(1) = 2;
        else if (downvtxUZ < max_dist && (downvtxXZ < max_dist || downvtxVZ < max_dist)) evt.m_allmichels[i].passable_matchtype.at(1) = 2;
        else if (downvtxVZ < max_dist && (downvtxXZ < max_dist || downvtxUZ < max_dist)) evt.m_allmichels[i].passable_matchtype.at(1) = 2;
        else evt.m_allmichels[i].passable_matchtype.at(1) = -1;

        double upclusXZ = evt.m_allmichels[i].up_to_clus_XZ;
        double upclusUZ = evt.m_allmichels[i].up_to_clus_UZ;
        double upclusVZ = evt.m_allmichels[i].up_to_clus_VZ;
        double downclusXZ = evt.m_allmichels[i].down_to_clus_XZ;
        double downclusUZ = evt.m_allmichels[i].down_to_clus_UZ;
        double downclusVZ = evt.m_allmichels[i].down_to_clus_VZ;

        std::vector<double> upclus = {upclusXZ, upclusUZ, upclusVZ};
        std::vector<double> downclus = {downclusXZ, downclusUZ,downclusVZ};

        std::sort(upclus.begin(),upclus.end());
        std::sort(downclus.begin(),downclus.end());

        //if(upclus[0] < max_dist && upclus[1] < max_dist) evt.m_allmichels[i].passable_matchtype.at(2) = 3;
        //if(downclus[0] < max_dist && downclus[1] < max_dist) evt.m_allmichels[i].passable_matchtype.at(3) = 4;    


        if ( upclusXZ < max_dist && (upclusUZ < max_dist || upclusVZ < max_dist)) evt.m_allmichels[i].passable_matchtype.at(2) = 3;
        else if (upclusUZ < max_dist && (upclusXZ < max_dist || upclusVZ < max_dist)) evt.m_allmichels[i].passable_matchtype.at(2) = 3;
        else if (upclusVZ < max_dist && (upclusXZ < max_dist || upclusUZ < max_dist)) evt.m_allmichels[i].passable_matchtype.at(2) = 3;
        else evt.m_allmichels[i].passable_matchtype.at(2) = -1;

        if (downclusXZ < max_dist && (downclusUZ < max_dist || downclusVZ < max_dist)) evt.m_allmichels[i].passable_matchtype.at(3) = 4;
        else if (downclusUZ < max_dist && (downclusXZ < max_dist || downclusVZ < max_dist)) evt.m_allmichels[i].passable_matchtype.at(3) = 4;
        else if (downclusVZ < max_dist && (downclusXZ < max_dist || downclusUZ < max_dist)) evt.m_allmichels[i].passable_matchtype.at(3) = 4;
        else evt.m_allmichels[i].passable_matchtype.at(3) = -1;

        if (evt.m_allmichels[i].passable_matchtype.at(1) == -1 && evt.m_allmichels[i].passable_matchtype.at(2) == -1 && evt.m_allmichels[i].passable_matchtype.at(3) == -1 && evt.m_allmichels[i].passable_matchtype.at(0) == -1) continue;

        std::vector<double> distances3D;

        if (evt.m_allmichels[i].passable_matchtype[0] == 1) distances3D.push_back(evt.m_allmichels[i].up_to_vertex_dist3D); //Distance between michel to vertex
        if (evt.m_allmichels[i].passable_matchtype[1] == 2) distances3D.push_back(evt.m_allmichels[i].down_to_vertex_dist3D); //distancebetween michel to vertex 
        if (evt.m_allmichels[i].passable_matchtype[2] == 3) distances3D.push_back(evt.m_allmichels[i].up_clus_michel_dist3D); //distnace between michel to cluster
        if (evt.m_allmichels[i].passable_matchtype[3] == 4) distances3D.push_back(evt.m_allmichels[i].down_clus_michel_dist3D); // distance between michel to cluster
        if (distances3D.empty()) distances3D = {9999.,9999.,9999.,9999.};
        //if (distances3D[0] == 9999.) continue; // Comment this line out if you are making efficiency plots
        //std::cout << "Passable Match Types for this Michel are " << evt.m_allmichels[i].passable_matchtype.at(0) << evt.m_allmichels[i].passable_matchtype.at(1) << evt.m_allmichels[i].passable_matchtype.at(2) << evt.m_allmichels[i].passable_matchtype.at(3) << std::endl;

        std::sort(distances3D.begin(), distances3D.end());

        if (distances3D[0] == 9999.) continue;
        else if (distances3D[0] == evt.m_allmichels[i].up_to_vertex_dist3D){
          evt.m_allmichels[i].best_XZ = evt.m_allmichels[i].up_to_vertex_XZ;
          evt.m_allmichels[i].best_UZ = evt.m_allmichels[i].up_to_vertex_UZ;
          evt.m_allmichels[i].best_VZ = evt.m_allmichels[i].up_to_vertex_VZ;

          evt.m_allmichels[i].BestMatch = 1;
          evt.m_allmichels[i].Best3Ddist = evt.m_allmichels[i].up_to_vertex_dist3D;

          //std::cout << "This  Michel is UPVTX and has true endpoint " << evt.m_allmichels[i].trueEndpoint << std::endl;
        }
        else if (distances3D[0] == evt.m_allmichels[i].down_to_vertex_dist3D){
          evt.m_allmichels[i].BestMatch = 2;
          evt.m_allmichels[i].best_XZ = evt.m_allmichels[i].down_to_vertex_XZ;
          evt.m_allmichels[i].best_UZ = evt.m_allmichels[i].down_to_vertex_UZ;
          evt.m_allmichels[i].best_VZ = evt.m_allmichels[i].down_to_vertex_VZ;
          evt.m_allmichels[i].Best3Ddist = evt.m_allmichels[i].down_to_vertex_dist3D;
          //std::cout << "This  Michel is DOWNVTX and has true endpoint " << evt.m_allmichels[i].trueEndpoint << std::endl;

        }
        else if (distances3D[0] == evt.m_allmichels[i].up_clus_michel_dist3D){
          evt.m_allmichels[i].BestMatch = 3;
          evt.m_allmichels[i].best_XZ = evt.m_allmichels[i].up_to_clus_XZ;
          evt.m_allmichels[i].best_UZ = evt.m_allmichels[i].up_to_clus_VZ;
          evt.m_allmichels[i].best_VZ = evt.m_allmichels[i].up_to_clus_UZ;
          evt.m_allmichels[i].Best3Ddist = evt.m_allmichels[i].up_clus_michvtx_dist3D;
          //std::cout << "This  Michel is UPCLUS and has true endpoint " << evt.m_allmichels[i].trueEndpoint << std::endl;

        }
        else if (distances3D[0] == evt.m_allmichels[i].down_clus_michel_dist3D){
          evt.m_allmichels[i].BestMatch = 4;
          evt.m_allmichels[i].Best3Ddist = evt.m_allmichels[i].down_clus_michvtx_dist3D;
          evt.m_allmichels[i].best_XZ = evt.m_allmichels[i].down_to_clus_XZ;
          evt.m_allmichels[i].best_UZ = evt.m_allmichels[i].down_to_clus_UZ;
          evt.m_allmichels[i].best_VZ = evt.m_allmichels[i].down_to_clus_VZ;
          //std::cout << "This  Michel is DOWNCLUS and has true endpoint " << evt.m_allmichels[i].trueEndpoint << std::endl;

        }
        else{
          evt.m_allmichels[i].BestMatch = 0;
          evt.m_allmichels[i].Best3Ddist = 9999.;
          evt.m_allmichels[i].best_XZ = 9999.;
          evt.m_allmichels[i].best_UZ = 9999.;
          evt.m_allmichels[i].best_VZ = 9999.;

          continue;
        }


        if (distances3D[1] == evt.m_allmichels[i].up_to_vertex_dist3D) evt.m_allmichels[i].SecondBestMatch = 1;
        else if (distances3D[1] == evt.m_allmichels[i].down_to_vertex_dist3D) evt.m_allmichels[i].SecondBestMatch = 2;
        else if (distances3D[1] == evt.m_allmichels[i].up_clus_michel_dist3D) evt.m_allmichels[i].SecondBestMatch = 3;
        else if (distances3D[1] == evt.m_allmichels[i].down_clus_michel_dist3D) evt.m_allmichels[i].SecondBestMatch = 4;
        else {evt.m_allmichels[i].SecondBestMatch = 0;} 

        if (evt.m_allmichels[i].best_XZ < evt.m_allmichels[i].best_UZ and evt.m_allmichels[i].best_XZ < evt.m_allmichels[i].best_VZ) evt.m_allmichels[i].bestview = 1;

        else if (evt.m_allmichels[i].best_UZ < evt.m_allmichels[i].best_XZ and evt.m_allmichels[i].best_UZ < evt.m_allmichels[i].best_VZ) evt.m_allmichels[i].bestview = 2; 

        else if (evt.m_allmichels[i].best_VZ < evt.m_allmichels[i].best_UZ and evt.m_allmichels[i].best_VZ < evt.m_allmichels[i].best_XZ) evt.m_allmichels[i].bestview = 3; 
        int matchtype = evt.m_allmichels[i].BestMatch;
        if (matchtype == 1 || matchtype == 3) evt.m_allmichels[i].recoEndpoint = 1;
        else if (matchtype == 2 || matchtype == 4) evt.m_allmichels[i].recoEndpoint = 2;
        //evt.m_allmichels[i].GetPionAngle(univ);
        nmichelspass.push_back(evt.m_allmichels[i]);  
        evt.m_nmichelspass.push_back(evt.m_allmichels[i]);  
      }

      evt.m_nmichels.clear();      
      //std::cout << "BEst 2D Distance Cut done looping over Michels " << std::endl;
      //evt.m_allmichels.clear(); //empty existing vector of Michels
      if (!nmichelspass.empty()){
        evt.selection = 1;
        //std::cout << "Event has a Selection Michel." << std::endl; 
        evt.m_nmichels.clear();
        evt.m_nmichels = nmichelspass; // replace vector of michels with the vector of michels that passed the above cut
        return bool(true);
      }
      else{ return bool(false);}
      //evt.m_nmichelspass = nmichelspass;
      //return true; 
      //return !nmichelspass.empty();
  }

 private:
  using Michel = typename std::remove_reference<decltype(
      std::declval<EVENT>().m_nmichels.front())>::type;

  double m_maxDistance;  // Maximum distance from the vertex that the best
                         // Michel can have in mm

  bool checkCut(const UNIVERSE& univ, EVENT& evt) const {
    return BestMichelDistance2DCut(univ, evt, m_maxDistance);
  }
};

template <class UNIVERSE, class EVENT>
class GetClosestMichel : public PlotUtils::Cut<UNIVERSE, EVENT> {
 public:
  GetClosestMichel(const int michelgroup, const double max_dist = kMAX3DDIST)
      : PlotUtils::Cut<UNIVERSE, EVENT>(
            "Getting Closest Michel for Michel Group " +
            std::to_string(michelgroup)),
        m_max_dist(max_dist) {}
  // michel group means getting Closest Michel in selection (evt.m_nmichelspass)
  // or sideband group (m_sidebandpass)

  // Update our michels and return whether there any satisfying 3D distance
  // cut.
  static bool GetClosestMichelCut(const UNIVERSE& univ, EVENT& evt,
                                  const double max_dist = kMAX3DDIST) {
    // evt.m_nmichels.clear();
    // if (michelgroup == 0) evt.m_nmichels = evt.m_nmichelspass;
    // else if (michelgroup == 1) evt.m_nmichels = evt.m_sidebandpass;
    int noverlay = 0.0;
    int nmichels = evt.m_nmichels.size();
    std::vector<Michel> closestMichel;
    if (nmichels == 0) return false;
    evt.m_bestdist = 9999.;  // setting some default value for best distance
    std::vector<double> allmichel3Ddist;
    for (int i = 0; i < nmichels; ++i) {
      double dist =
          evt.m_nmichels[i].Best3Ddist;  // getting the minimum pion range
                                         // (vertex to Michel/Clus distance)
      allmichel3Ddist.push_back(dist);

    }  // Get all the distances for the michels that pass the 2D dist cut

    std::sort(allmichel3Ddist.begin(), allmichel3Ddist.end());

    for (int i = 0; i < nmichels; ++i) {
      double dist = evt.m_nmichels[i].Best3Ddist;
      int order = i + 1;
      if (dist == allmichel3Ddist[i]) evt.m_nmichels[i].OrderOfMichel = order;

      // std::cout << "Printing Order of Michels " << order << " , dist: "<<
      // dist << std::endl;
    }
    // std::cout << "Printing Closest Distance in event : " <<
    // allmichel3Ddist[0] << std::endl;
    for (int i = 0; i < nmichels; ++i) {
      evt.m_nmichels[i].GetPionAngle(univ);
      double dist = evt.m_nmichels[i].Best3Ddist;
      evt.m_nmichelspass.clear();
      if (dist < 2600.) {
        // std::cout << evt.m_nmichelspass.size() << std::endl;
        evt.m_nmichelspass.push_back(evt.m_nmichels[i]);
      }
      if (evt.m_nmichels[i].OrderOfMichel == 1) {
        // closestMichel.push_back(evt.m_nmichels[i]);
        // if (dist > 1000.) return false; // mimicing aaron's cuts to remove
        // high tpi events > 350 MeV. 1000 - 1200 mm in range is approximately
        // the bin with the most events in that tpi bins. .
        evt.m_bestdist = dist;
        evt.m_idx = i;
        if (evt.m_nmichels[i].overlay_fraction > 0.5)
          evt.ClosestMichelsIsOverlay = 1;
        evt.m_best_XZ = evt.m_nmichels[i].best_XZ;
        evt.m_best_UZ = evt.m_nmichels[i].best_UZ;
        evt.m_best_VZ = evt.m_nmichels[i].best_VZ;
        evt.m_matchtype = evt.m_nmichels[i].BestMatch;
        int bmatch = evt.m_nmichels[i].BestMatch;
        if (bmatch == 1 || bmatch == 3) {
          evt.best_x = evt.m_nmichels[i].m_x1;
          evt.best_y = evt.m_nmichels[i].m_y1;
          evt.best_z = evt.m_nmichels[i].m_z1;
        } else if (bmatch == 2 || bmatch == 4) {
          evt.best_x = evt.m_nmichels[i].m_x2;
          evt.best_y = evt.m_nmichels[i].m_y2;
          evt.best_z = evt.m_nmichels[i].m_z2;
        }
        evt.b_truex = evt.m_nmichels[i].true_initialx;
        evt.b_truey = evt.m_nmichels[i].true_initialy;
        evt.b_truez = evt.m_nmichels[i].true_initialz;
        closestMichel.push_back(evt.m_nmichels[i]);
      } else {
        // if (evt.m_nmichels[i].OrderOfMichel == 2)
        // closestMichel.push_back(evt.m_nmichels[i]);
        evt.m_idx = -1;
      }
    }
    if (closestMichel.empty()) return false;
//    closestMichel[0].GetPionAngle(univ);
    double lowtpiinevent = closestMichel[0].pionKE;
    // if (closestMichel[0].Best3Ddist > 1200.) continue; //Mimicking Aaron's
    // cuts to remove high Tpi events > 350 MeV.

    // std::cout << "Closest Michel Pion KE is " << lowtpiinevent << std::endl;
    // if (closestMichel[0].BestMatch == 1 || closestMichel[0].BestMatch == 2)
    // univ.PrintTrueArachneLink();
    evt.lowTpi = lowtpiinevent;
    evt.pT_reco = univ.GetMuonPT();
    evt.q3_reco = univ.Getq3();
    evt.eavail_reco = univ.NewEavail();
    evt.m_nmichels.clear();
    evt.m_nmichels = closestMichel;
    evt.m_bestthetaangle = evt.m_nmichels[evt.m_idx].best_angle; 
    // std::cout << "Get Closest Distance Michel Done looping over Michels.
    // Closest Distance is: " << closestMichel[0].Best3Ddist << std::endl;
    if (univ.GetMuonPT() < .20 and univ.NewEavail() < 50. and
        !evt.m_nmichels.empty()) {
      double vtx_t = univ.GetVertex().T() / pow(10, 3);  // mus

      // univ.PrintTrueArachneLink();
      //  univ.PrintDataArachneLink();
      //  std::cout << "Available Energy: " << univ.NewEavail() << " Muon pT
      //  Reco: " << univ.GetMuonPT()  << " Primary Vtx time: " << vtx_t << "
      //  Best Michel at time: " << evt.m_nmichels[0].time << " range: " <<
      //  evt.m_nmichels[0].Best3Ddist << " energy: " <<
      //  evt.m_nmichels[0].energy << " Matched to : " <<
      //  evt.m_nmichels[0].BestMatch << std::endl;
    }
    // if (closestMichel[0].Best3Ddist > 60.) return false;
    // if (closestMichel[0].Best3Ddist < 60 || closestMichel[0].Best3Ddist >
    // 1338.) return false; if (evt.m_nmichelspass.size() > 1) return false; //
    // Lets look at single pi+ events only
    if (closestMichel[0].Best3Ddist > 2600.)
      return false;
    else
      return true;
    // return !evt.m_nmichels.empty();
  };

 private:
  using Michel = typename std::remove_reference<decltype(
      std::declval<EVENT>().m_nmichels.front())>::type;
  int michelgroup;
  double m_max_dist;
  bool checkCut(const UNIVERSE& univ, EVENT& evt) const {
    return GetClosestMichelCut(univ, evt, m_max_dist);
  };
};

// Construct the michel EVENT object with all passing Michels
template <class UNIVERSE, class EVENT>
EVENT GetPassingMichels(const UNIVERSE& univ) {
  // Throw away all the return bools. We only want the evt, which is being
  // updated by reference.
  EVENT evt;
  hasMichel<UNIVERSE, EVENT>::hasMichelCut(univ, evt);
  BestMichelDistance2D<UNIVERSE, EVENT>::BestMichelDistance2DCut(univ, evt);
  GetClosestMichel<UNIVERSE, EVENT>::GetClosestMichelCut(univ, evt);
  return evt;
}

}  // namespace LowRecoilPion

#endif  // LOWRECOILPIONCUTS_H
