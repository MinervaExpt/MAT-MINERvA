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

  // Get Quality Michels -- first pass
  //
  // Quality = (1) is fitted, (2) has a vtx time diff < 0.4, AND (3) lies
  // within the tracker or ecal.
  //
  // This is the function that first makes Michels and populates the
  // MichelEvent with the good ones.
  static EVENT GetQualityMichels(const UNIVERSE& univ, const EVENT& _evt = EVENT()) {
    EVENT evt = _evt;
    for (int i = 0; i < univ.GetNMichels(); ++i) {
      Michel current_michel(univ, i);
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
    evt.nclusters = univ.GetNClusters();
    evt.pT_reco = univ.GetMuonPT();
    evt.q3_reco = univ.Getq3_lowrecoil();
    evt.eavail_reco = univ.NewEavail();
    return evt;
  }

 private:
  using Michel = typename std::remove_reference<decltype(
      std::declval<EVENT>().m_nmichels.front())>::type;

  bool checkCut(const UNIVERSE& univ, EVENT& evt) const {
    if (univ.GetNMichels() < 1) return false;
    if (univ.GetFittedMichelsOnly() < 1) return false;

    // if this function has already been called, return
    // TODO: Should this be IsVertical? I suspect so.
    if (univ.ShortName() != "cv" && !evt.m_allmichels.empty())
      return true;

    // Update our MichelEvent. In most cases, we'll be filling it for the first
    // time. Fills evt.m_allmichels, among other things.
    evt = GetQualityMichels(univ, evt);

    return !evt.m_allmichels.empty();
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

  
  // This function will get an integer for the best match type of the Michel.
  // It compares distance between Michel and whatever it's best match is to find
  // the Best type of Michel for a single Michel.
  //
  // Also update our MichelEvent evt.
  static bool BestMichelDistance2DCut(const UNIVERSE& univ, EVENT& evt,
                                      const double max_dist = kMAX2DDIST) {
    if (evt.m_allmichels.size() == 0)
     return false;

    // Update our MichelEvent evt AND fill a vector of passing michels
    std::vector<Michel> nmichelspass;
    for (unsigned int i = 0; i < evt.m_allmichels.size(); i++) {
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

      // TODO this can all just be a big OR.
      if (upvtxXZ < max_dist && (upvtxUZ < max_dist || upvtxVZ < max_dist))
        evt.m_allmichels[i].passable_matchtype.at(0) = 1;
      else if (upvtxUZ < max_dist && (upvtxXZ < max_dist || upvtxVZ < max_dist))
        evt.m_allmichels[i].passable_matchtype.at(0) = 1;
      else if (upvtxVZ < max_dist && (upvtxXZ < max_dist || upvtxUZ < max_dist))
        evt.m_allmichels[i].passable_matchtype.at(0) = 1;
      else
        evt.m_allmichels[i].passable_matchtype.at(0) = -1;

      // TODO this can all just be a big OR.
      if (downvtxXZ < max_dist &&
          (downvtxUZ < max_dist || downvtxVZ < max_dist))
        evt.m_allmichels[i].passable_matchtype.at(1) = 2;
      else if (downvtxUZ < max_dist &&
               (downvtxXZ < max_dist || downvtxVZ < max_dist))
        evt.m_allmichels[i].passable_matchtype.at(1) = 2;
      else if (downvtxVZ < max_dist &&
               (downvtxXZ < max_dist || downvtxUZ < max_dist))
        evt.m_allmichels[i].passable_matchtype.at(1) = 2;
      else
        evt.m_allmichels[i].passable_matchtype.at(1) = -1;

      double upclusXZ = evt.m_allmichels[i].up_to_clus_XZ;
      double upclusUZ = evt.m_allmichels[i].up_to_clus_UZ;
      double upclusVZ = evt.m_allmichels[i].up_to_clus_VZ;
      double downclusXZ = evt.m_allmichels[i].down_to_clus_XZ;
      double downclusUZ = evt.m_allmichels[i].down_to_clus_UZ;
      double downclusVZ = evt.m_allmichels[i].down_to_clus_VZ;

      std::vector<double> upclus = {upclusXZ, upclusUZ, upclusVZ};
      std::vector<double> downclus = {downclusXZ, downclusUZ, downclusVZ};

      std::sort(upclus.begin(), upclus.end());
      std::sort(downclus.begin(), downclus.end());

      // TODO this can all just be a big OR.
      if (upclusXZ < max_dist && (upclusUZ < max_dist || upclusVZ < max_dist))
        evt.m_allmichels[i].passable_matchtype.at(2) = 3;
      else if (upclusUZ < max_dist &&
               (upclusXZ < max_dist || upclusVZ < max_dist))
        evt.m_allmichels[i].passable_matchtype.at(2) = 3;
      else if (upclusVZ < max_dist &&
               (upclusXZ < max_dist || upclusUZ < max_dist))
        evt.m_allmichels[i].passable_matchtype.at(2) = 3;
      else
        evt.m_allmichels[i].passable_matchtype.at(2) = -1;

      // TODO this can all just be a big OR.
      if (downclusXZ < max_dist &&
          (downclusUZ < max_dist || downclusVZ < max_dist))
        evt.m_allmichels[i].passable_matchtype.at(3) = 4;
      else if (downclusUZ < max_dist &&
               (downclusXZ < max_dist || downclusVZ < max_dist))
        evt.m_allmichels[i].passable_matchtype.at(3) = 4;
      else if (downclusVZ < max_dist &&
               (downclusXZ < max_dist || downclusUZ < max_dist))
        evt.m_allmichels[i].passable_matchtype.at(3) = 4;
      else
        evt.m_allmichels[i].passable_matchtype.at(3) = -1;

      // If we don't have acceptable up and down vtx and cluster distances for
      // this Michel, it fails.
      // TODO condense? all? loop?
      if (evt.m_allmichels[i].passable_matchtype.at(1) == -1 &&
          evt.m_allmichels[i].passable_matchtype.at(2) == -1 &&
          evt.m_allmichels[i].passable_matchtype.at(3) == -1 &&
          evt.m_allmichels[i].passable_matchtype.at(0) == -1)
        continue;

      // Record 3D distance
      std::vector<double> distances3D;

      // TODO This could be condensed? For loop? case-switch?
      if (evt.m_allmichels[i].passable_matchtype[0] == 1)
        distances3D.push_back(
            evt.m_allmichels[i]
                .up_to_vertex_dist3D);  // Distance between michel to vertex
      if (evt.m_allmichels[i].passable_matchtype[1] == 2)
        distances3D.push_back(
            evt.m_allmichels[i]
                .down_to_vertex_dist3D);  // distancebetween michel to vertex
      if (evt.m_allmichels[i].passable_matchtype[2] == 3)
        distances3D.push_back(
            evt.m_allmichels[i]
                .up_clus_michel_dist3D);  // distnace between michel to cluster
      if (evt.m_allmichels[i].passable_matchtype[3] == 4)
        distances3D.push_back(
            evt.m_allmichels[i].down_clus_michel_dist3D);  // distance between
                                                           // michel to cluster
      if (distances3D.empty()) distances3D = {9999., 9999., 9999., 9999.};

      // Sort these 3D distance in ascending order.
      std::sort(distances3D.begin(), distances3D.end());

      // must have at least one good 3D distance between michel and cluster and
      // michel and vertex. Else this michel fails.
      if (distances3D[0] == 9999.)
        continue;

      // Set some MichelEvent info about this best match.
      // There's a continue here, but I think it's just an assert/sanity check.
      // TODO: condense? case-switch?
      if (distances3D[0] == evt.m_allmichels[i].up_to_vertex_dist3D) {
        evt.m_allmichels[i].best_XZ = evt.m_allmichels[i].up_to_vertex_XZ;
        evt.m_allmichels[i].best_UZ = evt.m_allmichels[i].up_to_vertex_UZ;
        evt.m_allmichels[i].best_VZ = evt.m_allmichels[i].up_to_vertex_VZ;
        evt.m_allmichels[i].BestMatch = 1;
        evt.m_allmichels[i].Best3Ddist =
            evt.m_allmichels[i].up_to_vertex_dist3D;

        // std::cout << "This  Michel is UPVTX and has true endpoint " <<
        // evt.m_allmichels[i].trueEndpoint << std::endl;
      } else if (distances3D[0] == evt.m_allmichels[i].down_to_vertex_dist3D) {
        evt.m_allmichels[i].BestMatch = 2;
        evt.m_allmichels[i].best_XZ = evt.m_allmichels[i].down_to_vertex_XZ;
        evt.m_allmichels[i].best_UZ = evt.m_allmichels[i].down_to_vertex_UZ;
        evt.m_allmichels[i].best_VZ = evt.m_allmichels[i].down_to_vertex_VZ;
        evt.m_allmichels[i].Best3Ddist =
            evt.m_allmichels[i].down_to_vertex_dist3D;
        // std::cout << "This  Michel is DOWNVTX and has true endpoint " <<
        // evt.m_allmichels[i].trueEndpoint << std::endl;
      } else if (distances3D[0] == evt.m_allmichels[i].up_clus_michel_dist3D) {
        evt.m_allmichels[i].BestMatch = 3;
        evt.m_allmichels[i].best_XZ = evt.m_allmichels[i].up_to_clus_XZ;
        evt.m_allmichels[i].best_UZ = evt.m_allmichels[i].up_to_clus_VZ;
        evt.m_allmichels[i].best_VZ = evt.m_allmichels[i].up_to_clus_UZ;
        evt.m_allmichels[i].Best3Ddist =
            evt.m_allmichels[i].up_clus_michvtx_dist3D;
        // std::cout << "This  Michel is UPCLUS and has true endpoint " <<
        // evt.m_allmichels[i].trueEndpoint << std::endl;
      } else if (distances3D[0] ==
                 evt.m_allmichels[i].down_clus_michel_dist3D) {
        evt.m_allmichels[i].BestMatch = 4;
        evt.m_allmichels[i].Best3Ddist =
            evt.m_allmichels[i].down_clus_michvtx_dist3D;
        evt.m_allmichels[i].best_XZ = evt.m_allmichels[i].down_to_clus_XZ;
        evt.m_allmichels[i].best_UZ = evt.m_allmichels[i].down_to_clus_UZ;
        evt.m_allmichels[i].best_VZ = evt.m_allmichels[i].down_to_clus_VZ;
      } else {
        evt.m_allmichels[i].BestMatch = 0;
        evt.m_allmichels[i].Best3Ddist = 9999.;
        evt.m_allmichels[i].best_XZ = 9999.;
        evt.m_allmichels[i].best_UZ = 9999.;
        evt.m_allmichels[i].best_VZ = 9999.;

        continue;
      }

      // Set some MichelEvent info about second best match (not used)
      if (distances3D[1] == evt.m_allmichels[i].up_to_vertex_dist3D)
        evt.m_allmichels[i].SecondBestMatch = 1;
      else if (distances3D[1] == evt.m_allmichels[i].down_to_vertex_dist3D)
        evt.m_allmichels[i].SecondBestMatch = 2;
      else if (distances3D[1] == evt.m_allmichels[i].up_clus_michel_dist3D)
        evt.m_allmichels[i].SecondBestMatch = 3;
      else if (distances3D[1] == evt.m_allmichels[i].down_clus_michel_dist3D)
        evt.m_allmichels[i].SecondBestMatch = 4;
      else {
        evt.m_allmichels[i].SecondBestMatch = 0;
      }

      // Set more best match info
      if (evt.m_allmichels[i].best_XZ < evt.m_allmichels[i].best_UZ and
          evt.m_allmichels[i].best_XZ < evt.m_allmichels[i].best_VZ) {
        evt.m_allmichels[i].bestview = 1;
      }
      else if (evt.m_allmichels[i].best_UZ < evt.m_allmichels[i].best_XZ and
               evt.m_allmichels[i].best_UZ < evt.m_allmichels[i].best_VZ) {
        evt.m_allmichels[i].bestview = 2;
      }

      else if (evt.m_allmichels[i].best_VZ < evt.m_allmichels[i].best_UZ and
               evt.m_allmichels[i].best_VZ < evt.m_allmichels[i].best_XZ) {
        evt.m_allmichels[i].bestview = 3;
      }
      int matchtype = evt.m_allmichels[i].BestMatch;
      if (matchtype == 1 || matchtype == 3) {
        evt.m_allmichels[i].recoEndpoint = 1;
      }
      else if (matchtype == 2 || matchtype == 4) {
        evt.m_allmichels[i].recoEndpoint = 2;
      }

      // Record this michel, which has passed the cuts
      nmichelspass.push_back(evt.m_allmichels[i]);
      evt.m_nmichelspass.push_back(evt.m_allmichels[i]);
    }

    evt.m_nmichels.clear();

    // Do we have at least one passing michel?
    if (nmichelspass.empty()){
      return false;
    }
    else {
      evt.selection = 1;
      evt.m_nmichels =
          nmichelspass;  // replace vector of michels with the vector of michels
                         // that passed the above cut
      return true;
    }
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

  // Here's the action function that makes the cut
  //
  // Return closestMichel[0].Best3Ddist <= kMAX3DDIST (2600.)
  //
  // ALSO: modify the passed EVENT object and return it as reference.
  static bool MichelRangeCut(const UNIVERSE& univ, EVENT& evt,
                             const double max_dist = kMAX3DDIST) {
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
    }
    for (int i = 0; i < nmichels; ++i) {
      evt.m_nmichels[i].GetPionAngle(univ);
      double dist = evt.m_nmichels[i].Best3Ddist;
      evt.m_nmichelspass.clear();
      if (dist < max_dist) {
        evt.m_nmichelspass.push_back(evt.m_nmichels[i]);
      }
      if (evt.m_nmichels[i].OrderOfMichel == 1) {
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
        evt.m_idx = -1;
      }
    }
    if (closestMichel.empty()) return false;
    double lowtpiinevent = closestMichel[0].pionKE;
    evt.lowTpi = lowtpiinevent;
    evt.pT_reco = univ.GetMuonPT();
    evt.q3_reco = univ.Getq3_lowrecoil();
    evt.eavail_reco = univ.NewEavail();
    evt.m_nmichels.clear();
    evt.m_nmichels = closestMichel;
    if (univ.GetMuonPT() < .20 and univ.NewEavail() < 50. and
        !evt.m_nmichels.empty()) {
      double vtx_t = univ.GetVertex().T() / pow(10, 3);  // mus
    }
    // Lets look at single pi+ events only
    return closestMichel[0].Best3Ddist <= max_dist;
  };

 private:
  using Michel = typename std::remove_reference<decltype(
      std::declval<EVENT>().m_nmichels.front())>::type;
  int michelgroup;
  double m_max_dist;
  bool checkCut(const UNIVERSE& univ, EVENT& evt) const {
    return MichelRangeCut(univ, evt, m_max_dist);
  };
};

}  // namespace LowRecoilPion

#endif  // LOWRECOILPIONCUTS_H
