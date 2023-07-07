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
  //
  // TODO Also sets some stuff in the MichelEvent, (e.g. muon PT, q3, eavail)
  // which does NOT need to be in here.
  static EVENT GetQualityMichels(const UNIVERSE& univ,
                                 const EVENT& _evt = EVENT()) {
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
    if (univ.ShortName() != "cv" && !evt.m_allmichels.empty()) return true;

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
  // It compares distance between Michel and whatever its best match is to find
  // the Best type of Michel for a single Michel.
  static EVENT Apply2DDistanceCut(const UNIVERSE& univ, const EVENT& evt,
                                  const double max_dist = kMAX2DDIST) {
    EVENT ret_evt = evt;

    // Update our MichelEvent ret_evt AND fill a vector of passing michels
    std::vector<Michel> nmichelspass;
    for (unsigned int i = 0; i < ret_evt.m_allmichels.size(); i++) {
      int upvtxmatch = 0;
      int downvtxmatch = 0;
      int upclusmatch = 0;
      int downclusmatch = 0;

      // For Vertex Match Check to see if 2D distance cut will
      double upvtxXZ = ret_evt.m_allmichels[i].up_to_vertex_XZ;
      double downvtxXZ = ret_evt.m_allmichels[i].down_to_vertex_XZ;
      double upvtxUZ = ret_evt.m_allmichels[i].up_to_vertex_UZ;
      double downvtxUZ = ret_evt.m_allmichels[i].down_to_vertex_UZ;
      double upvtxVZ = ret_evt.m_allmichels[i].up_to_vertex_VZ;
      double downvtxVZ = ret_evt.m_allmichels[i].down_to_vertex_VZ;

      std::vector<double> upvtx = {upvtxXZ, upvtxUZ, upvtxVZ};
      std::vector<double> downvtx = {upvtxXZ, upvtxUZ, upvtxVZ};

      std::sort(upvtx.begin(), upvtx.end());
      std::sort(downvtx.begin(), downvtx.end());

      // TODO this can all just be a big OR.
      if (upvtxXZ < max_dist && (upvtxUZ < max_dist || upvtxVZ < max_dist))
        ret_evt.m_allmichels[i].passable_matchtype.at(0) = 1;
      else if (upvtxUZ < max_dist && (upvtxXZ < max_dist || upvtxVZ < max_dist))
        ret_evt.m_allmichels[i].passable_matchtype.at(0) = 1;
      else if (upvtxVZ < max_dist && (upvtxXZ < max_dist || upvtxUZ < max_dist))
        ret_evt.m_allmichels[i].passable_matchtype.at(0) = 1;
      else
        ret_evt.m_allmichels[i].passable_matchtype.at(0) = -1;

      // TODO this can all just be a big OR.
      if (downvtxXZ < max_dist &&
          (downvtxUZ < max_dist || downvtxVZ < max_dist))
        ret_evt.m_allmichels[i].passable_matchtype.at(1) = 2;
      else if (downvtxUZ < max_dist &&
               (downvtxXZ < max_dist || downvtxVZ < max_dist))
        ret_evt.m_allmichels[i].passable_matchtype.at(1) = 2;
      else if (downvtxVZ < max_dist &&
               (downvtxXZ < max_dist || downvtxUZ < max_dist))
        ret_evt.m_allmichels[i].passable_matchtype.at(1) = 2;
      else
        ret_evt.m_allmichels[i].passable_matchtype.at(1) = -1;

      double upclusXZ = ret_evt.m_allmichels[i].up_to_clus_XZ;
      double upclusUZ = ret_evt.m_allmichels[i].up_to_clus_UZ;
      double upclusVZ = ret_evt.m_allmichels[i].up_to_clus_VZ;
      double downclusXZ = ret_evt.m_allmichels[i].down_to_clus_XZ;
      double downclusUZ = ret_evt.m_allmichels[i].down_to_clus_UZ;
      double downclusVZ = ret_evt.m_allmichels[i].down_to_clus_VZ;

      std::vector<double> upclus = {upclusXZ, upclusUZ, upclusVZ};
      std::vector<double> downclus = {downclusXZ, downclusUZ, downclusVZ};

      std::sort(upclus.begin(), upclus.end());
      std::sort(downclus.begin(), downclus.end());

      // TODO this can all just be a big OR.
      if (upclusXZ < max_dist && (upclusUZ < max_dist || upclusVZ < max_dist))
        ret_evt.m_allmichels[i].passable_matchtype.at(2) = 3;
      else if (upclusUZ < max_dist &&
               (upclusXZ < max_dist || upclusVZ < max_dist))
        ret_evt.m_allmichels[i].passable_matchtype.at(2) = 3;
      else if (upclusVZ < max_dist &&
               (upclusXZ < max_dist || upclusUZ < max_dist))
        ret_evt.m_allmichels[i].passable_matchtype.at(2) = 3;
      else
        ret_evt.m_allmichels[i].passable_matchtype.at(2) = -1;

      // TODO this can all just be a big OR.
      if (downclusXZ < max_dist &&
          (downclusUZ < max_dist || downclusVZ < max_dist))
        ret_evt.m_allmichels[i].passable_matchtype.at(3) = 4;
      else if (downclusUZ < max_dist &&
               (downclusXZ < max_dist || downclusVZ < max_dist))
        ret_evt.m_allmichels[i].passable_matchtype.at(3) = 4;
      else if (downclusVZ < max_dist &&
               (downclusXZ < max_dist || downclusUZ < max_dist))
        ret_evt.m_allmichels[i].passable_matchtype.at(3) = 4;
      else
        ret_evt.m_allmichels[i].passable_matchtype.at(3) = -1;

      // CUT
      // If we don't have acceptable up and down vtx and cluster distances for
      // this Michel, it fails.
      // TODO condense? all? loop?
      if (ret_evt.m_allmichels[i].passable_matchtype.at(1) == -1 &&
          ret_evt.m_allmichels[i].passable_matchtype.at(2) == -1 &&
          ret_evt.m_allmichels[i].passable_matchtype.at(3) == -1 &&
          ret_evt.m_allmichels[i].passable_matchtype.at(0) == -1)
        continue;

      // Record 3D distance
      std::vector<double> distances3D;

      // TODO This could be condensed? For loop? case-switch?
      if (ret_evt.m_allmichels[i].passable_matchtype[0] == 1)
        distances3D.push_back(
            ret_evt.m_allmichels[i]
                .up_to_vertex_dist3D);  // Distance between michel to vertex
      if (ret_evt.m_allmichels[i].passable_matchtype[1] == 2)
        distances3D.push_back(
            ret_evt.m_allmichels[i]
                .down_to_vertex_dist3D);  // distancebetween michel to vertex
      if (ret_evt.m_allmichels[i].passable_matchtype[2] == 3)
        distances3D.push_back(
            ret_evt.m_allmichels[i]
                .up_clus_michel_dist3D);  // distnace between michel to cluster
      if (ret_evt.m_allmichels[i].passable_matchtype[3] == 4)
        distances3D.push_back(
            ret_evt.m_allmichels[i]
                .down_clus_michel_dist3D);  // distance between
                                            // michel to cluster
      if (distances3D.empty()) distances3D = {9999., 9999., 9999., 9999.};

      // Sort these 3D distance in ascending order.
      std::sort(distances3D.begin(), distances3D.end());

      // CUT
      // must have at least one good 3D distance between michel and cluster and
      // michel and vertex. Else this michel fails.
      if (distances3D[0] == 9999.) continue;

      // Set some MichelEvent info about this best match.
      // There's a continue here, but I think it's just an assert/sanity check.
      // TODO: condense? case-switch?
      if (distances3D[0] == ret_evt.m_allmichels[i].up_to_vertex_dist3D) {
        ret_evt.m_allmichels[i].best_XZ =
            ret_evt.m_allmichels[i].up_to_vertex_XZ;
        ret_evt.m_allmichels[i].best_UZ =
            ret_evt.m_allmichels[i].up_to_vertex_UZ;
        ret_evt.m_allmichels[i].best_VZ =
            ret_evt.m_allmichels[i].up_to_vertex_VZ;
        ret_evt.m_allmichels[i].BestMatch = 1;
        ret_evt.m_allmichels[i].Best3Ddist =
            ret_evt.m_allmichels[i].up_to_vertex_dist3D;

        // std::cout << "This  Michel is UPVTX and has true endpoint " <<
        // ret_evt.m_allmichels[i].trueEndpoint << std::endl;
      } else if (distances3D[0] ==
                 ret_evt.m_allmichels[i].down_to_vertex_dist3D) {
        ret_evt.m_allmichels[i].BestMatch = 2;
        ret_evt.m_allmichels[i].best_XZ =
            ret_evt.m_allmichels[i].down_to_vertex_XZ;
        ret_evt.m_allmichels[i].best_UZ =
            ret_evt.m_allmichels[i].down_to_vertex_UZ;
        ret_evt.m_allmichels[i].best_VZ =
            ret_evt.m_allmichels[i].down_to_vertex_VZ;
        ret_evt.m_allmichels[i].Best3Ddist =
            ret_evt.m_allmichels[i].down_to_vertex_dist3D;
        // std::cout << "This  Michel is DOWNVTX and has true endpoint " <<
        // ret_evt.m_allmichels[i].trueEndpoint << std::endl;
      } else if (distances3D[0] ==
                 ret_evt.m_allmichels[i].up_clus_michel_dist3D) {
        ret_evt.m_allmichels[i].BestMatch = 3;
        ret_evt.m_allmichels[i].best_XZ = ret_evt.m_allmichels[i].up_to_clus_XZ;
        ret_evt.m_allmichels[i].best_UZ = ret_evt.m_allmichels[i].up_to_clus_VZ;
        ret_evt.m_allmichels[i].best_VZ = ret_evt.m_allmichels[i].up_to_clus_UZ;
        ret_evt.m_allmichels[i].Best3Ddist =
            ret_evt.m_allmichels[i].up_clus_michvtx_dist3D;
        // std::cout << "This  Michel is UPCLUS and has true endpoint " <<
        // ret_evt.m_allmichels[i].trueEndpoint << std::endl;
      } else if (distances3D[0] ==
                 ret_evt.m_allmichels[i].down_clus_michel_dist3D) {
        ret_evt.m_allmichels[i].BestMatch = 4;
        ret_evt.m_allmichels[i].Best3Ddist =
            ret_evt.m_allmichels[i].down_clus_michvtx_dist3D;
        ret_evt.m_allmichels[i].best_XZ =
            ret_evt.m_allmichels[i].down_to_clus_XZ;
        ret_evt.m_allmichels[i].best_UZ =
            ret_evt.m_allmichels[i].down_to_clus_UZ;
        ret_evt.m_allmichels[i].best_VZ =
            ret_evt.m_allmichels[i].down_to_clus_VZ;
      } else {
        ret_evt.m_allmichels[i].BestMatch = 0;
        ret_evt.m_allmichels[i].Best3Ddist = 9999.;
        ret_evt.m_allmichels[i].best_XZ = 9999.;
        ret_evt.m_allmichels[i].best_UZ = 9999.;
        ret_evt.m_allmichels[i].best_VZ = 9999.;

        continue;
      }

      // Set some MichelEvent info about second best match (not used)
      if (distances3D[1] == ret_evt.m_allmichels[i].up_to_vertex_dist3D)
        ret_evt.m_allmichels[i].SecondBestMatch = 1;
      else if (distances3D[1] == ret_evt.m_allmichels[i].down_to_vertex_dist3D)
        ret_evt.m_allmichels[i].SecondBestMatch = 2;
      else if (distances3D[1] == ret_evt.m_allmichels[i].up_clus_michel_dist3D)
        ret_evt.m_allmichels[i].SecondBestMatch = 3;
      else if (distances3D[1] ==
               ret_evt.m_allmichels[i].down_clus_michel_dist3D)
        ret_evt.m_allmichels[i].SecondBestMatch = 4;
      else {
        ret_evt.m_allmichels[i].SecondBestMatch = 0;
      }

      // Set more best match info
      if (ret_evt.m_allmichels[i].best_XZ < ret_evt.m_allmichels[i].best_UZ and
          ret_evt.m_allmichels[i].best_XZ < ret_evt.m_allmichels[i].best_VZ) {
        ret_evt.m_allmichels[i].bestview = 1;
      } else if (ret_evt.m_allmichels[i].best_UZ <
                     ret_evt.m_allmichels[i].best_XZ and
                 ret_evt.m_allmichels[i].best_UZ <
                     ret_evt.m_allmichels[i].best_VZ) {
        ret_evt.m_allmichels[i].bestview = 2;
      }

      else if (ret_evt.m_allmichels[i].best_VZ <
                   ret_evt.m_allmichels[i].best_UZ and
               ret_evt.m_allmichels[i].best_VZ <
                   ret_evt.m_allmichels[i].best_XZ) {
        ret_evt.m_allmichels[i].bestview = 3;
      }
      int matchtype = ret_evt.m_allmichels[i].BestMatch;
      if (matchtype == 1 || matchtype == 3) {
        ret_evt.m_allmichels[i].recoEndpoint = 1;
      } else if (matchtype == 2 || matchtype == 4) {
        ret_evt.m_allmichels[i].recoEndpoint = 2;
      }

      // Record this michel, which has passed the cuts
      nmichelspass.push_back(ret_evt.m_allmichels[i]);
      ret_evt.m_nmichelspass.push_back(ret_evt.m_allmichels[i]);
    }

    ret_evt.m_nmichels.clear();

    if (!nmichelspass.empty()) {
      ret_evt.selection = 1;
      ret_evt.m_nmichels =
          nmichelspass;  // replace vector of michels with the vector of michels
                         // that passed the above cut
    }

    return ret_evt;
  }

 private:
  using Michel = typename std::remove_reference<decltype(
      std::declval<EVENT>().m_nmichels.front())>::type;

  double m_maxDistance;  // Maximum distance from the vertex that the best
                         // Michel can have in mm

  bool checkCut(const UNIVERSE& univ, EVENT& evt) const {
    if (evt.m_allmichels.size() == 0) return false;

    // Update michels. Only those that pass 2D cut remain.
    evt = Apply2DDistanceCut(univ, evt, m_maxDistance);

    return evt.m_nmichels.empty()
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
  static EVENT Apply3DDistanceCut(const UNIVERSE& univ, const EVENT& evt,
                                  const double max_dist = kMAX3DDIST) {
    EVENT ret_evt = evt;
    int nmichels = ret_evt.m_nmichels.size();
    int noverlay = 0.0;
    std::vector<Michel> closestMichel;
    ret_evt.m_bestdist = 9999.;  // setting some default value for best distance
    std::vector<double> allmichel3Ddist;
    for (int i = 0; i < nmichels; ++i) {
      double dist =
          ret_evt.m_nmichels[i].Best3Ddist;  // getting the minimum pion range
                                             // (vertex to Michel/Clus distance)
      allmichel3Ddist.push_back(dist);
    }  // Get all the distances for the michels that pass the 2D dist cut

    std::sort(allmichel3Ddist.begin(), allmichel3Ddist.end());

    for (int i = 0; i < nmichels; ++i) {
      double dist = ret_evt.m_nmichels[i].Best3Ddist;
      int order = i + 1;
      if (dist == allmichel3Ddist[i])
        ret_evt.m_nmichels[i].OrderOfMichel = order;
    }
    for (int i = 0; i < nmichels; ++i) {
      ret_evt.m_nmichels[i].GetPionAngle(univ);
      double dist = ret_evt.m_nmichels[i].Best3Ddist;
      ret_evt.m_nmichelspass.clear();
      if (dist < max_dist) {
        ret_evt.m_nmichelspass.push_back(ret_evt.m_nmichels[i]);
      }
      if (ret_evt.m_nmichels[i].OrderOfMichel == 1) {
        ret_evt.m_bestdist = dist;
        ret_evt.m_idx = i;
        if (ret_evt.m_nmichels[i].overlay_fraction > 0.5)
          ret_evt.ClosestMichelsIsOverlay = 1;
        ret_evt.m_best_XZ = ret_evt.m_nmichels[i].best_XZ;
        ret_evt.m_best_UZ = ret_evt.m_nmichels[i].best_UZ;
        ret_evt.m_best_VZ = ret_evt.m_nmichels[i].best_VZ;
        ret_evt.m_matchtype = ret_evt.m_nmichels[i].BestMatch;
        int bmatch = ret_evt.m_nmichels[i].BestMatch;
        if (bmatch == 1 || bmatch == 3) {
          ret_evt.best_x = ret_evt.m_nmichels[i].m_x1;
          ret_evt.best_y = ret_evt.m_nmichels[i].m_y1;
          ret_evt.best_z = ret_evt.m_nmichels[i].m_z1;
        } else if (bmatch == 2 || bmatch == 4) {
          ret_evt.best_x = ret_evt.m_nmichels[i].m_x2;
          ret_evt.best_y = ret_evt.m_nmichels[i].m_y2;
          ret_evt.best_z = ret_evt.m_nmichels[i].m_z2;
        }
        ret_evt.b_truex = ret_evt.m_nmichels[i].true_initialx;
        ret_evt.b_truey = ret_evt.m_nmichels[i].true_initialy;
        ret_evt.b_truez = ret_evt.m_nmichels[i].true_initialz;
        closestMichel.push_back(ret_evt.m_nmichels[i]);
      } else {
        ret_evt.m_idx = -1;
      }
    }

    ret_evt.m_nmichels.clear();
    ret_evt.m_nmichels = closestMichel;

    if (ret_evt.m_nmichels.empty()) return ret_evt;

    ret_evt.lowTpi = ret_evt.m_nmichels[0].pionKE;
    ret_evt.pT_reco = univ.GetMuonPT();
    ret_evt.q3_reco = univ.Getq3_lowrecoil();
    ret_evt.eavail_reco = univ.NewEavail();

    return ret_evt;
  };

 private:
  using Michel = typename std::remove_reference<decltype(
      std::declval<EVENT>().m_nmichels.front())>::type;
  int michelgroup;
  double m_max_dist;
  bool checkCut(const UNIVERSE& univ, EVENT& evt) const {
    if (evt.m_nmichels.size() == 0) return false;

    evt = Apply3DDistanceCut(univ, evt, m_max_dist);

    if (evt.m_nmichels.empty())
      return false;
    else
      return evt.m_nmichels[0].Best3Ddist <= max_dist;
  };
};

}  // namespace LowRecoilPion

#endif  // LOWRECOILPIONCUTS_H
