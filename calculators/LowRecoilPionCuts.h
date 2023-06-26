#ifndef LOWRECOILPIONCUTS_H
#define LOWRECOILPIONCUTS_H

#include "PlotUtils/Cut.h"

namespace LowRecoilPion {

const double kMAX3DDIST = 2600.;

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
    evt.q3_reco = univ.Getq3();
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
