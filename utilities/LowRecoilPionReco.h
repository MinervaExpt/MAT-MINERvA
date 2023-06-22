#ifndef LowRecoilPionReco_h
#define LowRecoilPionReco_h

namespace LowRecoilPion {

struct Cluster {
  double ecalo;   // Passive corrected energy (calorimetric)
  double energy;  // MeV
  double pos;     // mm
  double time;    // microsec
  double zpos;    // mm
  int cluster_idx;
  int ismuon;  // 1 is on muon track, 0 is off muon track
  int subdet;
  int view;         // 1 = x, 2 = u, 3 = v
  bool is_quality;  // not filled in initialization

  Cluster(double _ecalo, double _energy, double _pos, double _time,
          double _zpos, int _cluster_idx, int _ismuon, int _subdet, int _view)
      : ecalo(_ecalo),
        energy(_energy),
        pos(_pos),
        time(_time),
        zpos(_zpos),
        cluster_idx(_cluster_idx),
        ismuon(_ismuon),
        subdet(_subdet),
        view(_view){};

  Cluster(){};

  template <class T>
  Cluster(const T& univ, const int& ci)
      : ecalo(univ.GetVecElem("cluster_ecalo", ci)),
        energy(univ.GetVecElem("cluster_energy", ci)),
        pos(univ.GetVecElem("cluster_pos", ci)),
        time(univ.GetVecElem("cluster_time", ci) / pow(10, 3)),
        zpos(univ.GetVecElem("cluster_z", ci)),
        cluster_idx(ci),
        ismuon(univ.GetVecElem("cluster_isMuontrack", ci)),
        subdet(univ.GetVecElem("cluster_subdet", ci)),
        view(univ.GetVecElem("cluster_view", ci)){};
};

}  // namespace LowRecoilPion

#endif
