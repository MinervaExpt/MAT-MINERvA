#ifndef LOWRECOILPIONRECO_H
#define LOWRECOILPIONRECO_H

#include <cmath>
#include <vector>

#include "TVector3.h"
#include "utilities/PlotUtilsPhysicalConstants.h"

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

template <class T>
class Michel {
 public:
  // Constructors
  Michel(const T& univ, const int ci);
  Michel(){};

  // Initialization functions
  void DoMoreInitializationAndProcessing();
  void DoMatching();  // fill in "best matching stuff"

  // Access
  void DoesMichelMatchVtx(const T& univ);   // Gets info for Vtx Match
  void DoesMichelMatchClus(const T& univ);  // Gets info for ClusterMatch
  void GetBestMatch();  // get type for best match out of all four saved matches
  void GetPionAngle(
      const T& univ);  // WRT beam; angle btwn best mich endpt and vtx

  // Generic Info
  std::vector<double> up_location;    // upstream location 0 X 1 U 2 V 3 Z
  std::vector<double> down_location;  // downstream location
  double m_x1 = 9999.;                // Michel Endpoint 1 x
  double m_x2 = 9999.;                // Michel Endpoint 2 x
  double m_y1 = 9999.;                // Michel Endpoint 1 y
  double m_y2 = 9999.;                // Michel Endpoint 2 y
  double m_u1 = 9999.;                // Michel Endpoint 1 u
  double m_u2 = 9999.;                // Michel Endpoint 2 u
  double m_v1 = 9999.;                // Michel Endpoint 1 v
  double m_v2 = 9999.;                // Michel Endpoint 2 v
  double m_z1 = 9999.;                // Mihel Endpoint z 1
  double m_z2 = 9999.;                // Michel Endpoiint z2
  double energy = -999.;              // Michel energy
  double time = -999.;                // Michel Time
  int is_fitted = -1;                 // Is the Michel fitted? 0 no. 1 yes.

  // Following are 2D distances (were stored in vectors but now as explicit data
  // members)
  double up_to_vertex_XZ = 9999.;
  double up_to_vertex_UZ = 9999.;
  double up_to_vertex_VZ = 9999.;
  double down_to_vertex_XZ = 9999.;
  double down_to_vertex_UZ = 9999.;
  double down_to_vertex_VZ = 9999.;
  double down_to_clus_XZ = 9999.;
  double down_to_clus_UZ = 9999.;
  double down_to_clus_VZ = 9999.;
  double up_to_clus_XZ = 9999.;
  double up_to_clus_UZ = 9999.;
  double up_to_clus_VZ = 9999.;

  // Michel End point to Vertex distance  (up = endpoint 1 and down = endpoint
  // 2) TODO: actually find out which end point is upstream or downstream
  double up_to_vertex_dist3D = 9999.;
  double down_to_vertex_dist3D = 9999.;
  // Maybe keep a vector of clusters that matched to each endpoint?
  std::vector<Cluster> cluster_to_up_match;
  std::vector<Cluster> cluster_to_down_match;
  // 3D distances between cluster and Michel
  double up_to_cluster_dist3D = -9999.;  // Distance between vertex and cluster
                                         // that was matched to endpoint 1
  double down_to_cluster_dist3D =
      -9999.;  // Distance between vertex and cluster matched to endpoint 2
  double up_clus_michel_dist3D =
      -9999.;  // Distance between Michel endpoint 1 and clusters
  double down_clus_michel_dist3D =
      -9999.;  // Distance between Michel endpoint 2 and clusters
  double up_clus_michvtx_dist3D =
      -9999.;  // Distance between the Michel end point 1 that matched to
               // clusters and the vertex - this will be used as pion range
  double down_clus_michvtx_dist3D =
      9999.;  // Distance between the Michel endpoint 2 that matched to clusters
              // and the vertex - this will be used as pion range
  double vtx_michel_timediff = 9999.;

  // Overlay fraction of the Michel, Default if Data. 0 if MC 1 if
  // Data... (maybe some events in between?)
  double overlay_fraction = -1.0;
  int nclusters = -1;  // number of (non-muon) clusters in the primary event
  int nnonmuclusters = 0.0;
  int vtx_endpoint = 0;   // 1 or 2 for which Michel end point is closest
  int clus_endpoint = 0;  // 1 or 2 for which Michel endpoint is closest

  // best matching stuff
  // enum *best_cluster_match; // just tells you which of the four matches are
  // the best match enum BestClusterMatch {kUpVtxMatch, kDownVtxMatch,
  // kUpClusMatch, kDownClusMatch, kNClusterMatches}; the following is in place
  // until i figure out how to use the enum.

  // 0 = null match, 1= kUpVtxMatch, 2 = kDownVtxMatch, 3 = kUpClusMatch, 4 =
  // kDownClusMatch,
  int BestMatch = 0;
  int SecondBestMatch = 0;
  int tuple_idx;  // index of this out of all the michels saved in the tuple
  double Best3Ddist = 9999.;  // Best 3D dist out of either vtx or clus match

  // best 2D distance for the best type of match
  double best_XZ = -9999.;
  double best_UZ = -9999.;
  double best_VZ = -9999.;
  int bestview = -1;  // 1 = XZ, 2 = UZ, 3 = VZ

  // Want to save the index of the clusters that the Michel best matched to.
  int xclus_idx = -1;
  int uclus_idx = -1;
  int vclus_idx = -1;

  // True initial position of the michel  TODO: initial the other Michel truth
  // member data here  (energy, time, momentum etc)
  double true_angle = -9999.;
  double true_phi = -9999.;
  double true_initialx = -9999.;
  double true_initialy = -9999.;
  double true_initialz = -9999.;
  double true_e = -9999.;
  double true_p = -9999.;
  double true_pdg = -1.0;
  int true_parentid = -1;
  int true_parentpdg = -1;
  double true_parent_energy = -9999.;
  double true_parent_p = -9999.;
  double true_parent_px = -9999.;
  double true_parent_py = -9999.;
  double true_parent_pz = -9999.;
  double true_parent_xi = -9999.;
  double true_parebt_yi = -9999.;
  double true_parent_zi = -9999.;
  double true_parent_xf = -9999.;
  double true_parent_yf = -9999.;
  double true_parent_zf = -9999.;

  // the following member data were created to investigate my weird convoluted
  // way of geting x and y values. Probably dont need them now. // TODO: check
  // and remove the following member data
  double best_angle = -9999.;
  double best_phi = -9999.;
  double reco_Epi = -9999.;
  double reco_KE = -9999.;
  double reco_ppi = -9999.;
  double reco_ppix = -9999.;
  double reco_ppiy = -9999.;
  double reco_ppiz = -9999.;
  double up_clus_x = -9999.;
  double up_clus_y = -9999.;
  double up_clus_z = -9999.;
  double down_clus_x = -9999.;
  double down_clus_y = -9999.;
  double down_clus_z = -9999.;
  double up_vtx_x = -9999.;
  double up_vtx_y = -9999.;
  double up_vtx_z = -9999.;
  double down_vtx_x = -9999.;
  double down_vtx_y = -9999.;
  double down_vtx_z = -9999.;
  double is_overlay = -1;
  double pionKE = 99999.;

  // Adding the following member data to determine the true endpoint position of
  // the Michel
  // 0 = Overlay Michel, 1 = Endpoint 1 is correct intial position of
  // Michel, 2 = Endpoint 2 is correct Initial Position of Michel
  int trueEndpoint = -1;

  // ~DeleteMichel(){delete this;};  // This is going to be the main destructor.

  // -1 is default/NULL like above. 1 = Endpoint 1 is better match, 2 =
  // Endpoint 2 is better match
  int recoEndpoint = -1;

  // Gives the order of the michel in the event. Starting from 1 being the
  // closest to the vertex
  int OrderOfMichel = -1;

  // This vector will contain a value for each match type -1 or the
  // matchtype 1, 2, 3, 4 for UpVtx, DownVTx, Upclus, DownClus
  // depending of that match type passes our 2D distance cut. (if
  // distance is large, then it'll pass all of them).
  std::vector<int> passable_matchtype{-1, -1, -1, -1};
};

}  // namespace LowRecoilPion

#include "utilities/LowRecoilPionReco.cxx"

#endif
