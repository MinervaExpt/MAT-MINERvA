#ifndef LOWRECOILPIONRECO_CXX
#define LOWRECOILPIONRECO_CXX

#include "utilities/LowRecoilPionReco.h"

namespace LowRecoilPion {

// Setting the Michel data members directly from CV universe. No calculations
// made at this stage
template <class T>
Michel<T>::Michel(const T& univ, int ci) {
  energy = univ.GetVecElem("FittedMichel_michel_energy", ci);
  time = univ.GetVecElem("FittedMichel_michel_time", ci) / pow(10, 3);
  is_fitted = univ.GetVecElemInt("FittedMichel_michel_fitPass", ci);
  up_location.push_back(univ.GetVecElem("FittedMichel_michel_x1", ci));
  up_location.push_back(univ.GetVecElem("FittedMichel_michel_u1", ci));
  up_location.push_back(univ.GetVecElem("FittedMichel_michel_v1", ci));
  up_location.push_back(univ.GetVecElem("FittedMichel_michel_z1", ci));
  down_location.push_back(univ.GetVecElem("FittedMichel_michel_x2", ci));
  down_location.push_back(univ.GetVecElem("FittedMichel_michel_u2", ci));
  down_location.push_back(univ.GetVecElem("FittedMichel_michel_v2", ci));
  down_location.push_back(univ.GetVecElem("FittedMichel_michel_z2", ci));
  m_x1 = univ.GetVecElem("FittedMichel_michel_x1", ci);
  m_y1 = univ.GetVecElem("FittedMichel_michel_y1", ci);
  m_u1 = univ.GetVecElem("FittedMichel_michel_u1", ci);
  m_v1 = univ.GetVecElem("FittedMichel_michel_v1", ci);
  m_z1 = univ.GetVecElem("FittedMichel_michel_z1", ci);
  m_x2 = univ.GetVecElem("FittedMichel_michel_x2", ci);
  m_y2 = univ.GetVecElem("FittedMichel_michel_y2", ci);
  m_u2 = univ.GetVecElem("FittedMichel_michel_u2", ci);
  m_z2 = univ.GetVecElem("FittedMichel_michel_z2", ci);
  m_v2 = univ.GetVecElem("FittedMichel_michel_v2", ci);
  nclusters = univ.GetInt("cluster_view_sz");
  overlay_fraction = univ.GetVecElem("FittedMichel_michel_datafraction", ci);

  if (univ.IsTruth()) {
    true_initialx = univ.GetVecElem(
        "truth_FittedMichel_reco_micheltrajectory_initialx", ci);
    true_initialy = univ.GetVecElem(
        "truth_FittedMichel_reco_micheltrajectory_initialy", ci);
    true_initialz = univ.GetVecElem(
        "truth_FittedMichel_reco_micheltrajectory_initialz", ci);
    is_overlay = univ.GetVecElemInt("FittedMichel_michel_isoverlay", ci);
    true_e =
        univ.GetVecElem("truth_FittedMichel_reco_micheltrajectory_energy", ci);
    true_pdg =
        univ.GetVecElemInt("truth_FittedMichel_reco_micheltrajectory_pdg", ci);
    true_parentpdg =
        univ.GetVecElemInt("truth_FittedMichel_true_primaryparent_pdg", ci);
    true_parentid =
        univ.GetVecElemInt("truth_FittedMichel_true_primaryparent_trackID", ci);
    true_p = univ.GetVecElem(
        "truth_FittedMichel_reco_micheltrajectory_momentum", ci);
    true_parent_energy =
        univ.GetVecElem("truth_FittedMichel_true_primaryparent_energy", ci);
    true_parent_p =
        univ.GetVecElem("truth_FittedMichel_true_primaryparent_momentum", ci);
    double true_parentp =
        univ.GetVecElem("truth_FittedMichel_true_primaryparent_momentum", ci);
    double true_parente =
        univ.GetVecElem("truth_FittedMichel_true_primaryparent_energy", ci);
    double mass = sqrt(pow(true_parente, 2) - pow(true_parentp, 2));
    pionKE = true_parente - mass;
    TVector3 parentfinalpos(
        univ.GetVecElem("truth_FittedMichel_true_primaryparent_finalx", ci),
        univ.GetVecElem("truth_FittedMichel_true_primaryparent_finaly", ci),
        univ.GetVecElem("truth_FittedMichel_true_primaryparent_finalz", ci));

    TVector3 parentinitialpos(
        univ.GetVecElem("truth_FittedMichel_true_primaryparent_initialx", ci),
        univ.GetVecElem("truth_FittedMichel_true_primaryparent_initialy", ci),
        univ.GetVecElem("truth_FittedMichel_true_primaryparent_initialz", ci));

    TVector3 dispvec = parentfinalpos - parentinitialpos;

    // TVector3 vtxtruepos(univ.GetTrueIntVtxX(), univ.GetTrueIntVtxY(),
    // univ.GetTrueIntVtxZ());

    // TVector3 positionvec = parentfinalpos - vtxtruepos;
    // true_angle = univ.thetaWRTBeam(positionvec.X(), positionvec.Y(),
    // positionvec.Z()); true_phi = univ.phiWRTBeam(positionvec.X(),
    // positionvec.Y(), positionvec.Z());
    dispvec.RotateX(MinervaUnits::numi_beam_angle_rad);

    double parent_px =
        univ.GetVecElem("truth_FittedMichel_true_primaryparent_momentumx", ci);
    double parent_py =
        univ.GetVecElem("truth_FittedMichel_true_primaryparent_momentumy", ci);
    double parent_pz =
        univ.GetVecElem("truth_FittedMichel_true_primaryparent_momentumz", ci);
    const double numi_beam_angle_rad = -0.05887;
    double pyp = -1.0 * sin(numi_beam_angle_rad) * parent_pz +
                 cos(numi_beam_angle_rad) * parent_py;
    double pzp = cos(numi_beam_angle_rad) * parent_pz +
                 sin(numi_beam_angle_rad) * parent_py;

    TVector3 truep(true_parent_px, true_parent_py, true_parent_pz);
    double true_theta = univ.thetaWRTBeam(
        parent_px, parent_py, parent_pz);  // Hopefully this is with respect to
                                           // the dang beam; //truep.Theta();
    // int sign = (parent_px < 0.) ? -1:1;

    // TVector3 unitvec = dispvec.Unit();
    truep.RotateX(MinervaUnits::numi_beam_angle_rad);
    true_phi = dispvec.Phi();    // univ.phiWRTBeam(parent_px, parent_py,
                                 // parent_pz); // Radians
    true_angle = true_theta;     //*TMath::RadToDeg();
    true_parent_px = parent_px;  // true_parent_p*sin(true_angle)*cos(true_phi);
    true_parent_py = pyp;        // true_parent_p*sin(true_angle)*sin(true_phi);
    true_parent_pz = pzp;        // true_parent_p*cos(true_angle);

    double end1diff = abs(
        true_initialz -
        m_z1);  // This gives a value for determining how close the
                // reconstructed endpoint of the michel is to the true intial
                // endpoint (the start point of where the michel decayed from)
    double end2diff =
        abs(true_initialz -
            m_z2);  // this is for endpoint 2. If you compare this to the
                    // endpoint that gets matched to a verted or cluster, you
                    // can determine which type of match ends up getting
                    // correcct matches or wrong matches.

    if (overlay_fraction > 0.5)
      trueEndpoint = 0;
    else if (true_parentpdg == 211 && end1diff < end2diff)
      trueEndpoint = 1;
    else if (true_parentpdg == 211 && end2diff < end1diff)
      trueEndpoint = 2;
  }
  if (is_fitted == 1) {        // Do theMatching for Fitted Michels
    DoesMichelMatchVtx(univ);  // GEts info for Vtx Match
    if (nclusters > 0) DoesMichelMatchClus(univ);  // Gets info for ClusterMatch
    // GetBestMatch();   // Needs to be commented out since this part is done by
    // the cut BestDistance2D.h GetPionAngle(univ);
  }
  // else (this->delete;)
}

// Currently, this is only called in the GetClosestMichelCut
//
// This function will get the angle between Michel endpoint that was matched
// and the vertex
// NONCONST -- SETS AND READS PROPERTIES OF this
template <class T>
void Michel<T>::GetPionAngle(const T& univ) {
  double vtx_x = univ.GetVertex().X();  // mm
  double vtx_y = univ.GetVertex().Y();  // mm
  double vtx_z = univ.GetVertex().Z();  // mm

  TVector3 vtx(vtx_x, vtx_y, vtx_z);
  TVector3 endpoint;

  if (this->BestMatch == 1 || this->BestMatch == 3)
    endpoint.SetXYZ(this->m_x1, this->m_y1, this->m_z1);
  else if (this->BestMatch == 2 || this->BestMatch == 4)
    endpoint.SetXYZ(this->m_x2, this->m_y2, this->m_z2);
  else
    endpoint.SetXYZ(9999., 9999., 9999.);

  TVector3 range = endpoint - vtx;
  double angle = univ.thetaWRTBeam(range.x(), range.y(), range.z());
  this->best_angle = angle;  // in Radians   //*TMath::RadToDeg(); // in Degrees

  double xyp =
      -1.0 * sin(MinervaUnits::numi_beam_angle_rad + 0.00000001) * range.z() +
      cos(MinervaUnits::numi_beam_angle_rad + 0.00000001) * range.y();
  double phiangle = std::atan2(xyp, range.x());
  double phi = phiangle;  // range.Phi();//univ.phiWRTBeam(range.x(), range.y(),
                          // range.z());
  this->best_phi = phi;   // in Radians

  double tpi = univ.GetTpiFromRange(this->Best3Ddist);
  double Epi = tpi + 139.57;
  // double ppi = sqrt(2*tpi*139.57 + tpi*tpi); //June 10 2023 Definition
  double ppi = sqrt(2 * tpi * 139.57 + tpi * tpi);
  double pz = ppi * cos(angle);
  double px = ppi * sin(angle) * cos(phi);
  double py = ppi * sin(angle) * sin(phi);
  this->reco_Epi = Epi;
  this->reco_KE = tpi;
  this->reco_ppi = ppi;
  this->reco_ppix = px;
  this->reco_ppiy = py;
  this->reco_ppiz = pz;
}

// NONCONST -- SETS AND READS PROPERTIES OF this
template <class T>
void Michel<T>::DoesMichelMatchVtx(const T& univ) {
  // Getting Vertex Information
  double vtx_x = univ.GetVertex().X();               // mm
  double vtx_y = univ.GetVertex().Y();               // mm
  double vtx_z = univ.GetVertex().Z();               // mm
  double vtx_t = univ.GetVertex().T() / pow(10, 3);  // mus
  double vtx_u = (0.5 * (vtx_x - sqrt(3.) * vtx_y));
  double vtx_v = (0.5 * (vtx_x + sqrt(3.) * vtx_y));

  // Initializing all the distance comparisons I will need to make
  double zdiff1 = vtx_z - this->m_z1;
  double zdiff2 = vtx_z - this->m_z2;
  double xdiff = 9999.;
  double udiff = 9999.;
  double vdiff = 9999.;
  double XZdist = 9999.;
  double UZdist = 9999.;
  double VZdist = 9999.;

  double michely1 = this->m_y1;
  double michely2 = this->m_y2;
  double michelx1 = this->m_x1;
  double michelx2 = this->m_z2;
  double michelz1 = this->m_z1;
  double michelz2 = this->m_z2;
  double timediff = (this->time) - vtx_t;
  this->vtx_michel_timediff = timediff;

  // 2D distance calculations for Endpoint 1
  xdiff = abs(vtx_x - this->m_x1);
  udiff = abs(vtx_u - this->m_u1);
  vdiff = abs(vtx_v - this->m_v1);

  XZdist = sqrt(xdiff * xdiff + zdiff1 * zdiff1);
  UZdist = sqrt(udiff * udiff + zdiff1 * zdiff1);
  VZdist = sqrt(vdiff * vdiff + zdiff1 * zdiff1);

  this->up_to_vertex_XZ = XZdist;
  this->up_to_vertex_UZ = UZdist;
  this->up_to_vertex_VZ = VZdist;

  // 2D Distance calculations for endpoint2
  xdiff = abs(vtx_x - this->m_x2);
  udiff = abs(vtx_u - this->m_u2);
  vdiff = abs(vtx_v - this->m_v2);
  XZdist = sqrt(xdiff * xdiff + zdiff2 * zdiff2);
  UZdist = sqrt(udiff * udiff + zdiff2 * zdiff2);
  VZdist = sqrt(vdiff * vdiff + zdiff2 * zdiff2);

  this->down_to_vertex_XZ = XZdist;
  this->down_to_vertex_UZ = UZdist;
  this->down_to_vertex_VZ = VZdist;

  // 3D distance calculations
  double xdiff1 = abs(vtx_x - this->m_x1);
  double xdiff2 = abs(vtx_x - this->m_x2);
  double ydiff1 = abs(vtx_y - this->m_y1);
  double ydiff2 = abs(vtx_y - this->m_y2);

  double dist1 = sqrt(zdiff1 * zdiff1 + xdiff1 * xdiff1 + ydiff1 * ydiff1);
  double dist2 = sqrt(zdiff2 * zdiff2 + xdiff2 * xdiff2 + ydiff2 * ydiff2);

  this->up_to_vertex_dist3D = dist1;
  this->down_to_vertex_dist3D = dist2;

  if (dist1 < dist2)
    this->vtx_endpoint = 1;
  else if (dist2 < dist1)
    this->vtx_endpoint = 2;
}

// NONCONST -- SETS AND READS PROPERTIES OF this
template <class T>
void Michel<T>::DoesMichelMatchClus(const T& univ) {
  // This is where the function for Cluster Matching goes

  // Inititalizing vertex variables needed for cluster matching
  int nclusters = univ.GetInt("cluster_view_sz");
  double vtx_x = univ.GetVertex().X();               // mm
  double vtx_y = univ.GetVertex().Y();               // mm
  double vtx_z = univ.GetVertex().Z();               // mm
  double vtx_t = univ.GetVertex().T() / pow(10, 3);  // mus
  double vtx_u = (0.5 * (vtx_x - sqrt(3.) * vtx_y));
  double vtx_v = (0.5 * (vtx_x + sqrt(3.) * vtx_y));

  double closestdistance1x = 9999.;
  double closestdistance1u = 9999.;
  double closestdistance1v = 9999.;
  double closestdistance1z = 9999.;

  double closestdistance2x = 9999.;
  double closestdistance2u = 9999.;
  double closestdistance2v = 9999.;
  double closestdistance2z = 9999.;

  double michelx1 = this->m_x1;
  double michelx2 = this->m_x2;
  double michelu1 = this->m_u1;
  double michelu2 = this->m_u2;
  double michelv1 = this->m_v1;
  double michelv2 = this->m_v2;
  double michelz1 = this->m_z1;
  double michelz2 = this->m_z2;
  double michely1 = this->m_y1;
  double michely2 = this->m_y2;

  double micheltime = this->time;

  int nonmuclus = 0.0; // not currently used
  std::vector<Cluster> endpoint1_clus;
  std::vector<Cluster> endpoint2_clus;

  // Get the closest distance for each view

  // want to save the index for each closest cluster
  int x1_idx = -1;
  int u1_idx = -1;
  int v1_idx = -1;
  int x2_idx = -1;
  int u2_idx = -1;
  int v2_idx = -1;
  const double minZ = 5980, maxZ = 8422;

  // First loop over Clusters
  for (int i = 0; i < nclusters; i++) {
    // skip processing a cluster whenever possible
    int subdet = univ.GetVecElem("cluster_subdet", i);
    if (subdet != 2 and subdet != 3) continue;
    int ismuon = univ.GetVecElem("cluster_isMuontrack",
                                 i);  // check to make sure cluster is not on
                                      // muon track, 0 is NOT muon, 1 is muon
    if (ismuon != 0) continue;
    double energy = univ.GetVecElem("cluster_energy", i);
    if (energy < 2.) continue;
    double pos = univ.GetVecElem("cluster_pos", i);
    if (abs(pos) > 1000) continue;
    double time = univ.GetVecElem("cluster_time", i) / pow(10, 3);
    double timediff = micheltime - time;
    if (timediff < 0.) continue;
    double zpos = univ.GetVecElem("cluster_z", i);
    int view = univ.GetVecElem("cluster_view", i);

    // nonmuclus++;
    // nnonmuclusters++;

    // if (zpos < minZ || zpos > maxZ) continue; //Require the matched clusters
    // also be in tracker

    // TODO: July 20, 2022 - Change if Kevin says this doesnt make sense

    // if (zpos <  5980 or zpos > 9038) continue;
    // Require match clusters to be in tracker and ecal

    double zdiff1 = abs(zpos - michelz1);
    double zdiff2 = abs(zpos - michelz2);

    // Calculating 2D distance in X view
    if (view == 1) {
      // Endpoint 1 calculations
      double xdiff1 = abs(pos - michelx1);

      double x2Ddistance1 = sqrt(xdiff1 * xdiff1 + zdiff1 * zdiff1);

      // Endpoint 2 Calculations
      double xdiff2 = abs(pos - michelx2);

      double x2Ddistance2 = sqrt(xdiff2 * xdiff2 + zdiff2 * zdiff2);

      if (x2Ddistance1 <= closestdistance1x) {
        closestdistance1x =
            x2Ddistance1;  // this is redundant if I just use index instead
        x1_idx = i;
      }
      if (x2Ddistance2 <= closestdistance2x) {
        closestdistance2x = x2Ddistance2;
        x2_idx = i;
      }
    } else if (view == 2)  // Calculating 2D distance in U view
    {
      // Endpoint 2 Calculations
      double udiff1 = abs(pos - michelu1);

      double u2Ddistance1 = sqrt(udiff1 * udiff1 + zdiff1 * zdiff1);

      // Endpoint 1 Calculations
      double udiff2 = abs(pos - michelu2);

      double u2Ddistance2 = sqrt(udiff2 * udiff2 + zdiff2 * zdiff2);

      if (u2Ddistance1 < closestdistance1u) {
        closestdistance1u = u2Ddistance1;
        u1_idx = i;
      }
      if (u2Ddistance2 < closestdistance2u) {
        closestdistance2u = u2Ddistance2;
        u2_idx = i;
      }
    } else if (view == 3)  // Calculating 2D dsitance in V view
    {
      // Endpoint 1 Calculations
      double vdiff1 = abs(pos - michelv1);

      double v2Ddistance1 = sqrt(vdiff1 * vdiff1 + zdiff1 * zdiff1);
      // Endpoint 2 Calculations
      double vdiff2 = abs(pos - michelv2);

      double v2Ddistance2 = sqrt(vdiff2 * vdiff2 + zdiff2 * zdiff2);

      if (v2Ddistance1 < closestdistance1v) {
        closestdistance1v = v2Ddistance1;
        v1_idx = i;
      }
      if (v2Ddistance2 <= closestdistance2v) {
        closestdistance2v = v2Ddistance2;
        v2_idx = i;
      }
    }
  }

  std::vector<int> closestidx = {u1_idx, u2_idx, x1_idx,
                                 x2_idx, v1_idx, v2_idx};

  // Now store the closest X, u, v clusters for each Michel Endpoint based on
  // the above closest distance

  // Closest cluster's index will be used to only

  std::vector<double> clusx1;
  std::vector<double> clusx2;

  std::vector<double> clusu1;
  std::vector<double> clusu2;

  std::vector<double> clusv1;
  std::vector<double> clusv2;
  // Second loop over Clusters
  for (int i = 0; i < nclusters; i++) {
    // Only look at clusters found to be closest to michel in previous cluster
    // loop
    std::vector<int>::iterator it =
        std::find(closestidx.begin(), closestidx.end(), i);
    if (it == closestidx.end()) continue;

    double energy = univ.GetVecElem("cluster_energy", i);
    double time = univ.GetVecElem("cluster_time", i) / pow(10, 3);
    double pos = univ.GetVecElem("cluster_pos", i);
    double zpos = univ.GetVecElem("cluster_z", i);
    int view = univ.GetVecElem("cluster_view", i);
    double timediff = micheltime - time;
    int ismuon = univ.GetVecElem("cluster_isMuontrack", i);

    if (ismuon != 0)
      continue;  // Checking to see if Cluster is on Muon Track or not. 0 is on.
                 // 1 is not.
    if (energy < 2.) continue;
    if (timediff < 0.) continue;
    if (zpos < 5980 or zpos > 9038) continue;
    double zdiff1 = abs(zpos - michelz1);
    double zdiff2 = abs(zpos - michelz2);

    // Saving clusters with distances equal to the closest clusters. Again, this
    // is probably not the best way. I need to rewrite this section to just use
    // clusters from the index I saved in the previous loop. That way I reduce
    // the number of clusters that I have to loop over.
    Cluster current_cluster = Cluster(univ, i);
    if (view == 1) {
      // Endpoint 1

      if (i == x1_idx) {  // if current index is same as the closest endpoint 1
                          // cluster in X View
        endpoint1_clus.push_back(current_cluster);
        this->cluster_to_up_match.push_back(
            current_cluster);  // this one is redundant

        clusx1.push_back(pos);
        clusx1.push_back(zpos);
      }
      if (i == x2_idx) {  // if current index is same as the closest endpoint 2
                          // cluster in X View
        endpoint2_clus.push_back(current_cluster);
        this->cluster_to_down_match.push_back(current_cluster);  // TODO remove

        clusx2.push_back(pos);
        clusx2.push_back(zpos);
      }

    } else if (view == 2) {
      if (i == u1_idx) {  // if current index is same as the closest endpoint 1
                          // cluster in U View
        endpoint1_clus.push_back(current_cluster);
        this->cluster_to_up_match.push_back(current_cluster);  // Remove?

        clusu1.push_back(pos);
        clusu1.push_back(zpos);
      }
      if (i == u2_idx)  // if current index is same as the closest endpoint 2
                        // cluster in U View
      {
        endpoint2_clus.push_back(current_cluster);
        this->cluster_to_down_match.push_back(current_cluster);  // remove?

        clusu2.push_back(pos);
        clusu2.push_back(zpos);
      }
    } else if (view == 3) {
      if (i == v1_idx) {  // if current index is same as the closest endpoint 1
                          // cluster in V View
        endpoint1_clus.push_back(current_cluster);
        this->cluster_to_up_match.push_back(current_cluster);  // remove?

        clusv1.push_back(pos);
        clusv1.push_back(zpos);
      }
      if (i == v2_idx) {  // if current index is same as the closest endpoint 2
                          // cluster in V View
        endpoint1_clus.push_back(current_cluster);
        this->cluster_to_down_match.push_back(current_cluster);  // remove?

        clusv2.push_back(pos);
        clusv2.push_back(zpos);
      }
    }
  }  // End of loop over clusters

  // This is vector of positions for each endpoint cluster match
  std::vector<double> matchclus1;  // index [0] = x, [1] = y, [2] = z
  std::vector<double> matchclus2;  // index [0] = x, [1] = y, [2] = z

  double XZdist1 = 9999.;
  double UZdist1 = 9999.;
  double VZdist1 = 9999.;

  // Check if our cluster vectors are empty and calculate 2D distances for
  // endpoint 1
  if (!clusx1.empty()) {
    double xdif = abs(this->m_x1 - clusx1[0]);
    double zdif = abs(this->m_z1 - clusx1[1]);
    XZdist1 = sqrt(xdif * xdif + zdif * zdif);
  }
  if (!clusu1.empty()) {
    double udif = abs(this->m_u1 - clusu1[0]);
    double zdif = abs(this->m_z1 - clusu1[1]);
    UZdist1 = sqrt(udif * udif + zdif * zdif);
  }
  if (!clusv1.empty()) {
    double vdif = abs(this->m_v1 - clusv1[0]);
    double zdif = abs(this->m_z1 - clusv1[1]);
    VZdist1 = sqrt(vdif * vdif + zdif * zdif);
  }

  this->up_to_clus_XZ = XZdist1;
  this->up_to_clus_UZ = UZdist1;
  this->up_to_clus_VZ = VZdist1;
  this->up_clus_y = michely1;
  this->up_clus_x = michelx1;
  this->up_clus_z = michelz1;

  double XZdist2 = 9999.;
  double UZdist2 = 9999.;
  double VZdist2 = 9999.;

  if (!clusx2.empty()) {
    double xdif = abs(this->m_x2 - clusx2[0]);
    double zdif = abs(this->m_z2 - clusx2[1]);
    XZdist2 = sqrt(xdif * xdif + zdif * zdif);
  }
  if (!clusu2.empty()) {
    double udif = abs(this->m_u2 - clusu2[0]);
    double zdif = abs(this->m_z2 - clusu2[1]);
    UZdist2 = sqrt(udif * udif + zdif * zdif);
  }
  if (!clusv2.empty()) {
    double vdif = abs(this->m_v2 - clusv2[0]);
    double zdif = abs(this->m_z2 - clusv2[1]);
    VZdist2 = sqrt(vdif * vdif + zdif * zdif);
  }

  this->down_to_clus_XZ = XZdist2;
  this->down_to_clus_UZ = UZdist2;
  this->down_to_clus_VZ = VZdist2;
  this->down_clus_y = michely2;
  this->down_clus_x = michelx2;
  this->down_clus_z = michelz2;

  // This is the convoluted system that calculates the Endpoint
  if (XZdist1 < UZdist1 && XZdist1 < VZdist1 &&
      UZdist1 < VZdist1) {  // XU views closest
    if (!clusu1.empty() && !clusx1.empty()) {
      double yclus = (1. / sqrt(3.)) * (clusx1[0] - 2 * clusu1[0]);
      matchclus1.push_back(clusx1[0]);
      matchclus1.push_back(yclus);      // y point of match 3D point
      matchclus1.push_back(clusx1[1]);  // setting the cluster 3D point z to be
                                        // of the closest view
    }
  } else if (XZdist1 < UZdist1 && XZdist1 < VZdist1 &&
             UZdist1 > VZdist1) {  // XV closest
    if (!clusv1.empty() && !clusx1.empty()) {
      double yclus = (1. / sqrt(3.)) * (2 * clusv1[0] - clusx1[0]);
      matchclus1.push_back(clusx1[0]);
      matchclus1.push_back(yclus);
      matchclus1.push_back(clusx1[1]);  // seting the cluster 3D point z to be
                                        // of the closest view
    }
  } else if (UZdist1 < XZdist1 && UZdist1 < VZdist1 &&
             VZdist1 < XZdist1) {  // UV closest
    if (!clusv1.empty() && !clusu1.empty()) {
      double yclus = (1. / sqrt(3.)) * (clusv1[0] - clusu1[0]);
      double xclus = clusu1[0] + clusv1[0];
      matchclus1.push_back(xclus);
      matchclus1.push_back(yclus);
      matchclus1.push_back(clusu1[1]);  // seting the cluster 3D point z to be
                                        // of the closest view
    }
  } else if (UZdist1 < XZdist1 && UZdist1 < VZdist1 &&
             VZdist1 > XZdist1) {  // UX closest
    if (!clusu1.empty() && !clusx1.empty()) {
      double yclus = (1. / sqrt(3.)) * (clusx1[0] - 2 * clusu1[0]);
      matchclus1.push_back(clusx1[0]);
      matchclus1.push_back(yclus);
      matchclus1.push_back(clusu1[1]);  // seting the cluster 3D point z to be
                                        // of the closest view
      this->up_clus_y = michely1;
    }
  } else if (VZdist1 < XZdist1 && VZdist1 < UZdist1 &&
             XZdist1 < UZdist1) {  // VX closest
    if (!clusv1.empty() && !clusx1.empty()) {
      double yclus = ((1. / sqrt(3.)) * (2 * clusv1[0] - clusx1[0]));
      matchclus1.push_back(clusx1[0]);
      matchclus1.push_back(yclus);
      matchclus1.push_back(clusv1[1]);  // seting the cluster 3D point z to be
                                        // of the closest view
      this->up_clus_y = michely1;
    }
  } else if (VZdist1 < XZdist1 && VZdist1 < UZdist1 &&
             XZdist1 > UZdist1) {  // VU closest
    if (!clusu1.empty() && !clusv1.empty()) {
      double xclus = (1. / sqrt(3.)) * (clusv1[0] - clusu1[0]);
      double yclus = clusu1[0] + clusv1[0];
      matchclus1.push_back(xclus);
      matchclus1.push_back(yclus);
      matchclus1.push_back(clusv1[1]);  // seting the cluster 3D point z to be
                                        // of the closest view
      this->up_clus_y = michely1;
      this->up_clus_x = michelx1;
    }
  }
  if (XZdist2 < UZdist2 && XZdist2 < VZdist2 && UZdist2 < VZdist2) {
    if (!clusu2.empty() && !clusx2.empty()) {
      double yclus = (1. / sqrt(3.)) * (clusx2[0] - 2 * clusu2[0]);
      matchclus2.push_back(clusx2[0]);
      matchclus2.push_back(yclus);
      matchclus2.push_back(clusx2[1]);  // seting the cluster 3D point z to be
                                        // of the closest view

      this->down_clus_y = michely2;
    }
  } else if (XZdist2 < UZdist2 && XZdist2 < VZdist2 && UZdist2 > VZdist2) {
    if (!clusv2.empty() && !clusx2.empty()) {
      double yclus = (1. / sqrt(3.)) * (2 * clusv2[0] - clusx2[0]);
      matchclus2.push_back(clusx2[0]);
      matchclus2.push_back(yclus);
      matchclus2.push_back(clusx2[1]);  // seting the cluster 3D point z to be
                                        // of the closest view
      this->down_clus_y = michely2;
    }
  } else if (UZdist2 < XZdist2 && UZdist2 < VZdist2 && VZdist2 < XZdist2) {
    if (!clusu2.empty() && !clusv2.empty()) {
      double xclus = clusu2[0] + clusv2[0];
      double yclus = (1. / sqrt(3.)) * (clusv2[0] - clusu2[0]);
      matchclus2.push_back(xclus);
      matchclus2.push_back(yclus);
      matchclus2.push_back(clusu2[1]);  // seting the cluster 3D point z to be
                                        // of the closest view
      this->down_clus_y = michely2;
      this->down_clus_y = michelx2;
    }
  } else if (UZdist2 < XZdist2 && UZdist2 < VZdist2 && VZdist2 > XZdist2) {
    if (!clusu2.empty() && !clusx2.empty()) {
      double yclus = (1. / sqrt(3.)) * (clusx2[0] - 2 * clusu2[0]);
      matchclus2.push_back(clusx2[0]);
      matchclus2.push_back(yclus);
      matchclus2.push_back(clusu2[1]);  // seting the cluster 3D point z to be
                                        // of the closest view
      this->down_clus_y = michely2;
    }
  } else if (VZdist2 < XZdist2 && VZdist2 < UZdist2 && XZdist2 < UZdist2) {
    if (!clusv2.empty() && !clusx2.empty()) {
      double yclus = ((1. / sqrt(3.)) * (2 * clusv2[0] - clusx2[0]));
      matchclus2.push_back(clusx2[0]);
      matchclus2.push_back(yclus);
      matchclus2.push_back(clusv2[1]);  // seting the cluster 3D point z to be
                                        // of the closest view
      this->down_clus_y = michely2;
    }
  } else if (VZdist2 < XZdist2 && VZdist2 < UZdist2 && XZdist2 > UZdist2) {
    if (!clusu2.empty() && !clusv2.empty()) {
      double xclus = (1. / sqrt(3.)) * (clusv2[0] - clusu2[0]);
      double yclus = clusu2[0] + clusv2[0];
      matchclus2.push_back(xclus);
      matchclus2.push_back(yclus);
      matchclus2.push_back(clusv2[1]);  // seting the cluster 3D point z to be
                                        // of the closest view
      this->down_clus_y = michely2;
      this->down_clus_x = michelx2;
    }
  }

  double clusx1diff = 9999.;
  double clusy1diff = 9999.;
  double clusz1diff = 9999.;

  double mclusx1diff = 9999.;
  double mclusy1diff = 9999.;
  double mclusz1diff = 9999.;

  double michvtx_x1diff = 9999.;
  double michvtx_x2diff = 9999.;
  double michvtx_y1diff = 9999.;
  double michvtx_y2diff = 9999.;
  double michvtx_z1diff = 9999.;
  double michvtx_z2diff = 9999.;

  if (!matchclus1.empty()) {
    clusx1diff = vtx_x - matchclus1[0];
    clusy1diff = vtx_y - matchclus1[1];
    clusz1diff = vtx_z - matchclus1[2];
    mclusx1diff = michelx1 - matchclus1[0];
    mclusy1diff = michely1 - matchclus1[1];
    mclusz1diff = michelz1 - matchclus1[2];
    michvtx_x1diff = michelx1 - vtx_x;
    michvtx_y1diff = michely1 - vtx_y;
    michvtx_z1diff = michelz1 - vtx_z;
  }
  double clusx2diff = 9999.;
  double clusy2diff = 9999.;
  double clusz2diff = 9999.;

  double mclusx2diff = 9999.;
  double mclusy2diff = 9999.;
  double mclusz2diff = 9999.;

  if (!matchclus2.empty()) {
    clusx2diff = vtx_x - matchclus2[0];
    clusy2diff = vtx_y - matchclus2[1];
    clusz2diff = vtx_z - matchclus2[2];
    mclusx2diff = michelx2 - matchclus2[0];
    mclusy2diff = michely2 - matchclus2[1];
    mclusz2diff = michelz2 - matchclus2[2];
    michvtx_x2diff = michelx2 - vtx_x;
    michvtx_y2diff = michely2 - vtx_y;
    michvtx_z2diff = michelz2 - vtx_z;
  }
  /// 2 types of 3D distance calculations for Cluster matching:
  // 1. Cluster to Vertex
  double dist1 =
      sqrt(pow(clusx1diff, 2) + pow(clusy1diff, 2) + pow(clusz1diff, 2));
  double dist2 =
      sqrt(pow(clusx2diff, 2) + pow(clusy2diff, 2) + pow(clusz2diff, 2));
  // 2. Cluster to Michel
  double mdist1 =
      sqrt(pow(mclusx1diff, 2) + pow(mclusy1diff, 2) + pow(mclusz1diff, 2));
  double mdist2 =
      sqrt(pow(mclusx2diff, 2) + pow(mclusy2diff, 2) + pow(mclusz2diff, 2));
  this->nnonmuclusters = nonmuclus;
  this->down_clus_michel_dist3D = mdist2;
  this->up_clus_michel_dist3D = mdist1;
  this->up_to_cluster_dist3D = dist1;
  this->down_to_cluster_dist3D = dist2;
  double michdist1 = sqrt(pow(michvtx_x1diff, 2) + pow(michvtx_y1diff, 2) +
                          pow(michvtx_z1diff, 2));
  double michdist2 = sqrt(pow(michvtx_x2diff, 2) + pow(michvtx_y2diff, 2) +
                          pow(michvtx_z2diff, 2));
  this->up_clus_michvtx_dist3D = michdist1;
  this->down_clus_michvtx_dist3D = michdist2;
  // Marks which endpoint is the closest match for this match type
  if (dist1 < dist2)
    this->clus_endpoint = 1;
  else if (dist1 > dist2)
    this->clus_endpoint = 2;
}

}  // namespace LowRecoilPion

#endif
