#ifndef LOWRECOILPIONFUNCTIONS_H
#define LOWRECOILPIONFUNCTIONS_H

virtual int GetNMichels() const {
  return GetInt("FittedMichel_michel_fitPass_sz");
}

ROOT::Math::XYZTVector GetVertex() const {
  ROOT::Math::XYZTVector result;
  result.SetCoordinates(GetVec<double>("vtx").data());
  return result;
}

// MeV
virtual double NewRecoilE() const {
  const double M_pi = 139.57039;  // in MeV
  double recoil = NewEavail() + M_pi;
  double newrecoil = recoil;
  return newrecoil;
}

// GeV^2
virtual double GetQ2Reco() const {
  const double M_mu = 105.6583 / 1000.;  // Converting to GeV
  double enu = GetEmu() / 1000. + (NewRecoilE()) / 1000.;
  double pmucos = (GetPmu() / 1000.) * cos(GetThetamu());
  double q2reco = 2. * enu * (GetEmu() / 1000. - pmucos) - (M_mu * M_mu);
  return q2reco;
}

// GeV
virtual double Getq3_lowrecoil() const {
  double eavail = NewRecoilE() / 1000.;
  double q2 = GetQ2Reco();
  double q3mec = sqrt(eavail * eavail + q2);
  return q3mec;  // Using Hang's q3 definition TODO: check to see if using
                 // Hang's definition is ok. //q3mec;
}

// GeV/c
double GetMuonPT() const { return GetPmu() / 1000. * sin(GetThetamu()); }

// MeV
virtual std::vector<double> GetTrackerECALMuFuzz() const {
  double trk_mufuzz = 0.0;
  double ecal_mufuzz = 0.0;
  int nfuzz = GetInt("muon_fuzz_per_plane_r80_planeIDs_sz");
  if (nfuzz == 0) return {0.0, 0.0};
  for (int i = 0; i < nfuzz; i++) {
    int planeID = GetVecElem("muon_fuzz_per_plane_r80_planeIDs", i);
    if (planeID < 1504968704 || planeID > 1709703168) continue;
    double fuzze = GetVecElem("muon_fuzz_per_plane_r80_energies", i);
    if (planeID > 1504968704 and planeID < 1560805376)
      trk_mufuzz += fuzze;
    else if (planeID > 1700003840 and planeID < 1709703168)
      ecal_mufuzz += fuzze;
  }

  return {trk_mufuzz, ecal_mufuzz};
}

// MeV
virtual double NewEavail() const {
  double recoiltracker =
      GetDouble("blob_recoil_E_tracker") - GetTrackerECALMuFuzz()[0];
  double recoilEcal =
      GetDouble("blob_recoil_E_ecal") - GetTrackerECALMuFuzz()[1];
  const double Eavailable_scale = 1.17;
  double eavail = recoiltracker + recoilEcal;
  return eavail * Eavailable_scale;
}

// MeV
virtual double GetTpiFromRange(double range) const {
//  double tpiest = 0.2142 * range + 2.864 * sqrt(range);
  double tpiest = 0.128706 * range + 3.42486 * sqrt(range);
  return tpiest;  // Tpi in MeV
}

int GetFittedMichelsOnly() const {
  std::vector<int> nfitted = GetVecInt("FittedMichel_michel_fitPass");
  int count = std::count(nfitted.begin(), nfitted.end(), 1);
  return count;
}

double GetNClusters() const {
  int nclusters = GetInt("cluster_view_sz");
  int nonmuclus = 0;
  double nclus = 1.0 * nclusters;
  return nclus;  // nonmuclus;
}

#endif  // LOWRECOILPIONFUNCTIONS_H
