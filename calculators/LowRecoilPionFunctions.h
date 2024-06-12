#ifndef LOWRECOILPIONFUNCTIONS_H
#define LOWRECOILPIONFUNCTIONS_H

#include "PlotUtils/PlotUtilsPhysicalConstants.h"

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
  double recoil = NewEavail() + MinervaUnits::M_pion;
  double newrecoil = recoil;
  return newrecoil;
}

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
  const double Eavailable_scale = 1.17; // determined from a study by Phil
  double eavail = recoiltracker + recoilEcal;
  return eavail * Eavailable_scale;
}

// mm --> MeV
virtual double GetTpiFromRange(double range) const {
  return NSFDefaults::tpi_from_michel_range_fit_p0_cv * range + 
      NSFDefaults::tpi_from_michel_range_fit_p1_cv * sqrt(range);
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
