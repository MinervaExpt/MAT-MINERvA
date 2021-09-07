//==============================================================================
/*! @brief Minerva Low recoil analyses functions Include
  this file inside of your User::CVUniverse class definition.

  Public-use functions here defined:

  - double GetEAvailable() const (MeV)
  - double GetCalorimetryQ0() const (MeV)

*/
//==============================================================================
#ifndef LOWRECOILFUNCTIONS_H
#define LOWRECOILFUNCTIONS_H

enum Subdet{nucl,tracker,ecal,hcal,od};

double GetEAvailable() const {
  const double Eavailable_scale = 1.17;
  std::vector<Subdet> eavail_subdets;
  eavail_subdets.push_back(tracker);
  eavail_subdets.push_back(ecal);
  double eavail = 0;
  for (std::vector<Subdet>::const_iterator it = eavail_subdets.begin();it!=eavail_subdets.end();++it) {
    eavail += std::max(0.0,GetSubdetE(*it)-GetSubdetFuzz(*it));
  }
  return eavail*Eavailable_scale;
}

double GetCalorimetryQ0() const {
  return GetDouble("MasterAnaDev_recoil_E");
}

double GetSubdetE(Subdet subdet) const {
  std::string eavail_prefix = "blob_recoil_E";
  return std::max(0.0,GetDouble(GetSubdetString(subdet,eavail_prefix).c_str()));
}

std::string GetSubdetString(Subdet subdet,const std::string& prefix) const {
  switch (subdet) {
  case nucl:
    return prefix+"_nucl";
  case tracker:
    return prefix+"_tracker";
  case ecal:
    return prefix+"_ecal";
  case hcal:
    return prefix+"_hcal";
  case od:
    return prefix+"_od";
  }
  return "";
}

double GetSubdetFuzz(Subdet subdet) const {
  return 0; // for now.
}

double GetLowRecoilQ3() const {
  double q0 = GetCalorimetryQ0();
  double E_lep = GetEmu();
  double E_nu = q0+E_lep;
  double q2 = calcqsquare_calorimetry(E_lep,GetThetamu(),GetPmu(),q0);
  return calcq3(q2,E_nu,E_lep);
}

double calcqsquare_calorimetry(double E_lep, double theta_lep, double p_lep, double q0) const {
  if (E_lep <= p_lep) throw std::runtime_error("Calculate Q^2 error: lepton mass is 0(imaginary), because Energy is equal to (less than) momentum.");
  double mass_square = E_lep*E_lep - p_lep *p_lep;
  double Enu = E_lep+q0;
  double q2 = 2*Enu *(E_lep - p_lep * std::cos(theta_lep))-mass_square;
  return std::max(0.,q2);
}

#endif /* LOWRECOILFUNCTIONS_H */
