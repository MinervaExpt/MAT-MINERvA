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
struct Q0Spline {
  std::vector<double> poly_x;
  std::vector<double> poly_y;
  double scale;
  bool set=false;
}

double GetEAvailable(double scale = 1.17) const {
  std::vector<Subdet> eavail_subdets = {tracker,ecal};
  double eavail = GetSubdetSum(eavail_subdets,true);
  return eavail*scale;
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

double GetSubdetSum(const std::vector<Subdet>& subdets, bool fuzz = false) {
  for (std::vector<Subdet>::const_iterator it = eavail_subdets.begin();it!=eavail_subdets.end();++it) {
    double e = GetSubdetE(*it) - (fuzz? GetSubdetFuzz(*it) : 0);
    eavail += std::max(0.0,e);
  }
}

double GetSplinedQ0(bool fuzz = true) const {
  std::vector<Subdet> q0_subdets = {tracker,ecal,hcal,od,nucl};
  double p_spline = GetSubdetSum(q0_subdets,fuzz);
  for (std::vector<Subdet>::const_iterator it = eavail_subdets.begin();it!=eavail_subdets.end();++it) {
    eavail += std::max(0.0,GetSubdetE(*it));
  }
  return Spline(eavail);
}

double Spline(double rawEnergy) const {
  if (!current_spline.set) LoadNukeCC_Nu_Tracker();
  rawEnergy*= 1e-3;
  const std::vector<double>& poly_x = current_spline.poly_x;
  const std::vector<double>& poly_y = current_spline.poly_y;
  size_t poly_n = poly_x.size()
  for(size_t i = 0; i != poly_n-1; ++i) {
    if(rawEnergy < poly_x[i+1]) {
      // interpolate between nodes
      rawEnergy = poly_y[i]+(poly_y[i+1]-poly_y[i])*(rawEnergy-poly_x[i])/(poly_x[i+1]-poly_x[i]);
      return (rawEnergy < 0.0)? 0.0 :current_spline.scale * rawEnergy * 1e3; // Change back to MeV
    }
  }
}

bool LoadNukeCC_Nu_Tracker() {
  double scale = 1.51324
  std::vector<std::pair<double,double>> points;
  points.push_back(0,0);
  points.push_back(1.05282,1.0433);
  points.push_back(1.41762,1.43065);
  points.push_back(1.89408,1.9389);
  points.push_back(2.38637,2.43997);
  points.push_back(2.88095,2.93705);
  points.push_back(3.37431,3.43602);
  points.push_back(5.07854,5.15973);
  points.push_back(10.6784,10.7691);
  points.push_back(15.931,15.9287);
  points.push_back(25.6898,25.6322);
  points.push_back(50,50);
  std::vector<double> x;
  std::vector<double> y;

  for (auto point : points) {
    x.push_back(point[0]);
    y.push_back(point[1]);
  }
  SetupSpline(x,y,scale);
}

bool SetupSpline(const std::vector<double>& x, const std::vector<double>& y,double scale) {
  if (!current_spline.set) return false;
  if (x.size() != y.size() || x.size() < 2) return false;
  //check whether x is monotonically increasing
  double p0 = x[0];
  for (size_t  i = 1 ; i< x.size();++i) {
    if (x[i]<p0) return false;
    p0=x[i];
  }
  current_spline.poly_x = x;
  current_spline.poly_y = y;
  current_spline.scale = scale;
  current_spline.set = true;
  return true;
}

Q0Spline current_spline;



#endif /* LOWRECOILFUNCTIONS_H */
