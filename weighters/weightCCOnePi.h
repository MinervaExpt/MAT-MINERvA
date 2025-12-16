#ifndef weightCCOnePi_h
#define weightCCOnePi_h


#include <TFile.h>
#include <TH2D.h>
#include <TString.h>

#include <cassert>
#include <iostream>
#include <stdexcept>

namespace PlotUtils {

class weightCCOnePi {
 public:
  weightCCOnePi();

  double get_weight(double tpi_true, double q2_true);

  ~weightCCOnePi();

 private:
  TFile* weight_file;

  TH2D* __h2d_Weight;

  double __tpi_max;
  double __q2_max;
  double __tpi_min;
  double __q2_min;

  void read(const std::string f);
};

PlotUtils::weightCCOnePi& weight_CCOnePi();
}  // namespace PlotUtils

#endif
