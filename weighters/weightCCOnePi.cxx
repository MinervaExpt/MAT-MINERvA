#include "weighters/weightCCOnePi.h"

using namespace PlotUtils;

weightCCOnePi::weightCCOnePi()
    : weight_file(nullptr),
      __h2d_Weight(NULL),
      __tpi_max(350),
      __q2_max(3.e6), 
      __tpi_min(0),
      __q2_min(0) {
  char* mparalocation = std::getenv("TOPDIR");
  std::string f = std::string(mparalocation) +
                  "/MAT-MINERvA/universes/RatioOutput_mixtpi_vs_q2_KevinMod.root";
  read(f);
}

void weightCCOnePi::read(const std::string f) {
  weight_file = TFile::Open(f.c_str(), "READONLY");
  assert(weight_file);

  __h2d_Weight = static_cast<TH2D*>(weight_file->Get("Ratio_mixtpi_vs_q2"));

}

weightCCOnePi::~weightCCOnePi() {
  delete __h2d_Weight;
}


// coherent pion weights
double weightCCOnePi::get_weight(double tpi_true, double q2_true) {
  if (tpi_true < __tpi_min || tpi_true > __tpi_max) return 1.0;
  if (q2_true < __q2_min || q2_true > __q2_max) return 1.0;

  return __h2d_Weight->GetBinContent(__h2d_Weight->FindBin(tpi_true,q2_true));
}


PlotUtils::weightCCOnePi& PlotUtils::weight_CCOnePi() {
  static PlotUtils::weightCCOnePi* _weight_CCOnePi =
      new PlotUtils::weightCCOnePi();
  return *_weight_CCOnePi;
}
