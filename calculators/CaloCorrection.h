//File: CaloCorrection.h
//Brief: A CaloCorrection corrects recoil energy (without vertex box, with OD, etc.) to
//       a better estimator for total hadronic energy.  It mimics Minerva::CalorimetryUtils
//       in doing this by interpolating a multiplicative constant between a series of points.
//Author: Andrew Olivier aolivier@ur.rochester.edu


//c++ includes
#include <unordered_map>
#include <vector>

namespace util
{
  class CaloCorrection
  {
    public:
      CaloCorrection(const std::string& caloFile, const std::string& tuningName = "Default");

      //Load one CaloCorrection from a plaintext file with the following format that originated
      //from https://nusoft.fnal.gov/minerva/minervadat/software_doxygen/HEAD/MINERVA/classCalorimetryUtils.html#261bfea9ea2c8880fbfa352fdfb9f60f:
      //IMPORT <path to file>: Path to another file to load.  Maximum recursion depth of 5.
      //                       May include environment variables. 
      //BEGIN <tuning name>: Begins a tuning made up of 0 or more POLYPOINTS.  All tunings except
      //                     tuningName are ignored.
      //POLYPOINT <x value> <y value>: A knot on a calorimetric polyline correction.  Between
      //                               <x value> input energy and the next point, multiply energy
      //                               by <y value>.
      //END: Ends a block of POLYPOINTS that started with BEGIN
      //SCALE, CALCONST: Ignored because they were applied in the Gaudi stage if at all.
      //# denotes a comment
      static std::unordered_map<std::string, CaloCorrection> parse(const std::string& caloFile);

      //Apply correction
      double eCorrection(const double rawRecoil) const;

    private:
      CaloCorrection();

      struct PolyPoint
      {
        double threshold;
        double correction;
      };

      double fScale; //Overall energy scale applied to all Clusters

      std::vector<PolyPoint> fPoints;
  };
}