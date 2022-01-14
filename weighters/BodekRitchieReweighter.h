//File: BodekRitchieReweighter.h
//Brief: Reweights MINERvA's Central Value Monte Carlo simulation to match
//       the Bodek-Ritchie prediction for a kinematic tail.  Used in Marvin's
//       low recoil analysis.
//Author: Andrew Olivier aolivier@ur.rochester.edu

#ifndef BODEKRITCHIEREWEIGHTER_H
#define BODEKRITCHIEREWEIGHTER_H

//MAT-MINERvA includes
#include "weighters/Reweighter.h"
#include "weighters/weightGenieBodekRitchieClass.h"

namespace PlotUtils
{
  //Map Rik's weighter interface to the MAT
  template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
  class BodekRitchieReweighter: public Reweighter<UNIVERSE, EVENT>
  {
    public:
      BodekRitchieReweighter(int mode): m_weighter(), m_mode(mode) {}
      virtual ~BodekRitchieReweighter() = default;

      virtual double GetWeight(const UNIVERSE& univ, const EVENT& event) const
      {
        return m_weighter.getWeight(m_mode, univ.GetInt("mc_er_nPart"), univ.GetInt("mc_intType"),
                                    univ.GetInt("mc_targetA"), univ.GetVecInt("mc_er_status"),
                                    univ.GetVecInt("mc_er_ID"), univ.GetVecDouble("mc_er_Px"),
                                    univ.GetVecDouble("mc_er_Py"), univ.GetVecDouble("mc_er_Pz"));
      }

      virtual std::string GetName() const { return "BodekRitchieTail"; }

      virtual bool DependsReco() const { return false; }

    private:
      weightGenieBodekRitchieClass m_weighter;
      int m_mode; //Options are 1 and 2 so far.  Mode 2 talks about weighting some tail down
                  //in addition to whatever else weightGenieBodekRitchieClass does.
  };
}

#endif //BODEKRITCHIEREWEIGHTER_H
