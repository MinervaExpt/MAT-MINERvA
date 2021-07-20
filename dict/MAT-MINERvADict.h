#ifndef MATMINERVADICT_H 
#define MATMINERVADICT_H 1

// Include files for PlotUtils dictionary.

/** @file PlotUtilsDict.h
 *  
 *
 *  @author Jeremy Wolcott <jwolcott@fnal.gov>
 *  @date   2012-11-25
 */
// ============================================================================
// PlotUtils
// ============================================================================

// here we need to include all the header files
// for the classes we want to make dictionaries of

#include <vector>

// PlotUtils weight classes
#include "weighters/weightRPA.h"
#include "weight_2p2h.h"
#include "weighters/weightLowQ2Pi.h"
#include "weighters/weightDIS.h"
#include "weighters/weightZExp.h"

//PlotUtils systematic universes classes (new sys framework)
#include "universes/MinervaUniverse.h"
#include "universes/GenieSystematics.h"
#include "universes/GeantHadronSystematics.h"
#include "universes/MnvTuneSystematics.h"
#include "universes/MinosEfficiencySystematics.h"
#include "universes/MuonSystematics.h"
#include "universes/MuonResolutionSystematics.h"
#include "universes/MichelSystematics.h"
#include "universes/AngleSystematics.h"
#include "universes/ResponseSystematics.h"


//TODO: Do I need this?
//#include "ErrorHandler.h"

// this garbage is necessary so that gccxml is able to create dictionaries for these custom containers
// (since it otherwise doesn't know which specific version of these templated classes to instantiate)
// see: http://root.cern.ch/root/roottalk/roottalk10/0035.html
// somehow std::map<>s seem to be instantiated somewhere else, so explicit instantiation is not necessary?
#ifdef __GCCXML__

// Use extern keyword because these functions are instantiated in FluxReweighter.cxx already

#endif
#endif // MATMINERVADICT_H

