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
#include "weightRPA.h"
#include "weight_2p2h.h"
#include "weightLowQ2Pi.h"
#include "weightDIS.h"
#include "weightZExp.h"

//PlotUtils systematic universes classes (new sys framework)
#include "MinervaUniverse.h"
#include "GenieSystematics.h"
#include "GeantHadronSystematics.h"
#include "MnvTuneSystematics.h"
#include "MinosEfficiencySystematics.h"
#include "MuonSystematics.h"
#include "MuonResolutionSystematics.h"
#include "MichelSystematics.h"
#include "AngleSystematics.h"
#include "ResponseSystematics.h"

//Utils
#include "TargetUtils.h"
#include "POTCounter.h"


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

