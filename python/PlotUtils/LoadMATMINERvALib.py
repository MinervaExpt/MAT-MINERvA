"""
  LoadPlotUtilsLib.py:
   The code necessary to load the libplotutils.so library
   so that PlotUtils C++ objects are available.
  
   Original author: J. Wolcott (jwolcott@fnal.gov)
                    November 2012
"""

# hms 2021-11-20 comment out classes that moved to MAT-MINERvA
# to use this you need to
# export PYTHONPATH=$WHEREMATIS/MAT/python:$WHEREMATIS/MAT/python/PlotUtils

import sys
import os

import ROOT

# List of classes to copy from the Reflex library into the "PlotUtils" namespace
CLASSES_TO_LOAD = [
# from MAT
    "MnvH1D",
    "MnvH2D",
    "MnvH3D",
    "MnvHist",
    "MnvHistoConstrainer",
    "MnvLatErrorBand",
    "MnvLatErrorBand2D",
    "MnvLatErrorBand3D",
    "MnvVertErrorBand",
    "MnvVertErrorBand2D",
    "MnvVertErrorBand3D",
    "MnvPlotter",
    "TargetUtils",
    "MnvNormalizer",
    "FluxReweighter",
    "ChainWrapper",
    "GridCanvas",
    "weight_2p2h",
    "weightLowQ2Pi",
    "weightRPA",
    "weightDIS",
    "weightZExp",
    "MinervaUniverse",
    "FluxUniverse",
    "GenieUniverse",
    "GenieRvx1piUniverse",
    "GenieFaCCQEUniverse",
    "MuonUniverseMinerva",
    "MuonUniverseMinos",
    "MuonResolutionUniverse",
    "BeamAngleXUniverse",
    "BeamAngleYUniverse",
    "Universe2p2h",
    "RPAUniverse",
    "ResponseUniverse",
    "LowQ2PionUniverse",
    "MinosEfficiencyUniverse",
    "flux_reweighter",
    #From MAT-MINERvA
    # "RecoilEnergyFunctions",
    #"TruthFunctions",
    #"WeightFunctions",
    #"AngleSystematics",
    #"CCQE3DFitsSystematics",
    #"DiffractiveUncSystematics",
    #"GeantHadronSystematics",
    #"GenericVerticalSystematic",
    #"GenieSystematics",
    #"MichelSystematics",
    "MinervaUniverse",
    #"MinosEfficiencySystematics",
    #"MnvTuneSystematics",
    #"MuonResolutionSystematics",
    #"MuonSystematics",
    #"RecoProtonKECutSystematics",
    #"ResponseSystematics",
    #"TargetMassSystematics",
    #"TrueProtonKECutSystematics",
    #"AnaBinning",
    #"HyperDimLinearizer",
    #"MnvNormalization",
    #"NSFDefaults",
    "POTCounter",
    #"ParticleResponseDefaults",
    #"PhysicsVariables",
    #"PlotUtilsPhysicalConstants",
    "TargetUtils",
    #"AMUDISReweighter",
    #"FSIReweighter",
    #"FluxAndCVReweighter",
    #"GENIEReweighter",
    #"GeantNeutronCVReweighter",
    #"LowQ2PiReweighter",
    #"LowRecoil2p2hReweighter",
    #"MINOSEfficiencyReweighter",
    #"MKReweighter",
    #"MinosMuonEfficiencyCorrection",
    #"MnvHadronReweight",
    #"Model",
    #"NeutronInelasticReweighter",
    #"RPAReweighter",
    #"Reweighter",
    #"SuSAFromValencia2p2hReweighter",
    #"weightCoherentPi",
    #"weightDIS",
    #"weightGenieBodekRitchieClass",
    #"weightLowQ2Pi",
    #"weightMK",
    #"weightNuclearScreening",
    #"weightRPA",
    #"weightRemoveUnphysical2p2hExtendedEventsClass",
    #"weightSusaGenieQEClass",
    #"weightSusaValenciaClass",
    #"weightZExp",
    #"weight_2p2h",
    #"weight_fsi",
    #"weight_fsi_absorption",
    #"weight_fsi_cai",
    #"weight_minervaq2qe"
]

# new code for ROOT6
if os.getenv("PLOTUTILSVERSION")=="ROOT6":
# replace with voodoo found at https://gitlab.cern.ch/clemenci/Gaudi/commit/47772dfbb429a1392e12143403314e0abb45322d
    try:
        import PyCintex # needed to enable Cintex if it exists
        del PyCintex
    except ImportError:
        pass
else:
    import PyCintex

import PlotUtils

# the ROOT interpreter gobbles up any arguments
# that might be in sys.argv, whether you like it or not.
# we need to tuck them away for later.
args = sys.argv[:]
sys.argv = sys.argv[:0]

if "PLOTUTILSROOT" in os.environ:
  #ROOT.gSystem.Load("libMAT")
  ROOT.gSystem.Load("libMAT-MINERvA")

	# copy the classes from the Reflex library
	# into the "PlotUtils" namespace for more
	# straightforward access.
  for cls in CLASSES_TO_LOAD:
    setattr(PlotUtils, cls, getattr(ROOT.PlotUtils, cls))
else:
	print >> sys.stderr, "Note: $PLOTUTILSROOT is not defined in the current environment.  PlotUtils libraries were not loaded."
	
# now that ROOT has done its thing,
# we can restore the arguments...
sys.argv = args
