# Almost all CMake files should start with this
# You should always specify a range with the newest
# and oldest tested versions of CMake. This will ensure
# you pick up the best policies.

cmake_minimum_required(VERSION 3.1...3.14)

# This is your project statement. You should always list languages;
# Listing the version is nice here since it sets lots of useful variables

project(PlotUtils VERSION 1.0 LANGUAGES CXX)

# If you set any CMAKE_ variables, that can go here.
# (But usually don't do this, except maybe for C++ standard)

# You need to tell CMake where to find the ROOT installation. This can be done in a number of ways:
#   - ROOT built with classic configure/make use the provided $ROOTSYS/etc/cmake/FindROOT.cmake
#   - ROOT built with CMake. Add in CMAKE_PREFIX_PATH the installation prefix for ROOT
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS} $ENV{BOOSTDIR})

set(CMAKE_CXX_FLAGS "-Wno-format -Wno-unused-value" )  # get rid of errors from Form and assert

#---Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
find_package(ROOT REQUIRED COMPONENTS RIO Net)

#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
include(${ROOT_USE_FILE})
link_directories(${ROOT_LIBRARY_DIR})

MESSAGE("ROOT_LIBRARY_DIR:${ROOT_LIBRARY_DIR}")

# put some special flags here forgot what the DBUILD_SHARED_LIBS=OFF does..

add_definitions( -DFORM -DMNVROOT6 -DPLOTUTILS_STANDALONE "-DBUILD_SHARED_LIBS=OFF")

# do this or mac will make it dylib
set(CMAKE_SHARED_LIBRARY_SUFFIX ".so")


# tell it about the source files
set(SRC_FILES
PlotUtils/BaseUniverse.cxx
PlotUtils/ChainWrapper.cxx
PlotUtils/Cutter.cxx
PlotUtils/ErrorHandler.cxx
PlotUtils/FluxReweighter.cxx
PlotUtils/FluxSystematics.cxx
PlotUtils/GridCanvas.cxx
PlotUtils/Hist2DWrapper.cxx
PlotUtils/HistFolio.cxx
PlotUtils/HistWrapper.cxx
PlotUtils/HistogramUtils.cxx
PlotUtils/MacroUtil.cxx
PlotUtils/MnvApplication.cxx
PlotUtils/MnvColors.cxx
PlotUtils/MnvFluxConstraint.cxx
PlotUtils/MnvH1D.cxx
PlotUtils/MnvH1DToCSV.cxx
PlotUtils/MnvH2D.cxx
PlotUtils/MnvH2DLog.cxx
PlotUtils/MnvH2DToCSV.cxx
PlotUtils/MnvH3D.cxx
PlotUtils/MnvLatErrorBand.cxx
PlotUtils/MnvLatErrorBand2D.cxx
PlotUtils/MnvLatErrorBand3D.cxx
PlotUtils/MnvNuclearModelWeight.cxx
PlotUtils/MnvPlotter.cxx
PlotUtils/MnvVertErrorBand.cxx
PlotUtils/MnvVertErrorBand2D.cxx
PlotUtils/MnvVertErrorBand3D.cxx
PlotUtils/NuclModUtils.cxx
PlotUtils/ROOTglob.cxx
PlotUtils/Table.cxx
PlotUtils/TreeWrapper.cxx
PlotUtils/Variable2DBase.cxx
PlotUtils/VariableBase.cxx
PlotUtils/genie_particle.cxx
PlotUtils/makeChainWrapper.cxx
universes/AngleSystematics.cxx
universes/CCQE3DFitsSystematics.cxx
universes/DiffractiveUncSystematics.cxx
universes/GeantHadronSystematics.cxx
universes/GenieSystematics.cxx
universes/MichelSystematics.cxx
universes/MinervaUniverse.cxx
universes/MinosEfficiencySystematics.cxx
universes/MnvTuneSystematics.cxx
universes/MuonResolutionSystematics.cxx
universes/MuonSystematics.cxx
universes/RecoProtonKECutSystematics.cxx
universes/ResponseSystematics.cxx
universes/TargetMassSystematics.cxx
universes/TrueProtonKECutSystematics.cxx
utilities/AnaBinning.cxx
utilities/HyperDimLinearizer.cxx
utilities/MnvNormalization.cxx
utilities/POTCounter.cxx
utilities/ParticleResponseDefaults.cxx
utilities/PhysicsVariables.cxx
utilities/TargetUtils.cxx
weighters/MinosMuonEfficiencyCorrection.cxx
weighters/MnvHadronReweight.cxx
weighters/weightCoherentPi.cxx
weighters/weightDIS.cxx
weighters/weightGenieBodekRitchieClass.cxx
weighters/weightLowQ2Pi.cxx
weighters/weightMK.cxx
weighters/weightNuclearScreening.cxx
weighters/weightRPA.cxx
weighters/weightRemoveUnphysical2p2hExtendedEventsClass.cxx
weighters/weightSusaGenieQEClass.cxx
weighters/weightSusaValenciaClass.cxx
weighters/weightZExp.cxx
weighters/weight_2p2h.cxx
weighters/weight_fsi.cxx
weighters/weight_fsi_absorption.cxx
weighters/weight_fsi_cai.cxx
weighters/weight_minervaq2qe.cxx
)

# tell it about the header files

set(H_FILES
PlotUtils/BaseUniverse.h
PlotUtils/ChainWrapper.h
PlotUtils/CrashOnROOTMessage.h
PlotUtils/Cut.h
PlotUtils/Cutter.h
PlotUtils/ErrorHandler.h
PlotUtils/Exceptions.h
PlotUtils/FluxReweighter.h
PlotUtils/FluxSystematics.h
PlotUtils/GridCanvas.h
PlotUtils/Hist2DWrapper.h
PlotUtils/HistFolio.h
PlotUtils/HistWrapper.h
PlotUtils/HistogramUtils.h
PlotUtils/MacroUtil.h
PlotUtils/MnvApplication.h
PlotUtils/MnvColors.h
PlotUtils/MnvFluxConstraint.h
PlotUtils/MnvH1D.h
PlotUtils/MnvH1DToCSV.h
PlotUtils/MnvH2D.h
PlotUtils/MnvH2DLog.h
PlotUtils/MnvH2DToCSV.h
PlotUtils/MnvH3D.h
PlotUtils/MnvLatErrorBand.h
PlotUtils/MnvLatErrorBand2D.h
PlotUtils/MnvLatErrorBand3D.h
PlotUtils/MnvNuclearModelWeight.h
PlotUtils/MnvPlotter.h
PlotUtils/MnvVertErrorBand.h
PlotUtils/MnvVertErrorBand2D.h
PlotUtils/MnvVertErrorBand3D.h
PlotUtils/NamedCategory.h
PlotUtils/NuclModUtils.h
PlotUtils/ROOTglob.h
PlotUtils/Table.h
PlotUtils/TreeWrapper.h
PlotUtils/Variable2DBase.h
PlotUtils/VariableBase.h
PlotUtils/genie_particle.h
PlotUtils/makeChainWrapper.h
universes/AngleSystematics.h
universes/CCQE3DFitsSystematics.h
universes/DiffractiveUncSystematics.h
universes/GeantHadronSystematics.h
universes/GenieSystematics.h
universes/MichelSystematics.h
universes/MinervaUniverse.h
universes/MinosEfficiencySystematics.h
universes/MnvTuneSystematics.h
universes/MuonResolutionSystematics.h
universes/MuonSystematics.h
universes/RecoProtonKECutSystematics.h
universes/ResponseSystematics.h
universes/TargetMassSystematics.h
universes/TrueProtonKECutSystematics.h
utilities/AnaBinning.h
utilities/HyperDimLinearizer.h
utilities/MnvNormalization.h
utilities/NSFDefaults.h
utilities/POTCounter.h
utilities/ParticleResponseDefaults.h
utilities/PhysicsVariables.h
utilities/PlotUtilsPhysicalConstants.h
utilities/TargetUtils.h
weighters/AMUDISReweighter.h
weighters/FSIReweighter.h
weighters/FluxAndCVReweighter.h
weighters/GENIEReweighter.h
weighters/GeantNeutronCVReweighter.h
weighters/LowQ2PiReweighter.h
weighters/LowRecoil2p2hReweighter.h
weighters/MINOSEfficiencyReweighter.h
weighters/MKReweighter.h
weighters/MinosMuonEfficiencyCorrection.h
weighters/MnvHadronReweight.h
weighters/Model.h
weighters/RPAReweighter.h
weighters/Reweighter.h
weighters/SuSAFromValencia2p2hReweighter.h
weighters/weightCoherentPi.h
weighters/weightDIS.h
weighters/weightGenieBodekRitchieClass.h
weighters/weightLowQ2Pi.h
weighters/weightMK.h
weighters/weightNuclearScreening.h
weighters/weightRPA.h
weighters/weightRemoveUnphysical2p2hExtendedEventsClass.h
weighters/weightSusaGenieQEClass.h
weighters/weightSusaValenciaClass.h
weighters/weightZExp.h
weighters/weight_2p2h.h
weighters/weight_fsi.h
weighters/weight_fsi_absorption.h
weighters/weight_fsi_cai.h
weighters/weight_minervaq2qe.h
)

# need to know about Boost on macs so check that you know where it is.

MESSAGE("boost:$ENV{BOOSTDIR}")
include_directories( $ENV{PLOTUTILSROOT} $ENV{PLOTUTILSROOT}/include  $ENV{PLOTUTILSROOT}/include/PlotUtils # $ENV{PLOTUTILSMINERVAROOT} $ENV{PLOTUTILSMINERVAROOT}/utilities #$ENV{PLOTUTILSMINERVAROOT}/universes $ENV{PLOTUTILSMINERVAROOT}/weighters
)

# now generate the dictionary
MESSAGE (${CMAKE_SOURCE_DIR})
ROOT_GENERATE_DICTIONARY(G__${PROJECT_NAME} ${H_FILES} LINKDEF ${CMAKE_SOURCE_DIR}/../dict/Linkdef.h)

# and then the shared library
# Output libname matches target name, with the usual extensions on your system

add_library(${PROJECT_NAME} SHARED ${SRC_FILES} G__${PROJECT_NAME}.cxx )
#add_library(${PROJECT_NAME} SHARED ${SRC_FILES} ${PROJECT_NAME}.cxx )

MESSAGE("ROOT_LIBRARIES:${ROOT_LIBRARIES}")

MESSAGE("PROJECT:${PROJECT_NAME}")

#target_link_libraries(${PROJECT_NAME} ${ROOT_LIBRARIES} MathMore)
target_link_libraries(${PROJECT_NAME} ${ROOT_LIBRARIES} )

# Link each target with other targets or add options, etc.

# Adding something we can run - Output name matches target name
#add_executable(TransWarpExtraction macros/TransWarpExtraction.cxx)

# Make sure you link your targets with this command. It can also link libraries and
# even flags, so linking a target that does not exist will not give a configure-time error.
#target_link_libraries(TransWarpExtraction PRIVATE ${PROJECT_NAME} )
