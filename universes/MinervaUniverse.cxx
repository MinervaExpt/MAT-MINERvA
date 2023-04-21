#ifndef MINERVAUNIVERSE_cxx
#define MINERVAUNIVERSE_cxx

#include "universes/MinervaUniverse.h"

#include <algorithm>

//c++ includes
#include <string>
#include <cstdlib> //std::getenv()
#include <regex>
#include <fstream>
#include <stack>
#include <cassert>
#include <sstream>
#include <algorithm>

using namespace PlotUtils;

namespace {
bool _is_nue_constraint_set = false;
bool _is_playlist_set = false;
bool _is_horn_current_set = false;
bool _is_analysis_nu_pdg_set = false;
bool _is_nonResPi_reweight_set = false;
bool _is_D2_pion_genie_tune_set = false;
bool _is_zExpansionFa_reweight_set = false;
bool _is_readout_volume_set = false;
bool _is_mhrw_neutronCV_set = false;
bool _is_mhrw_elastics_set  = false;
bool _is_mhrw_filename_set  = false; 
}  // namespace

// Initialize variables
std::string MinervaUniverse::m_playlist = "minervame1a";
std::string MinervaUniverse::m_processing = "NX";
bool MinervaUniverse::m_is_FHC = false;
bool MinervaUniverse::m_is_RHC = false;
int MinervaUniverse::m_analysis_nu_pdg = 14;
bool MinervaUniverse::m_use_nuE_constraint = false;
int MinervaUniverse::m_n_flux_universes = 200;
bool MinervaUniverse::m_use_nonResPi_reweight = false;
bool MinervaUniverse::m_use_D2_pion_genie_tune = false;
bool MinervaUniverse::m_use_zExpansionFa_reweight = false;

double MinervaUniverse::m_muon_momentum_cv_offset = 0.0;
bool MinervaUniverse::_is_muon_momentum_cv_offset_set = false;
bool MinervaUniverse::_has_muon_error_been_thrown = false;

bool MinervaUniverse::m_mhrw_load_renorm            = true;
bool MinervaUniverse::m_use_mhrw_neutronCV_reweight = true;
bool MinervaUniverse::m_use_mhrw_elastics           = true;
double MinervaUniverse::m_readout_vol_minZ          = NSFDefaults::TrackerFace;
double MinervaUniverse::m_readout_vol_maxZ          = NSFDefaults::TrackerBack;
double MinervaUniverse::m_readout_vol_apothem       = NSFDefaults::StandardApothem;


double MinervaUniverse::m_true_proton_ke_cut_central = NSFDefaults::TrueProtonKECutCentral;
double MinervaUniverse::m_reco_proton_ke_cut_central = NSFDefaults::RecoProtonKECutCentral;

std::string MinervaUniverse::m_readout_volname      = "";
std::string MinervaUniverse::m_mhrw_process_name    = "Renorm_Kine_Truth";
std::string MinervaUniverse::m_mhrw_directory       = Form("%s/data/mhrwKineRenorm",getenv("PLOTUTILSROOT"));
std::string MinervaUniverse::m_treename= "MasterAnaDev";
// CTOR TODO:depends on the above conversation about whether all those bools
// need to be static or not

//! Playlists
std::string MinervaUniverse::GetPlaylist() {
  if (!_is_playlist_set) {
    std::cout << "WARNING: YOU ARE USING DEFAULT CVU PLAYLIST (minervame1a). "
              << "PLEASE SET THE VARIABLE BY CVU::SetPlaylist(playlist) "
              << "BEFORE USING UNIVERSES." << std::endl;
  }
  return m_playlist;
}

bool MinervaUniverse::IsPlaylistME(std::string playlist) {
  if (playlist.find("me") == std::string::npos &&
      playlist.find("ME") == std::string::npos) {
    return false;
  } else {
    return true;
  }
}

bool MinervaUniverse::SetPlaylist(std::string playlist) {
  if (_is_playlist_set) {
    std::cout << "WARNING: YOU ATTEMPTED SETTING THE PLAYLIST A SECOND TIME. "
              << "THIS IS NOT ALLOWED FOR CONSISTENCY." << std::endl;
    return false;
  } else {
    m_playlist = playlist;
    _is_playlist_set = true;
    if (!SetHornCurrent(playlist)) return false;
  }
  return true;
}

void MinervaUniverse::SetTreeName(std::string treename){m_treename=treename;}
std::string MinervaUniverse::GetTreeName(){return m_treename;}

//! MINERvA Software Processing
bool MinervaUniverse::SetProcessingEroica() {
  if (m_processing == "Eroica") {
    std::cout << "WARNING: THE PROCESSING IS ALREADY SET AS EROICA. "
              << "EXECUTING THIS AGAIN IS NOT ALLOWED FOR CONSISTENCY."
              << std::endl;
    return false;
  } else {
    m_processing = "Eroica";
  }
  return true;
}

bool MinervaUniverse::IsProcessingNX() {
  if (m_processing == "NX") {
    return true;
  } else if (m_processing == "Eroica") {
    return false;
  } else {
    std::cout << "WARNING: YOU HAVE AN UNRECOGNIZED PROCESSING SET!"
              << std::endl;
    return false;
  }
}

//! Horn current
bool MinervaUniverse::SetHornCurrent(std::string& playlist) {
  if (_is_horn_current_set) {
    std::cout
        << "WARNING: YOU ATTEMPTED SETTING THE HORN CURRENT A SECOND TIME. "
        << "You really shouldn't be trying to do this anyway"
        << "Horn current is set when you input a playlist"
        << "THIS IS NOT ALLOWED FOR CONSISTENCY." << std::endl;
    return false;
  } else {
    // case-independent
    std::transform(playlist.begin(), playlist.end(), playlist.begin(),
                   ::tolower);
    if (playlist.compare("minerva1") == 0 ||
        playlist.compare("minerva7") == 0 ||
        playlist.compare("minerva9") == 0 ||
        playlist.compare("minerva13") == 0 ||
        playlist.compare("minerva13c") == 0 ||
        playlist.compare("minerva13e") == 0 ||
        playlist.compare("minervame1a") == 0 ||
        playlist.compare("minervame1b") == 0 ||
        playlist.compare("minervame1c") == 0 ||
        playlist.compare("minervame1d") == 0 ||
        playlist.compare("minervame1e") == 0 ||
        playlist.compare("minervame1f") == 0 ||
        playlist.compare("minervame1g") == 0 ||
        playlist.compare("minervame1l") == 0 ||
        playlist.compare("minervame1m") == 0 ||
        playlist.compare("minervame1n") == 0 ||
        playlist.compare("minervame1o") == 0 ||
        playlist.compare("minervame1p") == 0)
      m_is_FHC = true;

    if (playlist.compare("minerva5") == 0 ||
        playlist.compare("minervame5a") == 0 ||
        playlist.compare("minervame6a") == 0 ||
        playlist.compare("minervame6b") == 0 ||
        playlist.compare("minervame6c") == 0 ||
        playlist.compare("minervame6d") == 0 ||
        playlist.compare("minervame6e") == 0 ||
        playlist.compare("minervame6g") == 0 ||
        playlist.compare("minervame6h") == 0 ||
        playlist.compare("minervame6i") == 0 ||
        playlist.compare("minervame6j") == 0)
      m_is_RHC = true;

    if (!m_is_FHC && !m_is_RHC)
      std::cout << "You have entered a nonstandard playlist. I can't "
                   "`SetHornCurrent`, and therefore I can't `SetPlaylist`!"
                << std::endl;
    _is_horn_current_set = true;
  }
  return true;
}

//! Analysis neutrino identity -- used for flux weights
int MinervaUniverse::GetAnalysisNuPDG() {
  if (!_is_analysis_nu_pdg_set) {
    std::cout << "WARNING: YOU ARE USING DEFAULT CVU NU PDG (14). "
              << "PLEASE SET THE VARIABLE BY CVU::SetAnalysisNuPDG(int pdg) "
              << "BEFORE USING UNIVERSES." << std::endl;
  }
  return m_analysis_nu_pdg;
}

bool MinervaUniverse::SetAnalysisNuPDG(int nu_pdg) {
  if (_is_analysis_nu_pdg_set) {
    std::cout << "WARNING: YOU ATTEMPTED SETTING THE NU PDG A SECOND TIME. "
                 "THIS IS NOT ALLOWED FOR CONSISTENCY."
              << std::endl;
    return false;
  } else {
    m_analysis_nu_pdg = nu_pdg;
    _is_analysis_nu_pdg_set = true;
  }
  return true;
}

bool MinervaUniverse::UseDeuteriumGeniePiTune() {
  if (!_is_D2_pion_genie_tune_set) {
    std::cout << "WARNING: YOU ARE USING DEFAULT DEUTERIUM GENIE PION TUNE (false). "
              << "PLEASE SET THE VARIABLE BY "
              << "CVU::SetDeuteriumGeniePiTune(bool use_D2_pion_genie_tune) "
              << "BEFORE USING UNIVERSES." << std::endl;
  }
  return m_use_D2_pion_genie_tune;
}

bool MinervaUniverse::SetDeuteriumGeniePiTune(bool use_tune) {
  if (_is_D2_pion_genie_tune_set) {
    std::cout
        << "WARNING: YOU ATTEMPTED SETTING DEUTERIUM GENIE PI TUNE A SECOND TIME. "
        << "THIS IS NOT ALLOWED FOR CONSISTENCY." << std::endl;
    return false;
  } else {
    m_use_D2_pion_genie_tune = use_tune;
    _is_D2_pion_genie_tune_set = true;
  }
  return true;
}

bool MinervaUniverse::UseNonResPiReweight() {
  if (!_is_nonResPi_reweight_set) {
    std::cout << "WARNING: YOU ARE USING DEFAULT NONRESPI REWEIGHT (false). "
              << "PLEASE SET THE VARIABLE BY "
              << "CVU::SetNonResPiReweight(bool use_nonResPi_reweight) "
              << "BEFORE USING UNIVERSES." << std::endl;
  }
  return m_use_nonResPi_reweight;
}

bool MinervaUniverse::SetNonResPiReweight(bool use_reweight) {
  if (_is_nonResPi_reweight_set) {
    std::cout
        << "WARNING: YOU ATTEMPTED SETTING NONRESPI REWEIGHT A SECOND TIME. "
        << "THIS IS NOT ALLOWED FOR CONSISTENCY." << std::endl;
    return false;
  } else {
    m_use_nonResPi_reweight = use_reweight;
    _is_nonResPi_reweight_set = true;
  }
  return true;
}

bool MinervaUniverse::UseZExpansionFaReweight() {
  // Debug message that can be activated if/when this feature becomes more
  // widely used
  if (!_is_zExpansionFa_reweight_set) {
    std::cout << "WARNING: YOU ARE USING DEFAULT ZEXPANSION REWEIGHT (false). "
              << "PLEASE SET THE VARIABLE BY "
              << "CVU::SetZExpansionFaReweight(bool use_zExpansionFa_reweight) "
              << "BEFORE USING UNIVERSES." << std::endl;
  }
  return m_use_zExpansionFa_reweight;
}

bool MinervaUniverse::SetZExpansionFaReweight(bool use_reweight) {
  if (_is_zExpansionFa_reweight_set) {
    std::cout
        << "WARNING: YOU ATTEMPTED SETTING ZEXPANSION REWEIGHT A SECOND TIME. "
        << "THIS IS NOT ALLOWED FOR CONSISTENCY." << std::endl;
    return false;
  } else {
    m_use_zExpansionFa_reweight = use_reweight;
    _is_zExpansionFa_reweight_set = true;
  }
  return true;
}

//
bool MinervaUniverse::UseNuEConstraint() {
  if (!_is_nue_constraint_set) {
    std::cout << "WARNING: YOU ARE USING DEFAULT CVU NUE CONSTRAINT (false). "
              << "PLEASE SET THE VARIABLE BY "
              << "CVU::SetNuEConstraint(bool use_nue_constraint) "
              << "BEFORE USING UNIVERSES." << std::endl;
  }
  return m_use_nuE_constraint;
}

bool MinervaUniverse::SetNuEConstraint(bool use_constraint) {
  if (_is_nue_constraint_set) {
    std::cout << "WARNING: YOU ATTEMPTED SETTING NUE CONSTRAINT A SECOND TIME. "
              << "THIS IS NOT ALLOWED FOR CONSISTENCY." << std::endl;
    return false;
  } else {
    m_use_nuE_constraint = use_constraint;
    _is_nue_constraint_set = true;
  }
  return true;
}

// This functionality should *not* be used by standard users. This will allow
// data preservation users in the future to offset the muon energy scale in data
// when it is no longer feasible to reproduce ntuples. So, the warning is
// inverted in that the user will be warned only if they _do_ use the method to
// specify an offset. n.b. 'offset' refers to CV only, and 'shift' refers to
// systematic universes.

int MinervaUniverse::GetNFluxUniverses() { return m_n_flux_universes; }

bool MinervaUniverse::SetNFluxUniverses(int n_flux_universes) {
  m_n_flux_universes = n_flux_universes;
  return true;
}

bool MinervaUniverse::SetMHRWeightFilename( std::string process_name, std::string directory ) {
  if( _is_mhrw_filename_set ) {  
    std::cout << "WARNING: YOU ATTEMPTED TO SET MHRW FILENAME A SECOND TIME. "
              << "THIS IS NOT ALLOWED FOR CONSISTENCY." << std::endl;
    return false ; 
  } else {
    m_mhrw_process_name = process_name;
    m_mhrw_directory    = directory;
    _is_mhrw_filename_set = true;
  }
  return true;
}

bool MinervaUniverse::SetReadoutVolume( std::string readout_vol ) {
  m_readout_volname = readout_vol;
  return true;
}

bool MinervaUniverse::SetReadoutVolume( double minZ, double maxZ, double apothem ) {
  if( _is_readout_volume_set ) {
    std::cout << "WARNING: YOU ATTEMPTED TO SET READOUT VOLUME A SECOND TIME. "
              << "THIS IS NOT ALLOWED FOR CONSISTENCY." << std::endl;
    return false ; 
  } else {
    m_readout_vol_minZ    = minZ;
    m_readout_vol_maxZ    = maxZ;
    m_readout_vol_apothem = apothem;
    _is_readout_volume_set = true;
  }
  return true;
}

bool MinervaUniverse::SetMHRWeightNeutronCVReweight( bool useNeutronCVReweight ) {
  if( _is_mhrw_neutronCV_set ) {
    std::cout << "WARNING: YOU ATTEMPTED TO SET MHRW NEUTRON CV REWEIGHT A SECOND TIME. "
              << "THIS IS NOT ALLOWED FOR CONSISTENCY." << std::endl;
    return false ; 
  } else {
    m_use_mhrw_neutronCV_reweight = useNeutronCVReweight;
    _is_mhrw_neutronCV_set = true;
  }
  return true;
}

bool MinervaUniverse::SetMHRWeightElastics( bool useElastics ) {
  if( _is_mhrw_elastics_set ) { 
    std::cout << "WARNING: YOU ATTEMPTED TO SET MHRW ELASTICS A SECOND TIME. "
              << "THIS IS NOT ALLOWED FOR CONSISTENCY." << std::endl;
    return false ; 
  } else {
    m_use_mhrw_elastics = useElastics;
    _is_mhrw_elastics_set = true;
  }
  return true;
}

//This is more for convenience.  The user shouldn't have to call this  
bool MinervaUniverse::SetupMHRWeighter() const {
  if( m_mhrw_load_renorm ) 
  {
    // Make sure to init the static MHRW 
    weight_hadron<PlotUtils::TreeWrapper*>(m_chw, GetPlaylist(), 
                                        m_readout_vol_minZ, m_readout_vol_maxZ, m_readout_vol_apothem, 
                                        m_use_mhrw_neutronCV_reweight, m_use_mhrw_elastics,
                                        m_mhrw_process_name, m_mhrw_directory);
    //If readout volume isn't empty, use that, and change the volume variables
    if( !m_readout_volname.empty() )
    {
      weight_hadron<PlotUtils::TreeWrapper*>().setReadoutVolume(m_readout_volname);

      double minZ = 0, maxZ = 0, apothem = 0;
      weight_hadron<PlotUtils::TreeWrapper*>().getBasicFiducial(minZ,maxZ,apothem);;
      SetReadoutVolume( minZ, maxZ, apothem );
      m_readout_vol_minZ    = minZ   ;  
      m_readout_vol_maxZ    = maxZ   ;  
      m_readout_vol_apothem = apothem;  
    }

    weight_hadron<PlotUtils::TreeWrapper*>().getTruthKineRenorm(); 
    m_mhrw_load_renorm = false; 
  }
  return true;
}

// These should probably move, but I'm not sure where, yet.
//! override these to implement beam angle universes

double MinervaUniverse::GetBeamAngleOffsetX() const {
  return NSFDefaults::beamThetaX_Central;
}

double MinervaUniverse::GetBeamAngleOffsetY() const {
  return NSFDefaults::beamThetaY_Central;
}

double MinervaUniverse::GetBatchPOT() const {
  double batch_pot = 0;

  // Anushree bash pot calculation
  double numi_pot = GetDouble("numi_pot");
  bool is_me = MinervaUniverse::IsPlaylistME(MinervaUniverse::GetPlaylist());
  std::string channel = MinervaUniverse::GetRecoilChannel();
  int batch_structure;
  int reco_vertex_batch;
  if (is_me) {
    batch_structure = GetInt("batch_structure");
    reco_vertex_batch = GetInt("reco_vertex_batch");
  } else {
    batch_structure = 0;
  }

  if (batch_structure == 0) {
    batch_pot = numi_pot / 6.;
  }
  if (batch_structure == 1) {
    if (reco_vertex_batch < 3) {
      batch_pot = numi_pot / 4.;
    } else {
      batch_pot = numi_pot / 8.;
    }
  }
  if (batch_structure == 2) {
    if (reco_vertex_batch < 5) {
      batch_pot = numi_pot / 5.;
    } else {
      batch_pot = numi_pot / 10.;
    }
  }
  if (batch_structure == 3) {
    batch_pot = numi_pot / 6.;
  }
  if (batch_structure == -1) {
    // std::cout << "There isn't a batch structure for this gate, "
    //          << "so I'm not going to guess the batch POT if the batch
    //          structure is 0. Entry = "
    //          << m_entry << std::endl;
    batch_pot = numi_pot / 6.;
  }
  return batch_pot;
}

std::string MinervaUniverse::GetRecoilChannel() const
  {
    int targetID = GetInt((m_treename +"_ANN_targetID").c_str());
    int targetZ = GetInt((m_treename +"_ANN_targetZ").c_str());
    int helicity = m_analysis_nu_pdg;
    bool m_useCCTrackerTuning = 0;
    std::string channel;

    std::string prefix = "NukeCC_";
    if(       -14 == m_analysis_nu_pdg)
      prefix += "AntiNu_";
    else if( 14 == m_analysis_nu_pdg )
      prefix += "Nu_";

    if(      1 == targetID ) channel = prefix + "Tgt1";
    else if( 2 == targetID ) channel = prefix + "Tgt2";
    else if( 3 == targetID ) channel = prefix + "Tgt3";
    else if( 4 == targetID ) channel = prefix + "Tgt4";
    else if( 5 == targetID ) channel = prefix + "Tgt5";
    else
      {
        if( m_useCCTrackerTuning == 0 )
          channel = ( -14 == m_analysis_nu_pdg ) ? "AntiCCInclusive" : "CCInclusive";
        else
          channel = prefix + "Tracker";
      }

    if( 1 <= targetID && targetID <= 5 )
      {
        if( 6  == targetZ ) channel += "_C";
        if( 26 == targetZ ) channel += "_Fe";
        if( 82 == targetZ ) channel += "_Pb";
      }
    return channel;
  }

#define FORBID(KEYWORD, STATE)\
  assert(line.find("KEYWORD") == std::string::npos && KEYWORD "keyword does not make sense in the " STATE "state!");
std::string& trim(std::string& toTrim)
  {
    const size_t firstWhitespace = toTrim.find_first_not_of(" ");
    if(firstWhitespace != std::string::npos) toTrim = toTrim.substr(firstWhitespace);

    toTrim.erase(toTrim.find_last_not_of(" ")+1, std::string::npos);
    
    return toTrim;
  }

  bool foundKeyword(const std::string& line, const std::string& keyword, size_t& nextToken)
  {
    nextToken = line.find(keyword);
    if(nextToken != std::string::npos) nextToken += keyword.length() + 1;
    return nextToken != std::string::npos;
  }

  std::string& replaceEnvVars(std::string& path)
  {
    std::regex envVar(R"(\$[[:alnum:]]+)");
    std::smatch results;

    while(std::regex_search(path, results, envVar))
    {
      auto toReplace = results.str();
      for(size_t found = toReplace.find("{"); found != std::string::npos; found = toReplace.find("{")) toReplace.erase(found, 1);
      for(size_t found = toReplace.find("}"); found != std::string::npos; found = toReplace.find("}")) toReplace.erase(found, 1);

      const char* env = std::getenv(toReplace.substr(1).c_str()); //Important to skip the $
      if(!env) env = "";

      for(auto& result: results) path.replace(std::distance(path.cbegin(), result.first), std::distance(result.first, result.second), env);
    }

    return path;
  }

  std::unordered_map<std::string, MinervaUniverse::CaloCorrection> MinervaUniverse::CaloCorrection::parse(const std::string& caloFile)
  {
    std::unordered_map<std::string, MinervaUniverse::CaloCorrection> corrections;
    auto currentCorr = corrections.end();

    //Model this parsing job as a state machine
    enum class State
    {
      Top, //When not in a BEGIN/END pair, look for BEGIN and IMPORT
      Loading //When between a BEGIN and an END, look for POLYPOINTs
    };

    //Theorem: Any recurive algorithm can be rewritten using a LIFO
    std::stack<std::string> filesLeft;
    filesLeft.push(caloFile);

    while(!filesLeft.empty())
    {
      std::ifstream currentFile(replaceEnvVars(filesLeft.top()));
      assert(currentFile.is_open() && "Failed to open a calorimetry spline file!");
      filesLeft.pop();
      //--depth; //TODO: This depth check doesn't seem correct because I keep going with the current file

      //int depth = 0;
      State parseState = State::Top;
      std::string line;
      while(std::getline(currentFile, line))
      {
        line = trim(line);
        if(line.empty() || line[0] == '#') continue;

        //Parser state machine.
        //Keywords have the following priorities because of how CalorimetryUtils
        //is written.
        size_t nextToken = std::string::npos;
        if(parseState == State::Top)
        {
          if(foundKeyword(line, "BEGIN", nextToken))
          {
            currentCorr = corrections.insert(std::make_pair(line.substr(nextToken), CaloCorrection())).first;
            parseState = State::Loading;
          }
          else if(foundKeyword(line, "IMPORT", nextToken)) 
          {
            filesLeft.push(line.substr(nextToken));
            //++depth; //TODO: This is just filesLeft.size(), not recursion depth.  Maybe I could use CHANGES in filesLeft.size() to track recursion depth?  When it decreases, --depth.  Keep a bool, addedFiles, and ++depth when it's true at the END of this loop.
            //if(depth > 5) throw std::runtime_error("Reached maximum recursion depth of 5 in parseCaloCorrections().  You've probably loaded a calorimetry polyline with recursive IMPORTs."); 
          }
          else
          {
            #ifndef NDEBUG
            std::cout << "In parseCaloCorrections(), ignoring a line that looks like " << line << " because it doesn't match any keywords for State::Top.\n";
            #endif //NDEBUG

            FORBID("END", "Top")
            FORBID("POLYPOINT", "Top")
          }
        }
        else if(parseState == State::Loading)
        {
          if(foundKeyword(line, "POLYPOINT", nextToken))
          {
            std::stringstream findPoint(line.substr(nextToken));
            double x, y;
            findPoint >> x >> y;
            currentCorr->second.fPoints.push_back({x, y});
          }
          else if(foundKeyword(line, "END", nextToken))
          {
            std::sort(currentCorr->second.fPoints.begin(), currentCorr->second.fPoints.end(), [](const auto& lhs, const auto& rhs) { return lhs.threshold < rhs.threshold; });
            parseState = State::Top;
          }
          else if(foundKeyword(line, "SCALE", nextToken))
          {
            std::stringstream findScale(line.substr(nextToken));
            findScale >> currentCorr->second.fScale;
          }
          else
          {
            #ifndef NDEBUG
            std::cout << "In parseCaloCorrections(), ignoring a line that looks like " << line << " because it doesn't match any keywords for State::Loading.\n";
            #endif //NDEBUG

            FORBID("BEGIN", "Loading")
            FORBID("IMPORT", "Loading")
          }
        }
        else assert(false && "Parser got into unknown state!");
      } //End loop over currentFile
    } //End while filesLeft

    return corrections;
  }

  MinervaUniverse::CaloCorrection::CaloCorrection(const std::string& caloFile, const std::string& tuningName)
  {
    //std::string caloTune;
    //if (tuningName == "NukeCCInclusiveTunings"){
     // caloTune = "NukeCC_AntiNu_Tgt5_Fe";
      //univ->GetRecoilChannel();
    //}
    const auto parsed = parse(caloFile);
    const auto found = parsed.find(tuningName);
    if(found != parsed.end())
    {
      fPoints = found->second.fPoints;
      fScale = found->second.fScale;
    }
    else
    {
      std::cerr << "Calorimetric splines found:\n";
      for(const auto& spline: parsed) std::cerr << spline.first << ".\n";

      throw std::runtime_error("Failed to find a calorimetric spline named " + tuningName + ".");
    }

    #ifndef NDEBUG
      std::cout << "Using a calorimetric correction named " << tuningName << " with scale " << fScale << " from " << caloFile << ":\n";
      for(const auto& point: fPoints) std::cout << point.threshold << ", " << point.correction << "\n";
    #endif
  }

  double MinervaUniverse::CaloCorrection::eCorrection(const double rawRecoil) const
  {
    const double scaledRecoil = rawRecoil * fScale;
    const auto upperPoint = std::lower_bound(std::next(fPoints.begin()), fPoints.end(), scaledRecoil, [](const auto& point, const auto recoil) { return point.threshold < recoil; });
    if(upperPoint != fPoints.end())
    {
      const auto lowerPoint = std::prev(upperPoint);
      return std::max(0., lowerPoint->correction + (upperPoint->correction - lowerPoint->correction) * (scaledRecoil - lowerPoint->threshold) / (upperPoint->threshold - lowerPoint->threshold) );
    }
    else return scaledRecoil; //scaledRecoil is after last point, so don't correct.
                              //This is both a fail-safe mechanism for very energetic
                              //neutrinos which are very rare in MINERvA's <6 GeV>
                              //beam and the Default spline.
  }

  MinervaUniverse::CaloCorrection::CaloCorrection(): fScale(0), fPoints{}
  {
  }

#endif  // MinervaUniverse_cxx
