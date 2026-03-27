#ifndef MAKEMACROUTIL_H
#define MAKEMACROUTIL_H

#include <algorithm>  // transform
#include <iostream>

#include "PlotUtils/MacroUtil.h"

namespace PlotUtils {
// TODO this is currently using Ever's and Ben's personal p3 lists.
// This should obviously be made official at some point. Committed. 
// TODO Python version of this function.
std::string GetPlaylistFile(std::string plist, const bool is_mc,
                            const bool do_test_playlist,
                            const bool use_xrootd,
                            const std::string& mad_version_tag) {
  if (do_test_playlist) {
    return is_mc ? "/minerva/app/users/bmesserl/MATAna/cc-ch-pip-ana/cache/"
                   "mc_ME1A_p3.txt"
                 : "/minerva/app/users/bmesserl/MATAna/cc-ch-pip-ana/cache/"
                   "data_ME1A_p3.txt";
  }
  std::transform(plist.begin(), plist.end(), plist.begin(), ::toupper);
  const std::string processing_date = "production_p3";
  const std::string is_mc_str = is_mc ? "mc" : "data";
  std::string topdir = is_mc ? "/minerva/data/users/granados/MAD_ana_plists/"
                             : "/minerva/data/users/granados/MAD_ana_plists/";
  topdir += processing_date;
  std::string playlist_file =
      use_xrootd ? Form("%s/%s_%s_xrootd_plist.txt", topdir.c_str(),
                        is_mc_str.c_str(), plist.c_str())
                 : Form("%s/%s_%s_plist.txt", topdir.c_str(),
                        is_mc_str.c_str(), plist.c_str());
  return playlist_file;
}

PlotUtils::MacroUtil MakeMacroUtil(const std::string plist,
                                   const bool do_mc,
                                   const bool do_truth = false,
                                   const bool do_data = false,
                                   const bool do_test_playlist = false,
                                   const bool use_xrootd = true,
                                   const std::string& mad_version_tag = "p3") {
  std::string mc_file_list = GetPlaylistFile(plist, true, do_test_playlist, use_xrootd, mad_version_tag);
  std::string data_file_list = GetPlaylistFile(plist, false, do_test_playlist, use_xrootd, mad_version_tag);

  if (do_mc && do_data) {
    return PlotUtils::MacroUtil("MasterAnaDev", mc_file_list, data_file_list, plist, do_truth) :
  } else if (do_mc) {
    return PlotUtils::MacroUtil( "MasterAnaDev", mc_file_list, plist, do_truth) :
  } else if (do_data) {
    PlotUtils::MacroUtil("MasterAnaDev", data_file_list, plist) :
  } else {
    std::cerr << "MakeMacroUtil: Data/MC not specified. Here's a mc-only MacroUtil.\n";
    return PlotUtils::MacroUtil( "MasterAnaDev", mc_file_list, plist, do_truth) :
  }
}
}  // namespace PlotUtils

#endif
