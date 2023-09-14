#ifndef MADMACROUTIL_CXX
#define MADMACROUTIL_CXX

#include "utilities/MacroUtil.h"

/*
MAD::MacroUtil::MacroUtil(const std::string& plist, const bool do_mc,
          const bool do_truth, const bool do_data,
          const bool do_test_playlist, const bool use_xrootd,
          const std::string mad_version_tag)
    : do_mc
      ? do_data ? PlotUtils::MacroUtil(
                      "MasterAnaDev",
                      GetPlaylistFile(plist, true, do_test_playlist,
                                      use_xrootd, mad_version_tag),
                      GetPlaylistFile(plist, false, do_test_playlist,
                                      use_xrootd, mad_version_tag),
                      plist,
                      do_truth)  // mc + data
                : PlotUtils::MacroUtil(
                      "MasterAnaDev",
                      GetPlaylistFile(plist, true, do_test_playlist,
                                      use_xrootd, mad_version_tag),
                      plist,
                      false)  // mc-only
      : PlotUtils::MacroUtil("MasterAnaDev",
                             GetPlaylistFile(plist, false, do_test_playlist,
                                             use_xrootd, mad_version_tag),
                             plist),  // data only
    m_do_mc(do_mc), m_do_truth(do_truth), m_do_data(do_data) {}
*/

#endif
