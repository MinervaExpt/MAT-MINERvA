#ifndef MADMACROUTIL_H 
#define MADMACROUTIL_H

#include "PlotUtils/MacroUtil.h"
#include <algorithm>  // transform

namespace MAD {
class MacroUtil : public PlotUtils::MacroUtil {
 public:
  // TODO this is currently using Ever's and Ben's personal p3 lists.
  // This should obviously be made official at some point.
  // Ben's idea for a path forward: commit official playlist_filelists.txts to
  // github and have this function read that. Additionally, we'll want a python
  // version of that function as well.
  std::string GetPlaylistFile(std::string plist, const bool is_mc,
                              const bool do_test_playlist,
                              const bool use_xrootd,
                              const std::string mad_version_tag) {
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

  // Construct a MacroUtil with MAD tuples from just a playlist name.
  //
  // Use the above GetPlaylistFile function above to get the file lists.
  //
  // Please focus on the function signature and excuse the monstrosity that is
  // the initializer list.
  MacroUtil(const std::string& plist, const bool do_mc = true,
            const bool do_truth = false, const bool do_data = false,
            const bool do_test_playlist = false, const bool use_xrootd = true,
            const std::string mad_version_tag = "p3");
  /*
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

  bool m_do_mc;
  bool m_do_truth;
  bool m_do_data;
};
}  // namespace MAD

#endif
