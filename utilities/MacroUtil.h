#ifndef MADMACROUTIL_H 
#define MADMACROUTIL_H

#include "PlotUtils/MacroUtil.h"
#include <algorithm>  // transform

namespace PlotUtils {
PlotUtils::MacroUtil MakeMacroUtil(const std::string& plist, const bool do_mc = true,
          const bool do_truth = false, const bool do_data = false,
          const bool do_test_playlist = false, const bool use_xrootd = true,
          const std::string mad_version_tag = "p3") {
}
} // namespace PlotUtils

/*
namespace MAD {
class MacroUtil : public PlotUtils::MacroUtil {
 public:

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

  //    : do_mc
  //      ? do_data ? PlotUtils::MacroUtil(
  //                      "MasterAnaDev",
  //                      GetPlaylistFile(plist, true, do_test_playlist,
  //                                      use_xrootd, mad_version_tag),
  //                      GetPlaylistFile(plist, false, do_test_playlist,
  //                                      use_xrootd, mad_version_tag),
  //                      plist,
  //                      do_truth)  // mc + data
  //                : PlotUtils::MacroUtil(
  //                      "MasterAnaDev",
  //                      GetPlaylistFile(plist, true, do_test_playlist,
  //                                      use_xrootd, mad_version_tag),
  //                      plist,
  //                      false)  // mc-only
  //      : PlotUtils::MacroUtil("MasterAnaDev",
  //                             GetPlaylistFile(plist, false, do_test_playlist,
  //                                             use_xrootd, mad_version_tag),
  //                             plist),  // data only
  //m_do_mc(do_mc), m_do_truth(do_truth), m_do_data(do_data) {}

  bool m_do_mc;
  bool m_do_truth;
  bool m_do_data;
};
}  // namespace MAD
*/

#endif
