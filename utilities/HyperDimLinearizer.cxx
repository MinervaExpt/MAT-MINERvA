#ifndef MNV_HYPERDIMLINEARIZER_cxx
#define MNV_HYPERDIMLINEARIZER_cxx 1
#include "utilities/HyperDimLinearizer.h"

#include <string>

namespace PlotUtils {

//! Get a set of bins
// ==========================================================================
// CTOR
// ==========================================================================
PlotUtils::HyperDimLinearizer::HyperDimLinearizer(std::vector<std::vector<double>> input, int type)
    : m_invec(input) {
    if (type == 0) {
        m_analysis_type = k2D;
    } else if (type == 1) {
        m_analysis_type = k1D;
    } else if (type == 2) {
        m_analysis_type = k2D_lite;
    } else if (type == 3) {
        m_analysis_type = k1D_lite;
    }
    // m_invec = input;
    int n_dim = input.size();
    // m_analysis_type = type;
    std::string n_lindim;
    if (m_analysis_type == k1D || m_analysis_type == k1D_lite) {
        n_lindim = "1D";
    } else {
        n_lindim = "2D";
    }
    std::cout << "Contructing class with " << n_dim << " dimensions of type " << m_analysis_type << " linearized to " << n_lindim << std::endl;
    if (m_analysis_type == k2D_lite || m_analysis_type == k1D_lite)
        std::cout << " Selected lightweight type. Under/overflow bins will not be counted." << std::endl;

    m_cell_size = {1};      // vector to put sizes of each cell for each bin in each axis. Start with x-axis which is always 1 bin.
    m_n_global_x_bins = 1;  // How  many bins you have on your linearized axis

    for (unsigned int i = 0; i < input.size(); i++) {
        // Figure out size of each axis in phase space depending on analysis type
        int tmp_el_size;
        if (m_analysis_type == k2D || m_analysis_type == k1D) {  // Counting under/overflow for each axis number of bins = (vector size - 1 + 2)
            tmp_el_size = input[i].size() + 1;
        } else {  // For lite types don't coun't under/overflow, number of bins = (vector size - 1)
            tmp_el_size = input[i].size() - 1;
        }
        std::cout << "Bin number " << i << "\t" << tmp_el_size << std::endl;
        m_el_size.push_back(tmp_el_size);

        if (m_analysis_type == k2D || m_analysis_type == k2D_lite) {
            if (i == 1)
                continue;  // Skip all these steps for y-axis if doing 2D
            if (i == 0)
                m_cell_size.push_back(0);  // this takes care of y-axis for 2D
        }
        m_n_global_x_bins *= tmp_el_size;  // Count the number of bins in linearized space and size of cells for each axis

        if (i < input.size() - 1)                                 // skip last axis since we start with 1 for x already.
            m_cell_size.push_back(tmp_el_size * m_cell_size[i]);  // cell size for axis i is c_i = n_(i-1)*c_(i-1)
    }                                                             // close loop over input
}  // close constructor
//    // Jank to make things build right for now. TODO: Make less jank
//    HyperDimLinearizer(input, m_analysis_type);
// }

PlotUtils::HyperDimLinearizer::HyperDimLinearizer(std::vector<std::vector<double>> input, EAnalysisType type)
    : m_invec(input),
      m_analysis_type(type) {
    // m_invec = input;
    int n_dim = input.size();
    // m_analysis_type = type;
    std::string n_lindim;
    if (type == k1D || type == k1D_lite) {
        n_lindim = "1D";
    } else {
        n_lindim = "2D";
    }
    std::cout << "Contructing class with " << n_dim << " dimensions of type " << type << " linearized to " << n_lindim << std::endl;
    if (type == k2D_lite || type == k1D_lite)
        std::cout << " Selected lightweight type. Under/overflow bins will not be counted." << std::endl;

    m_cell_size = {1};      // vector to put sizes of each cell for each bin in each axis. Start with x-axis which is always 1 bin.
    m_n_global_x_bins = 1;  // How  many bins you have on your linearized axis

    for (unsigned int i = 0; i < input.size(); i++) {
        // Figure out size of each axis in phase space depending on analysis type
        int tmp_el_size;
        if (type == k2D || type == k1D) {  // Counting under/overflow for each axis number of bins = (vector size - 1 + 2)
            tmp_el_size = input[i].size() + 1;
        } else {  // For lite types don't coun't under/overflow, number of bins = (vector size - 1)
            tmp_el_size = input[i].size() - 1;
        }
        std::cout << "Bin number " << i << "\t" << tmp_el_size << std::endl;
        m_el_size.push_back(tmp_el_size);

        if (type == k2D || type == k2D_lite) {
            if (i == 1)
                continue;  // Skip all these steps for y-axis if doing 2D
            if (i == 0)
                m_cell_size.push_back(0);  // this takes care of y-axis for 2D
        }
        m_n_global_x_bins *= tmp_el_size;  // Count the number of bins in linearized space and size of cells for each axis

        if (i < input.size() - 1)                                 // skip last axis since we start with 1 for x already.
            m_cell_size.push_back(tmp_el_size * m_cell_size[i]);  // cell size for axis i is c_i = n_(i-1)*c_(i-1)
    }                                                             // close loop over input
}  // close constructor

// ==========================================================================
// Getter values between spaces
// ==========================================================================

// Transform values in phase space (in units of each axis) to linearized bin values.
//     Return ranges from [0, ((# lin bins) - 1)] inclusively
//     For lite types, additional under/overflow bins placed at high end of range from [(# lin bins), ((# lin bins) + (3^n_dim) - 1)] inclusively
//     For 2D types, y axis is binned normally; with 0 being underflow, ((# y bins) + 1) being overflow. Return 0 for y value on 1D types.
std::pair<int, int> PlotUtils::HyperDimLinearizer::GetBin(std::vector<double> values) {
    int global_x = 0;  // Returned linearized bin
    int y_bin = 0;
    
    // For non-lite types, just add each number of cells
    if (m_analysis_type == k2D || m_analysis_type == k1D) {  
        for (unsigned int i = 0; i < values.size(); i++) {
            global_x += Get1DBin(values[i], i) * m_cell_size[i];
        }
    // For lite types under/overflow get special treatment, so need to do things a lil differently
    } else if (m_analysis_type == k2D_lite || m_analysis_type == k1D_lite) {
        int tmp_global_x = 0;
        bool underover_bool = false;
        int flow_x = 0;     // Bin number in the under/overflow bins at end of linearized x axis
        int flow_cell = 1;  // Cell size for over/underflow bins (gets incremented in the loop)
        for (unsigned int i = 0; i < values.size(); i++) {
            if (i == 1 && m_analysis_type == k2D_lite)  // Skip y axis if doing 2D to avoid under/overflow issues
                continue;
            int tmp_bin = Get1DBin(values[i], i);
            if (tmp_bin == 0) {  // If underflow
                // flow_x += 0;  // (# of flow cells - 1) * (flow cell size), but # of flow cells is always 0 here
                underover_bool = true;
            } else if (tmp_bin == m_el_size[i] + 1) {  // If overflow
                flow_x += 2 * flow_cell;               
                underover_bool = true;
            } else {                                             // If normally binned
                flow_x += 1 * flow_cell;                         // (# of flow cells - 1) * (flow cell size)
                tmp_global_x += (tmp_bin - 1) * m_cell_size[i];  // Add to global_x while we're here, tmp_bin - 1 since we don't have underflow like it was before
            }
            flow_cell *= 3;  // Bump up to next size of flow cells, accounts for underflow, normalflow, & overflow for each axis
        }
        if (!underover_bool) {  // Event was normally binned on all axes
            global_x = tmp_global_x;
        } else {                                    // If any number of axes were in under/overflow
            global_x = m_n_global_x_bins + flow_x;  // all under/overflow gets put out at the end.
        }
    }

    if (m_analysis_type == k2D || m_analysis_type == k2D_lite)  // Get the y-bin if doing 2D
        y_bin = Get1DBin(values[1], 1);
    std::pair<int, int> lin_bin = std::make_pair(global_x, y_bin);  // Make them a pair. If doing 1D, should get 0 there
    return lin_bin;
}

// Find out what bin a value is in for a given axis (indexed by 'el')
int PlotUtils::HyperDimLinearizer::Get1DBin(double value, int el) {
    int b = 0;
    for (unsigned int i = 0; i < m_invec[el].size(); i++) {  // loop over bin boundaries
        if (value < m_invec[el][i]) {
            break;
        }
        b += 1;  // didn't find the bin, add 1. Underflow is 0 and overflow is size()+1
    }
    return b;
}

// Transform bin space coordinate to phase space coordinates (in units of bin index in phase space)
std::vector<int> PlotUtils::HyperDimLinearizer::GetValues(int x_linbin, int y_bin) {  //  Default ybin to 0 to maintain behaviour from Dan's version
    std::vector<int> ps_bin_coords;                                        // Phase space bin coordinates, returned here

    if (x_linbin > m_n_global_x_bins) {  // TODO: what should you do if under/overflow for lite types cases?
        std::cout << "HyperDimLinearizer::GetValues: WARNING: requested values out of range of linearized space. Result might not make sense." << std::endl;
        std::cout << "                                        requested value: " << x_linbin << ", \t range of linearized space: " << m_n_global_x_bins << std::endl;
        // TODO do lite types?
    }

    int mod_bin = x_linbin - 1;                                                   // Place holder, minus one because it won't count the x axis bins correctly otherwise.
    for (unsigned int i = 0; i < m_invec.size(); i++) {                           // Loop over coordinate axes
        if ((m_analysis_type == k2D || m_analysis_type == k2D_lite) && i == 1) {  // If doing 2D, put in y-bin and skip the other steps for y-axis.
            ps_bin_coords.push_back(y_bin);
            continue;
        }
        int val = (mod_bin / m_cell_size[i]) % m_el_size[i];  // How many cells (which translates to bins in phase space) on this axis do you have?
        if ((m_analysis_type == k2D || m_analysis_type == k1D)) {
            ps_bin_coords.push_back(val);
        } else {  // For lite types, need to index phase space bins up one since under/overflow is binned differently
            ps_bin_coords.push_back(val + 1);
        }
        mod_bin += -(val * m_cell_size[i]);  // Trim off cells you just counted. Puts you on the "left side" of a cell for next axis, prevents rounding issues from division
    }

    return ps_bin_coords;
}

double PlotUtils::HyperDimLinearizer::GetBinVolume(int x_linbin) {
    double ps_bin_vol = 1.;
    // Given the linearized bin number, get corresponding bin index in phase space coordinates
    std::vector<int> ps_coords = GetValues(x_linbin);  // For 2D type, this puts 0 for y-axis.
    for (int i = 0; i < ps_coords.size(); i++) {
        if ((m_analysis_type == k2D || m_analysis_type == k2D_lite) && i == 1)  // Skip y-axis if doing 2D type
            continue;

        int var_bin = ps_coords[i]; // Bin index for that axis
        if (var_bin == 0 || var_bin == m_invec[i].size()) {  // Break if under/overflow
            std::cout << "HyperDimLinearizer::GetBinVolume: WARNING bin requested is under/overflow. Returning 1." << std::endl;
            return 1.;
        }
        std::vector<double> var_binning = m_invec[i];  // Get the binning for that axis
        double bin_width = var_binning[var_bin] - var_binning[var_bin - 1];
        if (bin_width < 0)
            std::cout << "HyperDimLinearizer::GetBinVolume: WARNING, negative binwidth, will affect bin volume" << std::endl;
        else if (bin_width == 0)
            std::cout << "HyperDimLinearizer::GetBinVolume: WARNING, 0 binwidth, will affect bin volume" << std::endl;
        ps_bin_vol *= bin_width;
    }
    return ps_bin_vol;
}

double PlotUtils::HyperDimLinearizer::GetBinVolume(std::vector<int> ps_bin_coords, bool IncludeX = true) { // Get the volume of the bin in phase space from the phase space bin coordinates
    double ps_bin_vol = 1.;
    int start = 0;
    if (!IncludeX)  // Skip x-axis if requested (user must scale by width)
        start = 1;
    for (int i = start; i < ps_bin_coords.size(); i++) {

        if ((m_analysis_type == k2D || m_analysis_type == k2D_lite) && i == 1)  // Skip y-axis if doing 2D type
            continue;

        int var_bin = ps_bin_coords[i];                      // Bin index for that axis
        if (var_bin == 0 || var_bin == m_invec[i].size()) {  // Break if under/overflow
            std::cout << "HyperDimLinearizer::GetBinVolume: WARNING bin requested" << var_bin << "is under/overflow. Returning 1." << m_invec[i].size() << std::endl;
            return 1.;
        }
        std::vector<double> var_binning = m_invec[i];  // Get the binning for that axis
        double bin_width = var_binning[var_bin] - var_binning[var_bin - 1];
        if (bin_width < 0)
            std::cout << "HyperDimLinearizer::GetBinVolume: WARNING, negative binwidth, will affect bin volume" << std::endl;
        else if (bin_width == 0)
            std::cout << "HyperDimLinearizer::GetBinVolume: WARNING, 0 binwidth, will affect bin volume" << std::endl;
        ps_bin_vol *= bin_width;
    }
    return ps_bin_vol;
}
// ==========================================================================
// Get Histograms, Not used in MAT implementation
// ==========================================================================

std::vector<TH2D *> PlotUtils::HyperDimLinearizer::Get2DHistos(PlotUtils::MnvH2D *result, bool IncludeSys = false) {
    //  std::cout <<"Entering Get2DHistos"  << std::endl;
    std::vector<TH2D *> expanded_results;
    if (m_analysis_type != k2D || m_analysis_type != k2D_lite) {  // This is only for 2D, so send it back if it's 1D type
        std::cout << "HyperDimLinearizer::Get2DHistos WARNING: you gave a 2D histogram, but have analysis type set to a 1D type. Returning blank list." << std::endl;
        return expanded_results;
    }

    TH2D mybigmap;
    if (!IncludeSys)
        mybigmap = result->GetCVHistoWithStatError();
    else
        mybigmap = result->GetCVHistoWithError();
    if (m_analysis_type == k2D) {
        //    std::cout << "Starting up get 2D histos with analysis type 0" << std::endl;
        // projected N dims (less Y) come in chunks of X bins (including under/over)
        int num_chunks = 1;
        const int num_x_bins = m_invec[0].size() - 1;
        const int num_y_bins = m_invec[1].size() - 1;
        //    std::cout << "Master Plot has x " << num_x_bins << "\ty\t" << num_y_bins << std::endl;
        for (unsigned int i = 2; i < m_invec.size(); i++)
            num_chunks *= (m_invec[i].size() + 1);
        //    std::cout << "I have number of chunks = " << num_chunks << std::endl;
        for (int i = 0; i < num_chunks; i++) {
            TH2D *tmp_bin = new TH2D(Form("Chunk_%d", i), Form("Chunk_%d", i), num_x_bins, &m_invec[0][0], num_y_bins, &m_invec[1][0]);
            int offset_x = (num_x_bins + 2) * i + 1;  // need to know low bin for chunk
            for (int j = 0; j < num_x_bins + 2; j++) {
                for (int k = 0; k < num_y_bins + 2; k++) {
                    double tmpval = mybigmap.GetBinContent(j + offset_x, k);
                    double tmperr = mybigmap.GetBinError(j + offset_x, k);
                    tmp_bin->SetBinContent(j, k, tmpval);
                    tmp_bin->SetBinError(j, k, tmperr);
                }                                                                     // end loop over subset y
            }                                                                         // end loop over subset x
            expanded_results.push_back((TH2D *)tmp_bin->Clone(Form("Clone_%d", i)));  // woohoo got a 2D result for one of the chunks!
        }                                                                             // end loop over chunks
    }
    return expanded_results;
}

std::vector<PlotUtils::MnvH2D *> PlotUtils::HyperDimLinearizer::Get2DMnvHistos(PlotUtils::MnvH2D *result, bool IncludeSys = false) {
    std::cout << "Entering Get2DMnvHistos" << std::endl;
    std::vector<PlotUtils::MnvH2D *> expanded_results;
    if (m_analysis_type != k2D || m_analysis_type != k2D_lite) {  // This is only for 2D, so send it back if it's 1D type
        std::cout << "HyperDimLinearizer::Get2DMnvHistos WARNING: you gave a 2D histogram, but have analysis type set to a 1D type. Returning blank list." << std::endl;
        return expanded_results;
    }

    std::vector<TH2D *> CV_vals = Get2DHistos(result, false);  // get CV
    std::cout << "I have " << CV_vals.size() << " CV histograms" << std::endl;
    for (uint i = 0; i < CV_vals.size(); i++)
        expanded_results.push_back(new PlotUtils::MnvH2D(*CV_vals[i]));

    std::vector<std::string> vertnames = result->GetVertErrorBandNames();
    std::vector<std::string> latnames = result->GetLatErrorBandNames();

    // Do vert first
    for (uint i = 0; i < vertnames.size(); i++) {
        std::cout << "Working on " << vertnames[i] << std::endl;
        std::vector<std::vector<TH2D *>> unihists;
        PlotUtils::MnvVertErrorBand2D *band = result->GetVertErrorBand(vertnames[i]);
        int bandsize = band->GetNHists();
        for (int uni = 0; uni < bandsize; uni++) {
            std::vector<TH2D *> tmpbandset = Get2DHistos(new PlotUtils::MnvH2D(*band->GetHist(uni)));  // Get the universe hist and spit out the N 2D results.
            unihists.push_back(tmpbandset);
        }
        // Have an N by uni matrix of TH2D. Now time to push back into the primary
        std::cout << "I have created a set of flux hists. This is size of the vector " << unihists.size() << "\t" << unihists[0].size() << std::endl;

        for (int j = 0; j < unihists[0].size(); j++) {  // unihists[0].size() is the number of projections needed
            std::vector<TH2D *> tmpband;
            for (int uni = 0; uni < bandsize; uni++) {
                tmpband.push_back(unihists[uni][j]);
            }
            expanded_results[j]->AddVertErrorBand(vertnames[i], tmpband);
        }
    }
    // now lat
    for (uint i = 0; i < latnames.size(); i++) {
        std::cout << "Working on " << latnames[i] << std::endl;
        std::vector<std::vector<TH2D *>> unihists;
        PlotUtils::MnvLatErrorBand2D *band = result->GetLatErrorBand(latnames[i]);
        int bandsize = band->GetNHists();
        for (int uni = 0; uni < bandsize; uni++) {
            std::vector<TH2D *> tmpbandset = Get2DHistos(new PlotUtils::MnvH2D(*band->GetHist(uni)));  // Get the universe hist and spit out the N 2D results.
            unihists.push_back(tmpbandset);
        }
        // Have an N by uni matrix of TH2D. Now time to push back into the primary

        std::cout << "I have created a set of flux hists. This is size of the vector " << unihists.size() << "\t" << unihists[0].size() << std::endl;

        for (uint j = 0; j < unihists[0].size(); j++) {  // unihists[0].size() is the number of projections needed
            std::vector<TH2D *> tmpband;
            for (int uni = 0; uni < bandsize; uni++) {
                tmpband.push_back(unihists[uni][j]);
            }
            expanded_results[j]->AddLatErrorBand(latnames[i], tmpband);
        }
    }
    return expanded_results;
}

std::vector<TH2D *> PlotUtils::HyperDimLinearizer::Get2DHistos(MnvH1D *result, bool IncludeSys = false) {
    std::vector<TH2D *> expanded_result;
    if (m_analysis_type == k2D || m_analysis_type == k2D_lite) { // This is only for 1D, so send it back if it's 2D type
        std::cout << "HyperDimLinearizer::Get2DHistos WARNING: you gave a 1D histogram, but have analysis type set to a 2D type. Returning blank list." << std::endl;
        return expanded_result;
    }

    TH1D result_hist;
    if (!IncludeSys)
        result_hist = result->GetCVHistoWithStatError();
    else
        result_hist = result->GetCVHistoWithError();

    const int n_x_bins = m_invec[0].size() - 1;
    const int n_y_bins = m_invec[1].size() - 1;

    int xy_cell_size;  // How big an x-y cell is
    if (m_cell_size.size() > 2)
        xy_cell_size = m_cell_size[2];  // If you have 3D or more linearized to 1D, x-y cell is the size of a z bin.
    else
        xy_cell_size = m_n_global_x_bins;  // If you only have 2D linearized to 1D, whole hist is one x-y cell

    const int n_xy_cells = m_n_global_x_bins / xy_cell_size;  // How many x-y cells are in your linearized histogram.
    int tmp_start = 1;                                       // Starting bin for each cell. Get's updated.
    for (int i = 0; i < n_xy_cells; i++) {                    // Loop over cells, making a 2D hist for each one
        TH2D *tmp_xy_cell = new TH2D(Form("XY_Cell_%d", i), Form("XY_Cell_%d", i), n_x_bins, &m_invec[0][0], n_y_bins, &m_invec[1][0]);
        for (int lin_bin = tmp_start; lin_bin < tmp_start + xy_cell_size; lin_bin++) { // Loop
            std::vector<int> ps_coords = GetValues(lin_bin);
            double tmp_val = result_hist.GetBinContent(lin_bin);
            double tmp_err = result_hist.GetBinError(lin_bin);
            tmp_xy_cell->SetBinContent(ps_coords[0], ps_coords[1], tmp_val);
            tmp_xy_cell->SetBinError(ps_coords[0], ps_coords[1], tmp_err);
        }
        expanded_result.push_back((TH2D *)tmp_xy_cell->Clone());
        tmp_start += xy_cell_size;
    }

    // for (int lin_bin = 1; lin_bin < m_n_global_x_bins + 1; lin_bin++) {
    //     std::vector<int> ps_coords = GetValues(lin_bin);
    // }
    return expanded_result;
}

// This is for 1D types
std::vector<PlotUtils::MnvH2D *> PlotUtils::HyperDimLinearizer::Get2DMnvHistos(PlotUtils::MnvH1D *result, bool IncludeSys = false) {
    std::cout << "Entering Get2DMnvHistos" << std::endl;
    std::vector<PlotUtils::MnvH2D *> expanded_result;
    if (m_analysis_type == k2D || m_analysis_type == k2D_lite) {  // This is only for 1D, so send it back if it's 2D type
        std::cout << "HyperDimLinearizer::Get2DMnvHistos WARNING: you gave a 1D histogram, but have analysis type set to a 2D type. Returning blank list." << std::endl;
        return expanded_result;
    }
    std::vector<TH2D *> CV_vals = Get2DHistos(result, false);  // get CV

    std::cout << "I have " << CV_vals.size() << " CV histograms" << std::endl;

    for (uint i = 0; i < CV_vals.size(); i++)
        expanded_result.push_back(new PlotUtils::MnvH2D(*CV_vals[i]));

    std::vector<std::string> vertnames = result->GetVertErrorBandNames();
    std::vector<std::string> latnames = result->GetLatErrorBandNames();

    // Do vert first
    for (uint i = 0; i < vertnames.size(); i++) {
        std::cout << "Working on " << vertnames[i] << std::endl;
        std::vector<std::vector<TH2D *>> unihists;
        PlotUtils::MnvVertErrorBand *band = result->GetVertErrorBand(vertnames[i]);
        int bandsize = band->GetNHists();
        for (int uni = 0; uni < bandsize; uni++) {
            std::vector<TH2D *> tmpbandset = Get2DHistos(new PlotUtils::MnvH1D(*band->GetHist(uni)));  // Get the universe hist and spit out the N 2D results.
            unihists.push_back(tmpbandset);
        }
        // Have an N by uni matrix of TH2D. Now time to push back into the primary
        std::cout << "I have created a set of flux hists. This is size of the vector " << unihists.size() << "\t" << unihists[0].size() << std::endl;

        for (int j = 0; j < unihists[0].size(); j++) {  // unihists[0].size() is the number of projections needed
            std::vector<TH2D *> tmpband;
            for (int uni = 0; uni < bandsize; uni++) {
                tmpband.push_back(unihists[uni][j]);
            }
            expanded_result[j]->AddVertErrorBand(vertnames[i], tmpband);
        }
    }
    // now lat
    for (uint i = 0; i < latnames.size(); i++) {
        std::cout << "Working on " << latnames[i] << std::endl;
        std::vector<std::vector<TH2D *>> unihists;
        PlotUtils::MnvLatErrorBand *band = result->GetLatErrorBand(latnames[i]);
        int bandsize = band->GetNHists();
        for (int uni = 0; uni < bandsize; uni++) {
            std::vector<TH2D *> tmpbandset = Get2DHistos(new PlotUtils::MnvH1D(*band->GetHist(uni)));  // Get the universe hist and spit out the N 2D results.
            unihists.push_back(tmpbandset);
        }
        // Have an N by uni matrix of TH2D. Now time to push back into the primary

        std::cout << "I have created a set of flux hists. This is size of the vector " << unihists.size() << "\t" << unihists[0].size() << std::endl;

        for (uint j = 0; j < unihists[0].size(); j++) {  // unihists[0].size() is the number of projections needed
            std::vector<TH2D *> tmpband;
            for (int uni = 0; uni < bandsize; uni++) {
                tmpband.push_back(unihists[uni][j]);
            }
            expanded_result[j]->AddLatErrorBand(latnames[i], tmpband);
        }
    }
    return expanded_result;
}

// This is for type 1 only, and only for 2D results
TH2D *PlotUtils::HyperDimLinearizer::Get2DHisto(PlotUtils::MnvH1D *result, bool IncludeSys = false) {
    if (m_invec.size() != 2)
        std::cout << "THIS ONLY WORKS FOR 2D RESULTS.\nIf you are a mapped 3D or more you need to use something different which might not exist." << std::endl;
    std::string myname = Form("Unmapped_%s", result->GetTitle());
    TH1D *mybigmap = NULL;
    if (!IncludeSys)
        mybigmap = new TH1D(result->GetCVHistoWithStatError());
    else
        mybigmap = new TH1D(result->GetCVHistoWithError());
    const int num_x_bins = m_invec[0].size() - 1;
    const int num_y_bins = m_invec[1].size() - 1;
    TH2D *my2D = new TH2D(myname.c_str(), myname.c_str(), num_x_bins, &m_invec[0][0], num_y_bins, &m_invec[1][0]);
    if (m_analysis_type == k1D) {
        std::cout << "Starting up get 2D histo with analysis type 1" << std::endl;
        for (int i = 0; i < num_x_bins + 2; i++)  // includes under/over
        {
            for (int j = 0; j < num_y_bins + 2; j++)  // includes under/over
            {
                double tmpval = mybigmap->GetBinContent(j * (num_x_bins + 2) + (i + 1));
                double tmperr = mybigmap->GetBinError(j * (num_x_bins + 2) + (i + 1));
                my2D->SetBinContent(i, j, tmpval);
                my2D->SetBinError(i, j, tmperr);
            }
        }
    }
    return my2D;
}

// This is for type 1 only, and only for 2D results
PlotUtils::MnvH2D *HyperDimLinearizer::Get2DMnvHisto(PlotUtils::MnvH1D *result, bool IncludeSys = false) {
    std::cout << "Entering Get2DMnvHisto" << std::endl;
    if (m_invec.size() != 2)
        std::cout << "THIS ONLY WORKS FOR 2D RESULTS.\nIf you are a mapped 3D or more you need to use something different which might not exist." << std::endl;

    TH2D *CV_vals = Get2DHisto(result, false);
    PlotUtils::MnvH2D *expanded_result = new PlotUtils::MnvH2D(*CV_vals);

    std::vector<std::string> vertnames = result->GetVertErrorBandNames();
    std::vector<std::string> latnames = result->GetLatErrorBandNames();

    for (uint i = 0; i < vertnames.size(); i++)  // Do vert first
    {
        std::cout << "Working on " << vertnames[i] << std::endl;
        std::vector<TH2D *> unihists;
        PlotUtils::MnvVertErrorBand *band = result->GetVertErrorBand(vertnames[i]);
        int bandsize = band->GetNHists();
        for (int uni = 0; uni < bandsize; uni++) {
            TH2D *tmpbandset = Get2DHisto(new PlotUtils::MnvH1D(*band->GetHist(uni)));  // Get the universe hist and spit out the N 2D results.
            unihists.push_back(tmpbandset);
        }
        expanded_result->AddVertErrorBand(vertnames[i], unihists);
    }

    for (uint i = 0; i < latnames.size(); i++) {
        std::cout << "Working on " << latnames[i] << std::endl;
        std::vector<TH2D *> unihists;
        PlotUtils::MnvLatErrorBand *band = result->GetLatErrorBand(latnames[i]);
        int bandsize = band->GetNHists();
        for (int uni = 0; uni < bandsize; uni++) {
            TH2D *tmpbandset = Get2DHisto(new PlotUtils::MnvH1D(*band->GetHist(uni)));  // Get the universe hist and spit out the N 2D results.
            unihists.push_back(tmpbandset);
        }
        expanded_result->AddLatErrorBand(latnames[i], unihists);
    }
    return expanded_result;
}

// This is for type 1 & 3 only, any dimension result, projected to an axis (default x)
std::vector<TH1D *> PlotUtils::HyperDimLinearizer::Get1DHistos(PlotUtils::MnvH1D *result, bool IncludeSys = false) {
    std::vector<TH1D *> expanded_result;
    if (m_analysis_type == k2D || m_analysis_type == k2D_lite) {  // This is only for 1D, so send it back if it's 2D type
        std::cout << "HyperDimLinearizer::Get2DHistos WARNING: you gave a 1D histogram, but have analysis type set to a 2D type. Returning blank list." << std::endl;
        return expanded_result;
    }

    TH1D result_hist;
    if (!IncludeSys)
        result_hist = result->GetCVHistoWithStatError();
    else
        result_hist = result->GetCVHistoWithError();

    const int n_x_bins = m_invec[0].size() - 1;
    const int n_y_bins = m_invec[1].size() - 1;

    const int x_cell_size = m_cell_size[1];
    const int n_x_cells = m_n_global_x_bins / m_cell_size[1];  // How many x cells are in your linearized histogram.
    int tmp_start = 1;                                         // Starting bin for each cell. Get's updated.
    for (int i = 0; i < n_x_cells; i++) {                     // Loop over cells, making a 2D hist for each one
        TH1D *tmp_x_cell = new TH1D(Form("X_Cell_%d", i), Form("X_Cell_%d", i), n_x_bins, &m_invec[0][0]);
        for (int lin_bin = tmp_start; lin_bin < tmp_start + x_cell_size; lin_bin++) {  // Loop
            std::vector<int> ps_coords = GetValues(lin_bin);
            double tmp_val = result_hist.GetBinContent(lin_bin);
            double tmp_err = result_hist.GetBinError(lin_bin);
            tmp_x_cell->SetBinContent(ps_coords[0], tmp_val);
            tmp_x_cell->SetBinError(ps_coords[0], tmp_err);
        }
        expanded_result.push_back((TH1D *)tmp_x_cell->Clone());
        tmp_start += x_cell_size;
    }
    return expanded_result;
}

std::vector<PlotUtils::MnvH1D *> PlotUtils::HyperDimLinearizer::Get1DMnvHistos(PlotUtils::MnvH1D *result, bool IncludeSys = false) {
    std::cout << "Entering Get1DMnvHistos" << std::endl;
    std::vector<PlotUtils::MnvH1D *> expanded_result;
    if (m_analysis_type == k2D || m_analysis_type == k2D_lite) {  // This is only for 1D, so send it back if it's 2D type
        std::cout << "HyperDimLinearizer::Get1DMnvHistos WARNING: you gave a 1D histogram, but have analysis type set to a 2D type. Returning blank list." << std::endl;
        return expanded_result;
    }
    std::vector<TH1D *> CV_vals = Get1DHistos(result, false);  // get CV

    std::cout << "I have " << CV_vals.size() << " CV histograms" << std::endl;

    for (uint i = 0; i < CV_vals.size(); i++)
        expanded_result.push_back(new PlotUtils::MnvH1D(*CV_vals[i]));

    std::vector<std::string> vertnames = result->GetVertErrorBandNames();
    std::vector<std::string> latnames = result->GetLatErrorBandNames();

    // Do vert first
    for (uint i = 0; i < vertnames.size(); i++) {
        std::cout << "Working on " << vertnames[i] << std::endl;
        std::vector<std::vector<TH1D *>> unihists;
        PlotUtils::MnvVertErrorBand *band = result->GetVertErrorBand(vertnames[i]);
        int bandsize = band->GetNHists();
        for (int uni = 0; uni < bandsize; uni++) {
            std::vector<TH1D *> tmpbandset = Get1DHistos(new PlotUtils::MnvH1D(*band->GetHist(uni)));  // Get the universe hist and spit out the N 2D results.
            unihists.push_back(tmpbandset);
        }
        // Have an N by uni matrix of TH2D. Now time to push back into the primary
        std::cout << "I have created a set of flux hists. This is size of the vector " << unihists.size() << "\t" << unihists[0].size() << std::endl;

        for (int j = 0; j < unihists[0].size(); j++) {  // unihists[0].size() is the number of projections needed
            std::vector<TH1D *> tmpband;
            for (int uni = 0; uni < bandsize; uni++) {
                tmpband.push_back(unihists[uni][j]);
            }
            expanded_result[j]->AddVertErrorBand(vertnames[i], tmpband);
        }
    }
    // now lat
    for (uint i = 0; i < latnames.size(); i++) {
        std::cout << "Working on " << latnames[i] << std::endl;
        std::vector<std::vector<TH1D *>> unihists;
        PlotUtils::MnvLatErrorBand *band = result->GetLatErrorBand(latnames[i]);
        int bandsize = band->GetNHists();
        for (int uni = 0; uni < bandsize; uni++) {
            std::vector<TH1D *> tmpbandset = Get1DHistos(new PlotUtils::MnvH1D(*band->GetHist(uni)));  // Get the universe hist and spit out the N 2D results.
            unihists.push_back(tmpbandset);
        }
        // Have an N by uni matrix of TH2D. Now time to push back into the primary

        std::cout << "I have created a set of flux hists. This is size of the vector " << unihists.size() << "\t" << unihists[0].size() << std::endl;

        for (uint j = 0; j < unihists[0].size(); j++) {  // unihists[0].size() is the number of projections needed
            std::vector<TH1D *> tmpband;
            for (int uni = 0; uni < bandsize; uni++) {
                tmpband.push_back(unihists[uni][j]);
            }
            expanded_result[j]->AddLatErrorBand(latnames[i], tmpband);
        }
    }
    return expanded_result;
}

// ==========================================================================
// Test
// ==========================================================================

void PlotUtils::HyperDimLinearizer::TestFunctionality() {
    std::cout << "Initializing funcationality test" << std::endl;

    // how many bins
    int n_bins = 1;
    std::cout << "Running with n-dimensions = " << m_invec.size() << std::endl;
    for (unsigned int i = 0; i < m_invec.size(); i++) {
        if (m_analysis_type == k2D && i == 1)
            continue;
        std::cout << "Coord " << i << " has " << m_el_size[i] << " bins " << std::endl;
        n_bins *= m_el_size[i];
    }
    std::cout << "This gives us a total of " << n_bins << " bins" << std::endl;

    for (int i = 0; i < n_bins; i++) {
        std::vector<int> coordinates = GetValues(i);
        std::cout << i << "\t";
        for (unsigned int j = 0; j < coordinates.size(); j++) {
            std::cout << coordinates[j] << "\t";
        }
        std::cout << std::endl;
    }
}

// These get you some values to check it worked as expected
PlotUtils::EAnalysisType PlotUtils::HyperDimLinearizer::GetAnalysisType() {
    return m_analysis_type;
}

std::vector<int> PlotUtils::HyperDimLinearizer::GetAxesSizes() {
    return m_el_size;
}

std::vector<int> PlotUtils::HyperDimLinearizer::GetCellSizes() {
    return m_cell_size;
}

int PlotUtils::HyperDimLinearizer::GetNLinBins() {
    if (m_analysis_type == k2D_lite || m_analysis_type == k1D_lite) {  // If doing lightweight, you'll also have under/overflow out at the end.
        // std::cout << "Lite analysis type, returning (# of global bins) + (# of under/overflow bins)" << std::endl;
        int flow_bins = 1;
        for (unsigned int i = 0; i < m_invec.size(); i++) {
            if (m_analysis_type == k2D_lite && i == 1)
                continue;
            flow_bins *= 3;
        }
        return m_n_global_x_bins + flow_bins;
    }

    return m_n_global_x_bins;
}

// ==========================================================================
// Helpers
// ==========================================================================

bool PlotUtils::HyperDimLinearizer::IsUnderflow(int lin_bin, int axis = -1) {
    if (m_analysis_type != 1) {
        std::cout << "WARNING: PlotUtils::HyperDimLinearizer::IsUnderflow is only configured for type 1 analyses" << std::endl;
        return false;
    }
    if (m_analysis_type == k2D_lite || m_analysis_type == k1D_lite) {
        std::cout << "Analysis type 2 or 3, so no underflow bins in normal bins." << std::endl;
        return false;
    }
    // Given a bin in linearized space, check if it's an underflow bin in phase space
    std::vector<int> ps_coords = GetValues(lin_bin);
    // For a given axis...
    if (axis >= 0) {
        if (ps_coords[axis] == 0) {
            return true;
        }
    }
    // In general...
    else {
        for (int i = 0; i < ps_coords.size(); i++) {
            if (ps_coords[i] == 0) {
                return true;
            }
        }
    }
    return false;
}

bool PlotUtils::HyperDimLinearizer::IsOverflow(int lin_bin, int axis = -1) {
    if (m_analysis_type != 1) {
        std::cout << "WARNING: HyperDimLinearizer::IsOverflow is only configured for type 1 analyses" << std::endl;
        return false;
    }
    if (m_analysis_type == k2D_lite || m_analysis_type == k1D_lite) {
        std::cout << "Analysis type 2 or 3, so no overflow bins in normal bins." << std::endl;
        return false;
    }
    // Given a bin in linearized space, check if it's an overflow bin in phase space
    std::vector<int> ps_coords = GetValues(lin_bin);
    // For a given axis...
    if (axis >= 0) {
        if (ps_coords[axis] == m_el_size[axis]) {
            return true;
        }
    }
    // In general...
    else {
        for (int i = 0; i < ps_coords.size(); i++) {
            if (ps_coords[i] == m_el_size[i]) {
                return true;
            }
        }
    }
    return false;
}

}  // namespace PlotUtils
#endif
