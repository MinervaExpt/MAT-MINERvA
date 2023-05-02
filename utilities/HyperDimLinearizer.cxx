#ifndef MNV_HYPERDIMLINEARIZER_cxx
#define MNV_HYPERDIMLINEARIZER_cxx 1
#include "utilities/HyperDimLinearizer.h"

#include <string>

namespace PlotUtils {

//! Get a set of bins
// ==========================================================================
// CTOR
// ==========================================================================
HyperDimLinearizer::HyperDimLinearizer(std::vector<std::vector<double>> input, int type) {
    if (type == 0) {
        m_analysis_type = k2D;
    } else if (type == 1) {
        m_analysis_type = k1D;
    } else if (type == 2) {
        m_analysis_type = k2D_lite;
    } else if (type == 3) {
        m_analysis_type = k1D_lite;
    }

    // Jank to make things build right for now. TODO: Make less jank
    HyperDimLinearizer(input, m_analysis_type);
}

HyperDimLinearizer::HyperDimLinearizer(std::vector<std::vector<double>> input, EAnalysisType type) {
    m_invec = input;
    int n_dim = input.size();
    m_analysis_type = type;

    std::cout << "Contructing class with " << n_dim << " dimensions of type " << type << std::endl;
    if (type == k2D_lite || type == k1D_lite)
        std::cout << " Selected lightweight type. Under/overflow bins will not be counted." << std::endl;

    m_cell_size = {1};      // vector to put sizes of each cell for each bin in each axis. Should start with x-axis which is always 1 bin.
    m_n_global_x_bins = 1;  // will tell you how many bins you have on your linearized axis

    for (unsigned int i = 0; i < input.size(); i++) {
        // Figure out size of each axis in phase space depending on analysis type
        int tmp_el_size;
        if (type == k2D || type == k1D) {  // number of bins = (vector size - 1 + 2) counting under/overflow for each axis
            tmp_el_size = input[i].size() + 1;
        } else {  // if (type == k2D_lite || type == k1D_lite) {  // number of bins = (vector size - 1) using global under/overflow
            tmp_el_size = input[i].size() - 1;
        }
        std::cout << "Bin number " << i << "\t" << tmp_el_size << std::endl;
        m_el_size.push_back(tmp_el_size);

        // Count the number of bins in linearized space and size of cells for each axis
        if (type == k1D || type == k1D_lite) {  // for type 1 and 3 fully linearized, this is simple
            m_n_global_x_bins *= tmp_el_size;
            if (i < input.size() - 1)  // skip last element since we already start with x-axis cell size
                m_cell_size.push_back(tmp_el_size * m_cell_size[i]);
        } else {         // Type 0 and 2 keep y axis so need to do some kerjiggering
            if (i != 1)  // For counting number of bins, skip y axis if doing type 0, 2
                m_n_global_x_bins *= tmp_el_size;
            if (i == 0) {
                m_cell_size.push_back(0);                             // y-axis should have cell size of 0 in type 0 and 2
                m_cell_size.push_back(tmp_el_size * m_cell_size[i]);  // This takes care of z-axis, so we skip i = 2 to avoid duplicate
            } else if (i > 1 && i < input.size() - 1) {               // When looking at higher than 3D
                m_cell_size.push_back(tmp_el_size * m_cell_size[i]);
            }
        }
    }
}

// ==========================================================================
// Getter values between spaces
// ==========================================================================

// Transform values in phase space (in units of each axis) to linearized bin values.
//     Return ranges from [0, ((# lin bins) - 1)] inclusively
//     For lite types, additional under/overflow bins placed at high end of range from [(# lin bins), ((# lin bins) + (3^n_dim) - 1)] inclusively
//     For 2D types, y axis is binned normally; with 0 being underflow, ((# y bins) + 1) being overflow. Return 0 for y value on 1D types.
std::pair<int, int> HyperDimLinearizer::GetBin(std::vector<double> values) {
    int global_x = 0;  // Returned linearized bin
    int y_bin = 0;

    if (m_analysis_type == k2D || m_analysis_type == k1D) {  // These include under/overflow bins
        for (unsigned int i = 0; i < values.size(); i++) {
            int tmp_bin = Get1DBin(values[i], i);  // Find the bin index on a given axis given a value in that axis
            global_x += tmp_bin * m_cell_size[i];  // Add that many cells of that axis to global_x in linearized space, if doing 2D, y should have cell size of 0
        }
    } else if (m_analysis_type == k2D_lite || m_analysis_type == k1D_lite) {  // These use more global under/overflow bins out at the end
        int tmp_global_x = 0;                                                 // Placeholder
        bool underover_bool = false;                                          // Switch turns on if you are in under/overflow on any axis (except y for 2D)
        int flow_x = 0;                                                       // Bin in the under/overflow bins at end of linearized x axis
        int flow_cell = 1;                                                    // Cell size for over/underflow bins (gets changed in the loop)

        for (unsigned int i = 0; i < values.size(); i++) {
            if (i == 1 && m_analysis_type == k2D_lite)  // Skip y axis if doing 2D, necessary to avoid putting y under/overflow in the global_x over/under
                continue;
            int tmp_bin = Get1DBin(values[i], i);
            if (tmp_bin == 0) {  // If underflow
                // flow_x += 0;                                 // (# of flow cells - 1) * (flow cell size), but # of flow cells is always 0 here
                underover_bool = true;
            } else if (tmp_bin == m_el_size[i] + 1) {  // If overflow
                flow_x += 2 * flow_cell;               // (# of flow cells - 1) * (flow cell size)
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
int HyperDimLinearizer::Get1DBin(double value, int el) {
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
std::vector<int> HyperDimLinearizer::GetValues(int x_linbin, int y_bin) {  //  Default ybin to 0 to maintain behaviour from Dan's version
    std::vector<int> ps_bin_coords;                                        // Phase space bin coordinates, returned here

    if (x_linbin >= m_n_global_x_bins) {  // TODO: what should you do if under/overflow for lite types cases?
        std::cout << "HyperDimLinearizer::GetValues: WARNING: requested values out of range of linearized space. Result might not make sense." << std::endl;
    }

    int mod_bin = x_linbin;                                                       // Place holder
    for (unsigned int i = 0; i < m_invec.size(); i++) {                           // Loop over coordinate axes
        if ((m_analysis_type == k2D || m_analysis_type == k2D_lite) && i == 1) {  // If doing 2D, put in y-bin and skip the other steps for y.
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

// ==========================================================================
// Get Histograms, Not used in MAT implementation
// ==========================================================================

// TODO: Make sure this all works for type 2,3
std::vector<TH2D *> HyperDimLinearizer::Get2DHistos(PlotUtils::MnvH2D *result, bool IncludeSys = false) {
    //  std::cout <<"Entering Get2DHistos"  << std::endl;
    std::vector<TH2D *> expanded_results;
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

// TODO: Make sure this works for type 2,3
std::vector<PlotUtils::MnvH2D *> HyperDimLinearizer::Get2DMnvHistos(PlotUtils::MnvH2D *result, bool IncludeSys = false) {
    std::cout << "Entering Get2DMnvHistos" << std::endl;
    std::vector<PlotUtils::MnvH2D *> expanded_results;
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

TH2D *HyperDimLinearizer::Get2DHisto(PlotUtils::MnvH1D *result, bool IncludeSys = false) {
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

// TODO: Make sure this works for type 2,3
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

// ==========================================================================
// Test
// ==========================================================================

void HyperDimLinearizer::TestFunctionality() {
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
EAnalysisType HyperDimLinearizer::GetAnalysisType() {
    return m_analysis_type;
}

std::vector<int> HyperDimLinearizer::GetAxesSizes() {
    return m_el_size;
}

std::vector<int> HyperDimLinearizer::GetCellSizes() {
    return m_cell_size;
}

int HyperDimLinearizer::GetNLinBins() {
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

bool HyperDimLinearizer::IsUnderflow(int lin_bin, int axis = -1) {
    if (m_analysis_type != 1) {
        std::cout << "WARNING: HyperDimLinearizer::IsUnderflow is only configured for type 1 analyses" << std::endl;
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

bool HyperDimLinearizer::IsOverflow(int lin_bin, int axis = -1) {
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
