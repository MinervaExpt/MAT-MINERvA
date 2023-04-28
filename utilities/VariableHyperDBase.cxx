#ifndef VARIABLEHYPERDBASE_CXX
#define VARIABLEHYPERDBASE_CXX

#include "VariableHyperDBase.h"

#include "utilities/HyperDimLinearizer.h"
// #include "PlotUtils/HyperDimLinearizer.h"

using namespace PlotUtils;

//==============================================================================
//
// Variable HyperD Base
//
//==============================================================================
//==============================================================================
// CTORS
//==============================================================================
// Default. Recommended if using a derived Variable class
template <class UNIVERSE>
VariableHyperDBase<UNIVERSE>::VariableHyperDBase(const EAnalysisType t2D_t1D) {
    m_analysis_type = t2D_t1D;  // default is 1D
                                // You will need to add variables, etc.
}

// Default but set the name too. Recommended if using a derived Variable class
template <class UNIVERSE>
VariableHyperDBase<UNIVERSE>::VariableHyperDBase(std::string name, const EAnalysisType t2D_t1D) {
    m_name = name;
    m_analysis_type = t2D_t1D;  // default is 1D
                                // You will need to add variables, etc.
}

// Constructor with vector of input variables, will likely have issues if using derived Variable class
template <class UNIVERSE>
VariableHyperDBase<UNIVERSE>::VariableHyperDBase(const std::vector<VariableBase<UNIVERSE> *> &d, const EAnalysisType t2D_t1D) {
    m_vars_vec = d;
    m_analysis_type = t2D_t1D;
    for (int i = 0; i < d.size(); i++) {
        if (d[i]->HasRecoBinning()) m_has_reco_binning = true;
    }

    Setup();
}

// Constructor with user designated name & vector of input variables, will likely have issues if using derived Variable class
template <class UNIVERSE>
VariableHyperDBase<UNIVERSE>::VariableHyperDBase(std::string name, const std::vector<VariableBase<UNIVERSE> *> &d, const EAnalysisType t2D_t1D) {
    m_name = name;
    m_vars_vec = d;
    m_analysis_type = t2D_t1D;
    for (int i = 0; i < d.size(); i++) {
        if (d[i]->HasRecoBinning()) m_has_reco_binning = true;
    }

    Setup(name);
}

//==============================================================================
// Set
//==============================================================================

// Change the name
template <class UNIVERSE>
std::string VariableHyperDBase<UNIVERSE>::SetName(const std::string name) {
    return m_name = name;
}

// Change the analysis type. Right now it is fixed to 1D, or "type 1"
template <class UNIVERSE>
void VariableHyperDBase<UNIVERSE>::SetAnalysisType(const EAnalysisType t2D_t1D) {
    m_analysis_type = t2D_t1D;
    std::cout << "VariableHyperDBase: WARNING you may be changing your analysis type to " << m_analysis_type << std::endl;

    // Reset everything given new analysis type.
    Setup(GetName());
}

// Add another variable, useful if using default constructor and derived Variable class
template <class UNIVERSE>
void VariableHyperDBase<UNIVERSE>::AddVariable(VariableBase<UNIVERSE> &var) {
    // If using derived, you'll need to get all the info you want out beforehand.
    m_vars_vec.emplace_back(new VariableBase<UNIVERSE>(var));

    if (var.HasRecoBinning()) m_has_reco_binning = true;

    // Reset everything given new vars.
    Setup();
}

// Setup variable. This is private and should only be used internally for now.
template <class UNIVERSE>
void VariableHyperDBase<UNIVERSE>::Setup(const std::string i_name)  // i_name defaults to "" and it will build a name for you instead
{
    m_dimension = m_vars_vec.size();

    std::string name;
    std::string lin_axis_label;

    std::vector<std::vector<double>> vars_bins;       // List of the bin edges lists for each variable to feed into hyperdim
    int n_lin_bins = 1;                               // Number of linearized bins
    std::vector<std::vector<double>> vars_reco_bins;  // Same for reco if you have separate reco binning
    int n_lin_reco_bins = 1;

    for (int i = 0; i < m_dimension; i++) {
        // Make vector of binnings, count number of bins
        std::vector<double> var_binning = m_vars_vec[i]->GetBinVec();
        vars_bins.push_back(var_binning);

        int tmp_bin_scale = 1;

        if (m_analysis_type == k1D_lite || m_analysis_type == k2D_lite)
            tmp_bin_scale = var_binning.size() - 1;  // Bin edges - 1 gives number of bins, exculding under/overflow
        else                                         // if k2D and k1D
            tmp_bin_scale = var_binning.size() + 1;  // Bin edges - 1 + 2, including under/overflow

        if (i == 1 && (m_analysis_type == k2D || m_analysis_type == k2D_lite))  // if doing 2D projection, skip y-axis
            n_lin_bins *= 1;
        else
            n_lin_bins *= tmp_bin_scale;

        // If recobinning is set (should happen at initialization or setup), do this too
        if (m_has_reco_binning) {
            // Make vector of binnings, count number of bins
            std::vector<double> var_reco_binning = m_vars_vec[i]->GetRecoBinVec();
            vars_reco_bins.push_back(var_reco_binning);
            int tmp_reco_bin_scale = 1;

            if (m_analysis_type == k1D_lite || m_analysis_type == k2D_lite)
                tmp_reco_bin_scale = var_reco_binning.size() - 1;  // Bin edges - 1 gives number of bins, exculding under/overflow
            else                                                   // if k2D and k1D
                tmp_reco_bin_scale = var_reco_binning.size() + 1;  // Bin edges - 1 + 2, including under/overflow

            if (i == 1 && (m_analysis_type == k2D || m_analysis_type == k2D_lite))  // if doing 2D projection, skip y-axis
                n_lin_reco_bins *= 1;
            else
                n_lin_reco_bins *= tmp_reco_bin_scale;
        }

        // Make the name and axis label. Won't use name if input name is defined.
        name += m_vars_vec[i]->GetName();
        if (i < (m_dimension - 1))
            name += "_";
        if (m_analysis_type == k2D && i == 1)  // if doing 2D, set y axis label and don't add it to the linearized axis label
        {
            m_y_axis_label = m_vars_vec[i]->GetAxisLabel();
            continue;
        }
        lin_axis_label += m_vars_vec[i]->GetAxisLabel();
        if (i < (m_dimension - 1))
            lin_axis_label += ", ";
    }

    // If there's an input name, it will use that one, otherwise it'll use the one it just made
    if (i_name.size() > 0)
        m_name = i_name;
    else
        m_name = name;
    // Store axis label, list of input binnings
    m_lin_axis_label = lin_axis_label;
    m_vars_binnings = vars_bins;

    // Initialize a hyperdim with that binning
    m_hyperdim = new HyperDimLinearizer(vars_bins, m_analysis_type);

    // Make a bin vector in linearized bin index space
    std::vector<double> lin_binning;
    // TODO: Does this need to start at 1 instead?
    for (int i = 0; i < n_lin_bins + 1; i++) lin_binning.push_back(i);
    m_lin_binning = lin_binning;

    // If you have reco binning, do that here too
    if (m_has_reco_binning) {
        m_vars_reco_binnings = vars_reco_bins;
        m_reco_hyperdim = new HyperDimLinearizer(vars_reco_bins, m_analysis_type);

        // TODO: Does this need to start at 1 instead?
        std::vector<double> lin_reco_binning;
        for (int i = 0; i < n_lin_reco_bins + 1; i++) lin_reco_binning.push_back(i);
        m_lin_reco_binning = lin_reco_binning;
    }
}

//==============================================================================
// Gets (truth binning, reco further down)
//==============================================================================

// Get name of linearized variable
template <class UNIVERSE>
std::string VariableHyperDBase<UNIVERSE>::GetName() const {
    return m_name;
}

// Get name of a component variable
template <class UNIVERSE>
std::string VariableHyperDBase<UNIVERSE>::GetName(int axis) const {
    return m_vars_vec[axis]->GetName();
}

// Get axis label of linearized variable
template <class UNIVERSE>
std::string VariableHyperDBase<UNIVERSE>::GetAxisLabel() const {
    return m_lin_axis_label;
}

// Get axis label of a component variable
template <class UNIVERSE>
std::string VariableHyperDBase<UNIVERSE>::GetAxisLabel(int axis) const {
    return m_vars_vec[axis]->GetAxisLabel();
}

// Get a vector of the axis labels of all component variables
template <class UNIVERSE>
std::vector<std::string> VariableHyperDBase<UNIVERSE>::GetAxisLabelVec() const {
    std::vector<std::string> axis_label_vec = {};
    for (const auto var : m_vars_vec) {
        axis_label_vec.push_back(var->GetAxisLabel());
    }
    return axis_label_vec;
}

// Get the number of linearized bins
template <class UNIVERSE>
int VariableHyperDBase<UNIVERSE>::GetNBins() const {
    return m_lin_binning.size() - 1;
    // Linearized binning has under/overflow bins in phase space spread throughout, so this number may seem bigger than expected.
}

// Get the number of bins in a given axis
template <class UNIVERSE>
int VariableHyperDBase<UNIVERSE>::GetNBins(int axis) const {
    return m_vars_vec[axis]->GetNBins();
}

// Get the linearized bin edges
template <class UNIVERSE>
std::vector<double> VariableHyperDBase<UNIVERSE>::GetBinVec() const {
    return m_lin_binning;
}

// Get the bin edges of a given axis
template <class UNIVERSE>
std::vector<double> VariableHyperDBase<UNIVERSE>::GetBinVec(int axis) const {
    return m_vars_vec[axis]->GetBinVec();
}

// Print linearized binning, should just be list of indexes
template <class UNIVERSE>
void VariableHyperDBase<UNIVERSE>::PrintBinning() const {
    std::cout << GetName() << " binning: ";
    for (const auto b : m_lin_binning) std::cout << b << " ";
    std::cout << "\n";

    if (m_analysis_type == k2D) {
        std::cout << " y axis: ";
        for (const auto b : m_vars_vec[1]->GetBinVec()) {
            std::cout << b << " ";
        }
        std::cout << "\n";
    }
}

// Print the binning of an axis
template <class UNIVERSE>
void VariableHyperDBase<UNIVERSE>::PrintBinning(int axis) const {
    m_vars_vec[axis]->PrintBinning();
}

// Get the hyperdim
template <class UNIVERSE>
PlotUtils::HyperDimLinearizer *VariableHyperDBase<UNIVERSE>::GetHyperDimLinearizer() const {
    return m_hyperdim;
}

// Get the phase space volume for a given linearized bin. If doing a 2D (type 0) analysis, this does not include the bin width of y axis.
// TODO: Maybe this belongs to hyperdim?
template <class UNIVERSE>
double VariableHyperDBase<UNIVERSE>::GetBinVolume(int lin_bin) const {
    double ps_bin_vol = 1.;
    // Given the linearized bin number, get corresponding bin index in phase space coordinates
    std::vector<int> ps_coords = m_hyperdim->GetValues(lin_bin);  // This automatically skips y axis for 2D
    for (int i = 0; i < ps_coords.size(); i++) {
        int var_bin = ps_coords[i];
        // Get the binning for that axis. Need to do some trickery to index properly if doing 2D.
        int var = i;
        if (m_analysis_type == k2D && i >= 1)
            var += 1;
        std::vector<double> var_binning = m_vars_vec[var]->GetBinVec();
        double bin_width = var_binning[var_bin] - var_binning[var_bin - 1];
        if (bin_width < 0)
            std::cout << "VariableHyperDBase: WARNING, negative binwidth, will affect bin volume" << std::endl;
        else if (bin_width == 0)
            std::cout << "VariableHyperDBase: WARNING, 0 binwidth, will affect bin volume" << std::endl;
        ps_bin_vol *= bin_width;
    }
    return ps_bin_vol;
}

template <class UNIVERSE>
std::vector<int> PlotUtils::VariableHyperDBase<UNIVERSE>::GetUnderflow(int axis) const {
    // TODO
    std::cout << "WARNING: VariableHyperDBase::GetUnderflow() not set up yet... " << std::endl;
    return std::vector<int>();
}

template <class UNIVERSE>
std::vector<int> PlotUtils::VariableHyperDBase<UNIVERSE>::GetOverflow(int axis) const {
    // TODO
    std::cout << "WARNING: VariableHyperDBase::GetOverflow() not set up yet... " << std::endl;
    return std::vector<int>();
}

//==============================================================================
// Reco Getters
//==============================================================================

// Check if there's reco binning set. If set to false reco getters will return true things
template <class UNIVERSE>
bool VariableHyperDBase<UNIVERSE>::HasRecoBinning() const {
    return m_has_reco_binning;
}

template <class UNIVERSE>
bool VariableHyperDBase<UNIVERSE>::HasRecoBinning(int axis) const {
    return m_vars_vec[axis]->HasRecoBinning();
}

template <class UNIVERSE>
int VariableHyperDBase<UNIVERSE>::GetNRecoBins() const {
    if (!m_has_reco_binning) return m_lin_binning.size() - 1;

    return m_lin_reco_binning.size() - 1;
    // Linearized binning has under/overflow bins in phase space spread throughout, so this number may seem bigger than expected.
}

template <class UNIVERSE>
int VariableHyperDBase<UNIVERSE>::GetNRecoBins(int axis) const {
    return m_vars_vec[axis]->GetNRecoBins();  // If variable has no reco binning, will return truth instead
}

// Get the linearized reco bin edges
template <class UNIVERSE>
std::vector<double> VariableHyperDBase<UNIVERSE>::GetRecoBinVec() const {
    if (!m_has_reco_binning) return m_lin_binning;

    return m_lin_reco_binning;
}

// Get the reco bin edges of a given axis
template <class UNIVERSE>
std::vector<double> VariableHyperDBase<UNIVERSE>::GetRecoBinVec(int axis) const {
    return m_vars_vec[axis]->GetRecoBinVec();  // If variable has no reco binning, will return truth instead
}

// Print linearized reco binning, should just be list of indexes
template <class UNIVERSE>
void VariableHyperDBase<UNIVERSE>::PrintRecoBinning() const {
    if (!m_has_reco_binning) {
        std::cout << GetName() << " reco binning same as true" << std::endl;
        PrintBinning();
    } else {
        std::cout << GetName() << " reco binning: ";
        for (const auto b : GetRecoBinVec()) {
            std::cout << b << " ";
        }
        std::cout << "\n";

        if (m_analysis_type == k2D) {
            std::cout << " y axis: ";
            for (const auto b : m_vars_vec[1]->GetRecBinVec()) {
                std::cout << b << " ";
            }
            std::cout << "\n";
        }
    }
}

// Print the reco binning of an axis
template <class UNIVERSE>
void VariableHyperDBase<UNIVERSE>::PrintRecoBinning(int axis) const {
    m_vars_vec[axis]->PrintRecoBinning();  // If variable has no reco binning, will return truth instead
}

// Get the reco hyperdim
template <class UNIVERSE>
PlotUtils::HyperDimLinearizer *VariableHyperDBase<UNIVERSE>::GetRecoHyperDimLinearizer() const {
    if (!m_has_reco_binning) {
        return m_hyperdim;
    }
    return m_reco_hyperdim;
}

// Get the phase space volume for a given linearized reco bin. If doing a 2D (type 0) analysis, this does not include the bin width of y axis.
// TODO: Maybe this belongs to hyperdim?
template <class UNIVERSE>
double VariableHyperDBase<UNIVERSE>::GetRecoBinVolume(int lin_bin) const {
    if (!m_has_reco_binning) {  // If there's no reco binning, just return truth instead
        return GetBinVolume(lin_bin);
    }
    // Given the linearized bin number, get corresponding bin number in phase space coordinates
    double ps_bin_vol = 1.;
    std::vector<int> ps_coords = m_reco_hyperdim->GetValues(lin_bin);  // This automatically skips y axis for 2D

    for (int i = 0; i < ps_coords.size(); i++) {
        // Get the binning for that axis. Need to do some trickery to index properly if doing 2D.
        int var = i;
        if (m_analysis_type == k2D && i >= 1)
            var += 1;
        int var_bin = ps_coords[i];
        std::vector<double> var_binning = m_vars_vec[var]->GetRecoBinVec();
        double bin_width = var_binning[var_bin] - var_binning[var_bin - 1];

        if (bin_width < 0)
            std::cout << "VariableHyperDBase: WARNING, negative binwidth, will affect bin volume" << std::endl;
        else if (bin_width == 0)
            std::cout << "VariableHyperDBase: WARNING, 0 binwidth, will affect bin volume" << std::endl;

        ps_bin_vol *= bin_width;
    }
    return ps_bin_vol;
}

template <class UNIVERSE>
std::vector<int> PlotUtils::VariableHyperDBase<UNIVERSE>::GetRecoUnderflow(int axis) const {
    if (m_analysis_type > 1) {
        std::cout << "VariableHyperDBase: analysis type selected uses global underflow bins. " << std::endl;
        return std::vector<int>();
    }
    // TODO: Get a vector of indices in of bins linearized space that are underflow bins in phase space
    std::cout << "WARNING: VariableHyperDBase::GetRecoUnderflow() not set up yet... " << std::endl;
    return std::vector<int>();
}

template <class UNIVERSE>
std::vector<int> PlotUtils::VariableHyperDBase<UNIVERSE>::GetRecoOverflow(int axis) const {
    if (m_analysis_type > 1) {
        std::cout << "VariableHyperDBase: analysis type selected uses global overflow bins. " << std::endl;
        return std::vector<int>();
    }
    // TODO: Get a vector of indices in of bins linearized space that are overflow bins in phase space
    std::cout << "WARNING: VariableHyperDBase::GetRecoOverflow() not set up yet... " << std::endl;
    return std::vector<int>();
}

//==============================================================================
// GetValues
//==============================================================================

// Return bin index value in linearized bin space, ie which bin an event goes in
template <class UNIVERSE>
double VariableHyperDBase<UNIVERSE>::GetRecoValue(const UNIVERSE &universe,
                                                  const int idx1,
                                                  const int idx2) const {
    std::vector<double> val_vec;
    for (int i = 0; i < m_dimension; i++) {
        val_vec.push_back(m_vars_vec[i]->GetRecoValue(universe, idx1, idx2));
    }
    if (!m_has_reco_binning) {
        return (m_hyperdim->GetBin(val_vec).first) + 0.0001;  // 0.0001 offset to so value isn't exactly on a bin edge and fillers can put it in that bin
    } else {
        return (m_reco_hyperdim->GetBin(val_vec).first) + 0.0001;  // If there's reco binning, use that hyperdim
    }
}

// Return bin index value in linearized bin space, ie which an bin event goes in
template <class UNIVERSE>
double VariableHyperDBase<UNIVERSE>::GetTrueValue(const UNIVERSE &universe,
                                                  const int idx1,
                                                  const int idx2) const {
    std::vector<double> val_vec;
    for (int i = 0; i < m_dimension; i++) {
        val_vec.push_back(m_vars_vec[i]->GetTrueValue(universe, idx1, idx2));
    }
    return (m_hyperdim->GetBin(val_vec).first) + 0.0001;  // 0.0001 offset to so fillers can put it in that bin
}

// Return phase space value of a single component variable, use axis = 1 to get y-value for "type 0" 2D analyses.
template <class UNIVERSE>
double VariableHyperDBase<UNIVERSE>::GetRecoValue(const int axis,
                                                  const UNIVERSE &universe,
                                                  const int idx1,
                                                  const int idx2) const {
    return m_vars_vec[axis]->GetRecoValue(universe, idx1, idx2);
}

// Return phase space value of a single component variable, use axis = 1 to get y-value for "type 0" 2D analyses.
template <class UNIVERSE>
double VariableHyperDBase<UNIVERSE>::GetTrueValue(const int axis,
                                                  const UNIVERSE &universe,
                                                  const int idx1,
                                                  const int idx2) const {
    return m_vars_vec[axis]->GetTrueValue(universe, idx1, idx2);
}

// Return phase space values of every component variable as a vector
template <class UNIVERSE>
std::vector<double> VariableHyperDBase<UNIVERSE>::GetRecoValueVec(const UNIVERSE &universe,
                                                                  const int idx1,
                                                                  const int idx2) const {
    std::vector<double> value_vec;
    for (int i = 0; i < m_dimension; i++) {
        value_vec.push_back(m_vars_vec[i]->GetRecoValue(universe, idx1, idx2));
    }
    return value_vec;
}

// Return phase space values of every component variable as a vector
template <class UNIVERSE>
std::vector<double> VariableHyperDBase<UNIVERSE>::GetTrueValueVec(const UNIVERSE &universe,
                                                                  const int idx1,
                                                                  const int idx2) const {
    std::vector<double> value_vec;
    for (int i = 0; i < m_dimension; i++) {
        value_vec.push_back(m_vars_vec[i]->GetTrueValue(universe, idx1, idx2));
    }
    return value_vec;
}

#endif  // VARIABLEHYPERDBASE_CXX
