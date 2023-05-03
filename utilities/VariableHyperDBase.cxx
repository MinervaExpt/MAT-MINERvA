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
VariableHyperDBase<UNIVERSE>::VariableHyperDBase(const EAnalysisType type) {
    m_analysis_type = type;  // Set to k2D, k1D, k2D_lite, k1D_lite
                             // You will need to add variables, etc.
}

// Default but set the name too. Recommended if using a derived Variable class
template <class UNIVERSE>
VariableHyperDBase<UNIVERSE>::VariableHyperDBase(std::string name, const EAnalysisType type) {
    m_name = name;
    m_analysis_type = type;  // Set to k2D, k1D, k2D_lite, k1D_lite
                             // You will need to add variables, etc.
}

// Constructor with vector of input variables, will likely have issues if using derived Variable class
template <class UNIVERSE>
VariableHyperDBase<UNIVERSE>::VariableHyperDBase(const std::vector<VariableBase<UNIVERSE> *> &d, const EAnalysisType type) {
    m_vars_vec = d;
    m_analysis_type = type;
    for (int i = 0; i < d.size(); i++) {
        if (d[i]->HasRecoBinning()) m_has_reco_binning = true;
    }

    Setup();
}

// Constructor with user designated name & vector of input variables, will likely have issues if using derived Variable class
template <class UNIVERSE>
VariableHyperDBase<UNIVERSE>::VariableHyperDBase(std::string name, const std::vector<VariableBase<UNIVERSE> *> &d, const EAnalysisType type) {
    m_name = name;
    m_vars_vec = d;
    m_analysis_type = type;
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

// Change the analysis type. Set to k2D, k1D, k2D_lite, k1D_lite
template <class UNIVERSE>
void VariableHyperDBase<UNIVERSE>::SetAnalysisType(const EAnalysisType type) {
    m_analysis_type = type;
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
void VariableHyperDBase<UNIVERSE>::Setup(const std::string i_name) {  // i_name defaults to "" and it will build a name for you instead
    m_dimension = m_vars_vec.size();

    std::string name;
    std::string lin_axis_label;

    std::vector<std::vector<double>> vars_bins;       // List of the bin edges lists for each variable to feed into hyperdim
    std::vector<std::vector<double>> vars_reco_bins;  // Same for reco if you have separate reco binning

    for (int i = 0; i < m_dimension; i++) {
        // Make vector of binnings, count number of bins
        vars_bins.push_back(m_vars_vec[i]->GetBinVec());

        // If recobinning is set (should happen at initialization or setup), do this for reco too
        if (m_has_reco_binning)
            vars_reco_bins.push_back(m_vars_vec[i]->GetRecoBinVec(););

        // Make the name and axis label. Won't use name if input name is defined.
        name += m_vars_vec[i]->GetName();
        if (i < (m_dimension - 1))
            name += "_";
        if ((m_analysis_type == k2D || m_analysis_type == k2D_lite) && i == 1) {  // if doing 2D, set y axis label and don't add it to the linearized axis label
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

    // TODO: If doing lite types still want to include the under/overflow bins?
    int n_lin_bins = m_hyperdim->GetNLinBins();
    std::vector<double> lin_binning;
    for (int i = 0; i < n_lin_bins + 1; i++)
        lin_binning.push_back(i);

    m_lin_binning = lin_binning;

    // If you have reco binning, do that here too
    if (m_has_reco_binning) {
        m_vars_reco_binnings = vars_reco_bins;
        m_reco_hyperdim = new HyperDimLinearizer(vars_reco_bins, m_analysis_type);
        int n_lin_reco_bins = m_reco_hyperdim->GetNLinBins();

        std::vector<double> lin_reco_binning;
        for (int i = 0; i < n_lin_reco_bins + 1; i++) 
            lin_reco_binning.push_back(i);
        
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
    // k2D and k1D have under/overflow bins in phase space spread throughout,
    // and k2D_lite and k1D_lite consolodate under/overflow and tacks them on the end, so this number may seem bigger than expected.
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

    if (m_analysis_type == k2D || m_analysis_type = k2D_lite) {
        std::cout << " y axis: ";
        for (const auto b : m_vars_vec[1]->GetBinVec())
            std::cout << b << " ";
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

// Get the phase space volume for a given linearized bin. If doing a 2D type, this does not include the bin width of y axis.
// TODO: Maybe this belongs to hyperdim?
template <class UNIVERSE>
double VariableHyperDBase<UNIVERSE>::GetBinVolume(int lin_bin) const {
    double ps_bin_vol = 1.;
    // Given the linearized bin number, get corresponding bin index in phase space coordinates
    std::vector<int> ps_coords = m_hyperdim->GetValues(lin_bin);  // For 2D type, this puts 0 for y-axis.
    for (int i = 0; i < ps_coords.size(); i++) {
        if ((m_analysis_type == k2D || m_analysis_type == k2D_lite) && i == 1)  // Skip y-axis if doing 2D type
            continue;

        int var_bin = ps_coords[i];
        if (var_bin == 0 || var_bin == (m_vars_vec[i]->GetNBins() + 1)) {  // Break if under/overflow
            std::cout << "VariableHyperDBase::GetBinVolume: WARNING bin requested is under/overflow. Returning 1." << std::endl;
            return 1.;
        }
        std::vector<double> var_binning = m_vars_vec[i]->GetBinVec();  // Get the binning for that axis
        double bin_width = var_binning[var_bin] - var_binning[var_bin - 1];
        if (bin_width < 0)
            std::cout << "VariableHyperDBase::GetBinVolume: WARNING, negative binwidth, will affect bin volume" << std::endl;
        else if (bin_width == 0)
            std::cout << "VariableHyperDBase::GetBinVolume: WARNING, 0 binwidth, will affect bin volume" << std::endl;
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
        for (const auto b : GetRecoBinVec()) 
            std::cout << b << " ";
        std::cout << "\n";

        if (m_analysis_type == k2D) {
            std::cout << " y axis: ";
            for (const auto b : m_vars_vec[1]->GetRecBinVec()) 
                std::cout << b << " ";
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

// Get the phase space volume for a given linearized reco bin. If doing a 2D type, this does not include the bin width of y-axis.
// TODO: Maybe this belongs to hyperdim?
template <class UNIVERSE>
double VariableHyperDBase<UNIVERSE>::GetRecoBinVolume(int lin_bin) const {
    if (!m_has_reco_binning) {  // If there's no reco binning, just return truth instead
        return GetBinVolume(lin_bin);
    }
    // Given the linearized bin number, get corresponding bin number in phase space coordinates
    double ps_bin_vol = 1.;
    std::vector<int> ps_coords = m_reco_hyperdim->GetValues(lin_bin);   // For 2D type, this puts 0 for y-axis
    for (int i = 0; i < ps_coords.size(); i++) {
        if ((m_analysis_type == k2D || m_analysis_type == k2D_lite) && i == 1)  // Skip y-axis if doing 2D type
            continue;

        int var_bin = ps_coords[i];
        if (var_bin == 0 || var_bin == (m_vars_vec[i]->GetNRecoBins() + 1)) {  // Break if under/overflow
            std::cout << "VariableHyperDBase::GetRecoBinVolume: WARNING bin requested is under/overflow. Returning 1." << std::endl;
            return 1.;
        }
        std::vector<double> var_binning = m_vars_vec[i]->GetBinVec();  // Get the binning for that axis
        double bin_width = var_binning[var_bin] - var_binning[var_bin - 1];
        if (bin_width < 0)
            std::cout << "VariableHyperDBase::GetRecoBinVolume: WARNING, negative binwidth, will affect bin volume" << std::endl;
        else if (bin_width == 0)
            std::cout << "VariableHyperDBase::GetRecoBinVolume: WARNING, 0 binwidth, will affect bin volume" << std::endl;
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
    for (int i = 0; i < m_dimension; i++) 
        val_vec.push_back(m_vars_vec[i]->GetRecoValue(universe, idx1, idx2));
    if (!m_has_reco_binning) 
        return ((m_hyperdim->GetBin(val_vec)).first) + 0.0001;  // 0.0001 offset to so value isn't exactly on a bin edge and fillers can put it in that bin
    return ((m_reco_hyperdim->GetBin(val_vec)).first) + 0.0001;  // If there's reco binning, use that hyperdim
}

// Return bin index value in linearized bin space, ie which an bin event goes in
template <class UNIVERSE>
double VariableHyperDBase<UNIVERSE>::GetTrueValue(const UNIVERSE &universe,
                                                  const int idx1,
                                                  const int idx2) const {
    std::vector<double> val_vec;
    for (int i = 0; i < m_dimension; i++) 
        val_vec.push_back(m_vars_vec[i]->GetTrueValue(universe, idx1, idx2));
    return ((m_hyperdim->GetBin(val_vec)).first) + 0.0001;  // 0.0001 offset to so fillers can put it in that bin
}

// Return phase space value of a single component variable, use axis = 1 to get y-value for 2D type analyses.
template <class UNIVERSE>
double VariableHyperDBase<UNIVERSE>::GetRecoValue(const int axis,
                                                  const UNIVERSE &universe,
                                                  const int idx1,
                                                  const int idx2) const {
    return m_vars_vec[axis]->GetRecoValue(universe, idx1, idx2);
}

// Return phase space value of a single component variable, use axis = 1 to get y-value for 2D type analyses.
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
    for (int i = 0; i < m_dimension; i++) 
        value_vec.push_back(m_vars_vec[i]->GetRecoValue(universe, idx1, idx2));
    
    return value_vec;
}

// Return phase space values of every component variable as a vector
template <class UNIVERSE>
std::vector<double> VariableHyperDBase<UNIVERSE>::GetTrueValueVec(const UNIVERSE &universe,
                                                                  const int idx1,
                                                                  const int idx2) const {
    std::vector<double> value_vec;
    for (int i = 0; i < m_dimension; i++) 
        value_vec.push_back(m_vars_vec[i]->GetTrueValue(universe, idx1, idx2));
    return value_vec;
}

#endif  // VARIABLEHYPERDBASE_CXX
