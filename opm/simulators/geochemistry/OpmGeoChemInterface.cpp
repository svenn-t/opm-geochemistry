/*
  Copyright 2025 Equinor ASA.

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/
#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/simulators/geochemistry/OpmGeoChemInterface.hpp>
#include <opm/simulators/geochemistry/Utility/HelpfulStringMethods.hpp>
#include <opm/simulators/geochemistry/Core/GeoChemPhases.hpp>

#include <fmt/format.h>

#include <fstream>
#include <stdexcept>

namespace {

// Print warning helper
auto printWarning (const std::string& file_name, const std::string& block, const char* keyword, bool ignore)
{
    std::string msg;
    if (ignore) {
        msg = fmt::format("Block \"{}\" in JSON file {} will be ignored! Use keyword {} instead!",
                            block, file_name, keyword);
    }
    else {
        msg = fmt::format("Block \"{}\" in JSON file {} will be overwritten by keyword {}!",
                            block, file_name, keyword);
    }
    Opm::OpmLog::warning(msg);
};

}

// ////
// PUBLIC METHODS
// ///
void OpmGeoChemInterface::initialize_from_opm_deck(const std::string& file_name,
                                                   const std::vector<std::string> species,
                                                   const std::vector<std::string> minerals,
                                                   const std::vector<std::string> ion_ex,
                                                   bool charge_balance,
                                                   std::pair<double, double> tol,
                                                   int splay_tree_resolution)
{
    // Append species from OPM deck to JSON
    nlohmann::json appended_json = opmDeckSpeciesToJSON_(file_name, charge_balance, species, minerals, ion_ex);

    // Parse JSON
    auto json_parsed = geoChemPhases_.resetFromJson(appended_json.dump());
    std::map<GeochemicalPhaseType, GeoChemPhaseData> phases_to_be_used;
    const int max_size = geoChemPhases_.selectPhases(geoChemPhases_.getPhaseNamesAndTypes(), phases_to_be_used);
    set_geochemical_database_modifications_from_json(appended_json.dump());

    // Init.
    ICS_full_ = InitChem::create_from_input_data(nullptr,
                                                 db_changes_made_by_user_,
                                                 phases_to_be_used,
                                                 max_size,
                                                 "TESTING",
                                                 species,
                                                 splay_tree_resolution);
    allocate_memory_for_solver_and_splay_tree_etc();

    // Solver tolerances: tol = {mbal, ph}
    if (json_parsed.contains("chemtol") || json_parsed.contains("CHEMTOL")) {
        printWarning(file_name, "CHEMTOL", "GEOCHEM", /*ignore=*/false);
    }
    set_solver_tolerances(*GCS_, tol.first, tol.second);

    // Check species defined in OPM deck vs. geochemical solver
    checkUserOrder_(species);
}

void OpmGeoChemInterface::initialize_json(const std::string& file_name,
                                          const std::vector<std::string> user_order)
{
    // Read input from JSON file
    auto json_parsed = geoChemPhases_.resetFromJson(file_name);
    std::map<GeochemicalPhaseType, GeoChemPhaseData> phases_to_be_used;
    const std::map<std::string, std::vector<std::string>> db_changes_made_by_user;
    const int max_size = geoChemPhases_.selectPhases(geoChemPhases_.getPhaseNamesAndTypes(), phases_to_be_used);

    // Init
    ICS_full_ = InitChem::create_from_input_data(nullptr,
                                                 db_changes_made_by_user,
                                                 phases_to_be_used,
                                                 max_size,
                                                 "TESTING",
                                                 user_order,
                                                 0);
    allocate_memory_for_solver_and_splay_tree_etc();

    // Solver option
    // TODO: Add to JSON parser instead of function below?
    // Add option to change ICS_full_->INTERPOLATE ?
    modify_selected_solver_options_json(*GCS_, json_parsed);  // HACK(?)

    // Check species defined in OPM deck vs. geochemical solver
    checkUserOrder_(user_order, file_name);
}

void OpmGeoChemInterface::initialize(const std::string& file_name,
                                     double temperature,
                                     double porosity,
                                     const std::vector<std::string> user_order)
{
    if (ICS_full_) throw MultipleInitializationException("Cannot initialize ICS_full_, it already exists...");

    std::ifstream inputStream(file_name, std::ios::binary);
    ptrInputReader_->read(inputStream);
    set_geochemical_database_modifications_from_user_input();

    geoChemPhases_.resetFromInputStream(inputStream);

    // Get all solutions entered, and mix all solutions into one solution (same for minerals and surfaces)
    // NB: Must set interpolate flag before allocating memory for solver, splay tree, etc.
    create_ICS_full(user_order);

    ICS_full_->INTERPOLATE_ = std::stoi(ptrInputReader_->get_simple_keyword_value("INTERPOLATE"));
    ICS_full_->PRINT_DEBUG_CHEM_ = std::stoi(ptrInputReader_->get_simple_keyword_value("DEBUG"));
    ICS_full_->possibly_change_inconsistent_options();

    allocate_memory_for_solver_and_splay_tree_etc();
    modify_selected_solver_options(*GCS_);

    Temp_ = temperature; // ICS_full_->Temp_;

    // We start by adding all phases initially in the reservoir, then we add the injected solutions.
    auto [reservoir_phases, names_of_injected_solutions] = geoChemPhases_.GetInjectReservoirSimple();
    AddPhases(reservoir_phases, "RESERVOIR");
    for (const auto& solution_name: names_of_injected_solutions)
    {
        AddAqueousSolution(solution_name);
    }

    for (int solutionIndex = 0; solutionIndex < static_cast<int>(ICS_sol_.size()); ++solutionIndex)
    {
        // Equilibration is done locally in each block (later).
        set_data_for_solution(solutionIndex, porosity);
    }
}

void OpmGeoChemInterface::calculate_initial_mineral_concentration(std::vector<double>& Cmin,
                                                                  double porosity,
                                                                  std::unordered_map<std::string, double> weight_mineral)
{
    assert(weight_mineral.size() == static_cast<std::size_t>(ICS_full_->size_min_));

    for (const auto& [name, val] : weight_mineral) {
        int minIdx = ICS_full_->get_mineral_index(name);
        if (minIdx < 0) {
            const std::string msg = fmt::format("Mineral = {} not found internally in geochemistry solver!", name);
            Opm::OpmLog::error(msg);
            throw std::runtime_error(msg);
        }
        ICS_full_->c_mineral_[minIdx] = val;
    }

    const double rock_density = ICS_full_->compute_rock_density(porosity);  // kg/m^3

    Cmin.resize(ICS_full_->size_min_);
    for (int bufferIdx = 0; bufferIdx < ICS_full_->size_min_; ++bufferIdx)
    {
        const auto name_of_buffer = ICS_full_->get_mineral_name(bufferIdx);
        const double input_concentration = ICS_full_->c_mineral_[bufferIdx];
        const double mol_weight = ICS_full_->SM_mineral_->mol_weight_[ICS_full_->pos_min_[bufferIdx]];

        if(is_gas_buffer(name_of_buffer))
        {
            Cmin[bufferIdx] = input_concentration / mol_weight;
        }
        else
        {
            const double mineral_wt_frac = input_concentration;
            const double fac = 1.0e-3*(rock_density/mol_weight)*(1.0-porosity)/porosity;
            Cmin[bufferIdx] = fac*mineral_wt_frac;
        }
    }
}

void OpmGeoChemInterface::set_surface_concentrations(double swat,
                                                     std::vector<double>& C_tot,
                                                     double& frac_DL,
                                                     std::unordered_map<std::string, double> C_io)
{
    // Set ion exchange concentration
    assert(C_io.size() == static_cast<std::size_t>(ICS_full_->size_io_));

    for (int i = 0; i < ICS_full_->size_io_; ++i) {
        const auto& name = ICS_full_->io_name_[i];
        ICS_full_->c_io_[i] = C_io.at(name);
    }

    return SetSurfaceConc(swat, C_tot, frac_DL);
}

void OpmGeoChemInterface::set_solver_tolerances(GCSolver& GCS_in,
                                                double mbal_tol,
                                                double ph_tol)
{
    GCS_in.options_.MBAL_CONV_CRITERION_ = mbal_tol;
    GCS_in.options_.PH_CONV_CRITERION_ = ph_tol;
}

std::vector<double>& OpmGeoChemInterface::get_log_a_mineral() const
{ return ICS_full_->log_a_mineral_; }

// ////
// PRIVATE METHODS
// ///
nlohmann::json OpmGeoChemInterface::appendUserSpeciesToJSON_(const std::string& file_name,
                                                             const std::vector<std::string> species_names)
{
    nlohmann::json json_parsed;
    try {
        std::ifstream file(file_name);
        if (file.good()) {
            std::ifstream in(file_name);
            if (!in) {
                throw std::runtime_error("Could not open JSON file: " + file_name);
            }
            in >> json_parsed;
        } else {
            // Fallback: parse as raw JSON string
            json_parsed = nlohmann::json::parse(file_name);
        }
    }
    catch (const nlohmann::json::parse_error& e) {
        std::cerr << "JSON parse error in appendUserSpeciesToJSON_: " << e.what() << std::endl;
    }

    // Make sure there is an aqueous solution section in the JSON file
    const std::string solution_section = PhaseKeyword(GeochemicalPhaseType::AQUEOUS_SOLUTION);
    if (!json_parsed.contains(solution_section)) {
        json_parsed[solution_section] = nlohmann::json::object();
    }

    if (!json_parsed[solution_section].contains(solution_section + " 0")) {
        json_parsed[solution_section][solution_section + " 0"] = nlohmann::json::object();
    }

    // Append species names with an (arbitrary) initial concentration
    for (const auto& species : species_names) {
        // Skip species already in JSON file
        if (json_parsed[solution_section][solution_section + " 0"].contains(species)) {
            continue;
        }

        // Generate key-value block for species
        // TODO: H with charge
        json_parsed[solution_section][solution_section + " 0"][species] = "1e-12";
    }

    return json_parsed;
}

nlohmann::json OpmGeoChemInterface::opmDeckSpeciesToJSON_(const std::string& file_name,
                                                          bool charge_balance,
                                                          const std::vector<std::string> species,
                                                          const std::vector<std::string> minerals,
                                                          const std::vector<std::string> ion_ex)
{
    // SPECIES keyword required!
    if (species.empty()) {
        const std::string msg = "SPECIES keyword required in geochemistry solver!";
        Opm::OpmLog::error(msg);
        throw std::runtime_error(msg);
    }

    nlohmann::json json_opm;
    if (!file_name.empty()) {
        // Read JSON file
        std::ifstream in(file_name);
        if (!in) {
            const std::string msg = fmt::format("Could not open JSON file: {}", file_name);
            Opm::OpmLog::error(msg);
            throw std::runtime_error(msg);
        }
        in >> json_opm;
    }

    // Append SOLUTION 0 species (i.e. aqueous species)
    const std::string solution_block = PhaseKeyword(GeochemicalPhaseType::AQUEOUS_SOLUTION);
    if (json_opm.contains(solution_block) || json_opm.contains(to_lower_case(solution_block))) {
        printWarning(file_name, solution_block, "SPECIES", /*ignore=*/false);
    }
    json_opm[solution_block] = nlohmann::json::object();
    json_opm[solution_block][solution_block + "0"] = nlohmann::json::object();
    for (const auto& elem : species) {
        if (elem == "H" && charge_balance) {
            json_opm[solution_block][solution_block + " 0"][elem] = "1.0 charge";
        }
        else {
            json_opm[solution_block][solution_block + " 0"][elem] = "1e-12";  // arbitrary value!
        }
    }

    // Append minerals
    const std::string mineral_block = PhaseKeyword(GeochemicalPhaseType::EQUILIBRIUM_MINERAL);
    if (!minerals.empty()) {
        if (json_opm.contains(mineral_block) || json_opm.contains(to_lower_case(mineral_block))) {
            printWarning(file_name, mineral_block, "MINERAL", /*ignore=*/false);
        }
        json_opm[mineral_block] = nlohmann::json::object();
        json_opm[mineral_block][mineral_block + " 0"] = nlohmann::json::object();
        for (const auto& elem : minerals) {
            json_opm[mineral_block][mineral_block + " 0"][elem] = "1";  // arbitrary value!
        }
    }
    else {
        if (json_opm.contains(mineral_block) || json_opm.contains(to_lower_case(mineral_block))) {
            printWarning(file_name, mineral_block, "MINERAL", /*ignore=*/true);
            json_opm.erase(mineral_block);
        }
    }

    // Append ion exchange
    const std::string iexchange_block = PhaseKeyword(GeochemicalPhaseType::EXCHANGE_SITES);
    if (!ion_ex.empty()) {
        if (json_opm.contains(iexchange_block) || json_opm.contains(to_lower_case(iexchange_block))) {
            printWarning(file_name, iexchange_block, "IONEX", /*ignore=*/false);
        }
        json_opm[iexchange_block] = nlohmann::json::object();
        json_opm[iexchange_block][iexchange_block + " 0"] = nlohmann::json::object();
        for (const auto& elem : ion_ex) {
            json_opm[iexchange_block][iexchange_block + " 0"][elem] = "";  // sets default value
        }
    }
    else {
        // Delete ion exchange from JSON
        if (json_opm.contains(iexchange_block)|| json_opm.contains(to_lower_case(iexchange_block))) {
            printWarning(file_name, iexchange_block, "IONEX", /*ignore=*/true);
            json_opm.erase(iexchange_block);
        }
    }

    return json_opm;
}

void OpmGeoChemInterface::checkUserOrder_(const std::vector<std::string> user_order,
                                          std::optional<std::string> file_name)
{
    int user_order_size = static_cast<int>(user_order.size());
    if (ICS_full_->size_aq_ != user_order_size){
        std::string msg;
        std::vector<std::string> species_not_present;
        bool add_or_remove;
        if (ICS_full_->size_aq_ > user_order_size) {
            for (int i = 0; i < ICS_full_->size_aq_; ++i) {
                const auto basis_name = ICS_full_->basis_species_name_[i];
                int pos_in_db = ICS_full_->SM_basis_->get_row_index(basis_name);
                const auto internal_name = to_upper_case(ICS_full_->SM_basis_->nick_name_[pos_in_db]);
                const auto it = std::find(user_order.begin(), user_order.end(), internal_name);
                if (it == user_order.end()) {
                    species_not_present.push_back(internal_name);
                }
            }
            add_or_remove = true;
        }
        else {
            for (int i = 0; i < user_order_size; ++i) {
                const auto user_species = user_order[i];
                if (ICS_full_->SM_basis_->get_row_index(user_species) < 0) {
                    species_not_present.push_back(user_species);
                }
            }
            add_or_remove = false;
        }
        if (file_name.has_value()) {
            msg = fmt::format("The following species must be {} SPECIES or {} \"{}\": ",
                              add_or_remove ? "added to" : "removed from",
                              add_or_remove ? "removed from" : "added to",
                              *file_name);
        }
        else {
            msg = fmt::format("The following species must be {} SPECIES : ",
                                add_or_remove ? "added to" : "removed from");
        }
        for (const auto& elem : species_not_present) {
            msg += fmt::format("{} ", elem);
        }
        Opm::OpmLog::error(msg);
        throw std::runtime_error(msg);
    }
}
