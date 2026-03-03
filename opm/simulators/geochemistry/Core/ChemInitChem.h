/*
 * MIT License
 *
 * Copyright (C) 2025 Aksel Hiorth
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
*/
#ifndef CHEM_INIT_CHEM_H
#define CHEM_INIT_CHEM_H

#include <algorithm>
#include <array>
#include <memory>
#include <string>
#include <utility>
#include <vector>
#include <map>

#include <opm/simulators/geochemistry/Core/ChemParam.h>
#include <opm/simulators/geochemistry/Core/ChemTable.h>
#include <opm/simulators/geochemistry/Core/MineralKinetics.hpp>
#include <opm/simulators/geochemistry/Core/SurfaceChemistry.hpp>

#include <opm/simulators/geochemistry/Common/Enums.hpp>
#include <opm/simulators/geochemistry/Common/Constants.hpp>
#include <opm/simulators/geochemistry/Common/CustomExceptions.hpp>
#include <opm/simulators/geochemistry/Core/GeoChemPhases.hpp>
#include <opm/simulators/geochemistry/IO/ErrorMsg.hpp>
#include <opm/simulators/geochemistry/IO/FileFunctions.hpp>
#include <opm/simulators/geochemistry/IO/GeoLogger.hpp>
#include <opm/simulators/geochemistry/IO/InputReader.hpp>
#include <opm/simulators/geochemistry/IO/ParseChemistry.hpp>
#include <opm/simulators/geochemistry/Utility/HelperMacros.hpp>
#include <opm/simulators/geochemistry/Utility/HelpfulStringMethods.hpp>
#include <opm/simulators/geochemistry/LinearAlgebra/LinearSolverRoutines.hpp>

#include <optional>

/* Class for initializing the chemical solver. */
class InitChem
{

// Just as a safety measure, explicitly prevent copy, move, and assignment.
// (Cannot have default copy ctor anyhow b/c unique_ptr is not copyable...)
  MAKE_NONCOPYABLE(InitChem)
  MAKE_NONMOVABLE(InitChem)

  public:

    static std::unique_ptr<InitChem> create_from_input_data
        (
            InitChem* ICS_full,
            const std::map<std::string, std::vector<std::string>>& db_changes_made_by_user,
            const std::map<GeochemicalPhaseType, GeoChemPhaseData>& phases_to_be_used,
            int max_size,
            const std::string& unique_identifier,
            const std::vector<std::string>& species_names_in_order=std::vector<std::string>(),
            int interpolate = 0
        );

    static InitChem* create_from_ICS_full
        (
            const InitChem& ICS_full,
            const double* c_bas,
            const double* c_min,
            const std::vector<double>& cb_cp
        );

    void possibly_change_inconsistent_options(bool enforce_charge_balance_when_using_diffuse_layer=false);

    int read_geochemical_database(const std::map<std::string, std::vector<std::string>>& databaseReadFromInput);
    void set_relative_position_and_mineral_conc(InitChem& ICS);

    [[nodiscard]] bool includes_ion_exchange() const;
    [[nodiscard]] bool includes_surface_complexes() const;

    // Also checks mineral weight fractions in ICS for consistency.
    [[nodiscard]] double compute_rock_density(double porosity) const;

    [[nodiscard]] std::vector<std::string> get_all_basis_species() const;
    [[nodiscard]] std::vector<std::string> get_all_secondary_species() const;
    [[nodiscard]] std::vector<std::string> get_all_minerals() const;
    [[nodiscard]] std::vector<std::string> get_basis_species_of_type(int component_type) const;
    [[nodiscard]] std::vector<std::string> get_secondary_species_of_type(int component_type) const;


    [[nodiscard]] std::string name() const;
    [[nodiscard]] int get_species_index(const std::string& specie_name) const;
    [[nodiscard]] int get_mineral_index(const std::string& mineral_name) const;
    [[nodiscard]] std::string get_basis_name(int pos) const;
    [[nodiscard]] std::string get_mineral_name(int pos) const;

    void print_species_to_screen() const;
    void write(const std::string& file_name) const;

  public:

    double Temp_{DefaultValues::Temperature};
    double Pres_{DefaultValues::Pressure};

    int PRINT_DEBUG_CHEM_{DebugInfoLevel::NONE};
    int INTERPOLATE_{0};  // Splay tree resolution

    bool REMOVE_ORGANIC_SPECIES{ true }; 
    bool charge_balance_{false};  // Calculate pH from charge balance?
    bool keep_ph_fixed_{ false }; // allow use to specify FIX to keep pH constant in simulation
    double WHTO_{1.0}; // Concentration of water (H2O)

    // NB:  It is the version in BasVec that is used. (They might be different!)
    SurfaceChemistryOptions surface_options_from_input_{};
    double SA_{1.0}; // Specific surface area (m^2/L PV) for surface complex. model.
    double d_DL_{0.0}; // Size of the electric double layer (nm).
    double gas_phase_frac_{ 1.0 }; //Volume fraction of a gas phase in contact with a solution 

    int max_size_{0};  // Ensure that enough space is allocated (only once!)
    int size_basis_{0};
    int size_aq_{0};
    int size_io_{0};
    int size_surf_{0};
    int size_mass_{0};
    int size_min_{0};
    int size_sup_min_{0};
    int size_rock_{0};
    int size_gas_phase_{ 0 };
    int size_gas_{0};  // Gas buffers in oil or gas phase

    // Index mappings: { local species } --> { corr. species in chem. databases (ChemTables) }
    std::vector<int> pos_;
    std::vector<int> pos_mass_;
    std::vector<int> pos_buffer_;
    std::vector<int> pos_rock_;
    std::vector<int> pos_min_;
    std::vector<int> pos_sup_min_;
    std::vector<int> pos_min_bas_;
    std::vector<int> pos_sup_bas_;
    std::vector<int> pos_gas_;
    std::vector<int> pos_gas_phase_;

    

    std::vector<int> phase_;

    std::vector<int> kcmap_;  // Index mapping: { local basis species }   -->   { ordering used by a solver }
    std::vector<int> krmap_;  // Index mapping  { local mineral species } -->   { ordering used by a solver }

    std::vector<std::string> basis_species_name_;
    std::vector<std::string> io_name_;
    std::vector<std::string> surf_name_;
    std::vector<std::string> mineral_name_;
    std::vector<std::string> sup_min_name_;
    std::vector<std::string> buffer_name_;
    std::vector<std::string> gas_phase_name_;

    // Used to initialize the actual geochemical database used
    std::vector<std::string> basis_spec_db_;
    std::vector<std::string> sec_spec_db_;
    std::vector<std::string> surf_spec_db_;
    std::vector<std::string> exch_spec_db_;
    std::vector<std::string> min_phas_db_;

    std::vector<std::string> species_to_remove_;
    std::vector<std::string> bspecies_always_included_{ "H", "H2O" };
	std::vector<std::string> species_added_after_input_file_;														 

    // Mol weights of water, oil, and gas in kg / mol
    std::array<double, 3> mol_weight_phase_{ 18.01528e-3,
                                             0.470*std::pow(0.8, 5.57),
                                             28.9647e-3 };

    std::vector<double> c_vchem_;
    std::vector<double> c_ads_;
    std::vector<double> c_io_;
    std::vector<double> c_surf_;

    std::vector<double> c_mineral_;
    std::vector<double> c_sup_min_;
    std::vector<double> c_buffer_;
    std::vector<double> gas_phase_pressure_;
    std::vector<double> log_a_mineral_;
    std::vector<double> log_af_sup_;
    std::vector<double> log_af_;  // (buffer)
    
    std::vector<double> Sg_;



    std::vector<std::optional<MineralRateLaw>> rate_{}; // Rate constants for minerals

    ChemParam CP_{};

    std::unique_ptr<ChemTable> SM_basis_{nullptr};  // Basis species database
    std::unique_ptr<ChemTable> SM_all_{nullptr};  // Aqueous complexes database
    std::unique_ptr<ChemTable> SM_mineral_{nullptr};  // Mineral database

    // Full databases
    std::unique_ptr<ChemTable> Basis_db_{nullptr};
    std::unique_ptr<ChemTable> Aq_db_{nullptr};
    std::unique_ptr<ChemTable> Mineral_db_{nullptr};

  private:

    std::string name_;
    std::map<std::string, std::string> basis_transformation_;

    // Helper functions for initialization
    void initialize_vectors(int max_size);
    void reorder_species(const std::vector<std::string>& speciesNamesInOrder);
    void set_phase_composition(const std::map<GeochemicalPhaseType, GeoChemPhaseData>& phases);
    void init_mineral_vector_in_correct_order();

    void check_gas_phase_for_equilibrium_phases(const std::vector<std::string>& gas_phase, const std::vector<std::string>& buffer);

    void init_chemistry_hkf();
    void init_chemistry_hkf(const ChemTable& SM_basis, const ChemTable& SM_all, const ChemTable& SM_mineral);
    void finalize_init_chemistry_hkf();

    void init_phase_mol_weights();
    void chem_set_pos_relative_to_full_db();
    void db_set_activity();
    void fix_hkf_units() const;

    // Modify database based on input?
    void add_basis_species(const std::vector<std::string>& buffer, int type);
    void add_secondary_species(const std::vector<std::string>& buffer);
    void add_surface_species(const std::vector<std::string>& buffer, std::vector<std::string>& bnames);
    void add_mineral_phases(const std::vector<std::string>& buffer);
    int add_mineral_basis_species(const std::vector<std::string>& mineral_names, std::vector<std::string>& bnames);
    std::vector<std::string> remove_organic_species();

    // functions for basis transformation
    void clean_transformation();
    void basis_transformation();

    // Helper functions to parse previously read input (stored in std::map's)
    std::vector<std::string> read_specie(const std::string& basis_specie_name, const std::string& line, double& conc) const;
    int read_solution(const GeoChemBlockContent& block_content);
    int read_rate(const GeoChemBlockContent& block_content);
    int read_ieq(const GeoChemBlockContent& block_content);
    int read_gas_phase(const GeoChemBlockContent& block_content);
    int read_io(const GeoChemBlockContent& block_content);
    int read_surface(const GeoChemBlockContent& block_content);

    // Other helper functions to set up and verify the database etc.
    [[nodiscard]] bool basis_specie_defined(std::size_t pos, const std::vector<std::string>& names) const;
    static void check_db(ChemTable* Test, ChemTable* Basis) ;
    void set_diffuse_layer_options_from_read_input();

    // See, e.g.: https://seanmiddleditch.github.io/enabling-make-unique-with-private-constructors/
    struct _constructor_tag { explicit _constructor_tag() = default; };

    InitChem(
        int debug_mode,
        bool charge_balance,
        const ChemTable& SM_basis,
        const ChemTable& SM_all,
        const ChemTable& SM_mineral,
        const std::vector<std::string>& basis_name,
        const std::vector<double>& basis_conc,
        const std::vector<int>& basis_offset_index,
        const std::vector<std::string>& mineral_name,
        const std::vector<double>& mineral_conc,
        const std::vector<int>& mineral_offset_index,
        const std::vector<double>& mineral_log_af,
        const std::vector<int>& pos_gas,
        const std::vector<int>& phase,
        const std::vector<MineralRateLaw>& rate_laws,
        const std::vector<std::string>& io_name,
        const std::vector<std::string>& surf_name,
        const SurfaceChemistryOptions& surf_options,
        double specific_surface_area,
        double size_of_the_diffuse_layer
    );

  public:
    InitChem(_constructor_tag,
             InitChem* ICS_full,
             const std::map<std::string, std::vector<std::string>>& db_changes_made_by_user,
             const std::map<GeochemicalPhaseType, GeoChemPhaseData>& phases_to_be_used,
             int max_size,
             const std::string& unique_identifier,
             const std::vector<std::string>& species_names_in_order=std::vector<std::string>(),
             int interpolate = 0);
};

/** Helper struct to store numerical data corresponding to a given InitChem object. */
struct InitChemData
{
    InitChemData(std::string solutionName, std::size_t size_aq, std::size_t size_minerals)
        : solutionName_(std::move(solutionName))
        // Potentially changing:
        , equilibrated_(false)
        , pH_(7.0)
        , sigma_(0.0)
        , psi_(0.0)
        , C_aq_(size_aq, 0.0)
        , C_ads_(size_aq, 0.0)
        , C_DL_(size_aq, 0.0)
        , C_minerals_(size_minerals, 0.0)
        , Log_a_minerals_(size_minerals, 0.0)
    {
    }

    void reset()
    {
        equilibrated_ = false;

        pH_ = 7.0;
        sigma_ = 0.0;
        sigma_ = 0.0;

        fill_zero(C_aq_);
        fill_zero(C_ads_);
        fill_zero(C_DL_);
        fill_zero(C_minerals_);
        fill_zero(Log_a_minerals_);
    }

    std::string solutionName_;

    bool equilibrated_{false};

    double pH_{7.0};
    double sigma_{0.0};  // C/m^2
    double psi_{0.0};

    std::vector <double> C_aq_;
    std::vector <double> C_ads_;
    std::vector <double> C_DL_;
    std::vector <double> C_minerals_;
    std::vector <double> Log_a_minerals_;
};

#endif  // CHEM_INIT_CHEM_H
