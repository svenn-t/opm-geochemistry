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
/**
* @file
*/
#ifndef GEO_CHEM_IF_IS_DEF_H
#define GEO_CHEM_IF_IS_DEF_H

#include <nlohmann/json.hpp>

#include <opm/simulators/geochemistry/Common/ChemGlobal.h>
#include <opm/simulators/geochemistry/Common/Constants.hpp>
#include <opm/simulators/geochemistry/Common/SerializeForTesting.hpp>
#include <opm/simulators/geochemistry/Core/SimulationSchedule.hpp>
#include <opm/simulators/geochemistry/Core/ChemBasVec.h>
#include <opm/simulators/geochemistry/Core/ChemGCSolver.h>
#include <opm/simulators/geochemistry/Core/ChemInitChem.h>
#include <opm/simulators/geochemistry/Core/ChemState.hpp>
#include <opm/simulators/geochemistry/Core/GeoChemPhases.hpp>
#include <opm/simulators/geochemistry/Core/GeoChemSolutionsManager.hpp>
#include <opm/simulators/geochemistry/IO/InputReader.hpp>
#include <opm/simulators/geochemistry/IO/NetCDF.hpp>
#include <opm/simulators/geochemistry/IO/ParseUtil.hpp>
#include <opm/simulators/geochemistry/SplayTree/SplayTree.h>
#include <opm/simulators/geochemistry/Utility/HelperMacros.hpp>
#include <opm/simulators/geochemistry/Utility/HelpfulStringMethods.hpp>
#include <opm/simulators/geochemistry/Utility/MiscUtilityFunctions.hpp>

/** Base class for all interface classes of the geochemical solver:
 *
 *      - EquilibriumSolver
 *      - OneDimensionalTransportSolver
 *
 * Instances of CGeoChemIF should no longer be created directly, instead we
 * always use one of the subclasses. The top-level class is only used to
 * share (mostly) common state and functions.
 */
class CGeoChemIF
{
  MAKE_NONCOPYABLE(CGeoChemIF)
  MAKE_NONMOVABLE(CGeoChemIF)

  public:
    //CGeoChemIF();
    explicit CGeoChemIF(std::unique_ptr<InitChem> initChem = nullptr);
    explicit CGeoChemIF(std::unique_ptr<InitChem> initChem, int INTERPOLATE);

    void printCurrentVchemToFile(const std::string& file_name);
    void printSplayTreeInfo() const;
    int getSolverCalls() const;

    [[nodiscard]] std::size_t numberOfSolutions() const;
    [[nodiscard]] std::size_t numberOfBasisSpecies() const;
    [[nodiscard]] std::size_t numberOfAqueousBasisSpecies() const;
    [[nodiscard]] std::size_t numberOfMinerals() const;
    [[nodiscard]] std::size_t numberOfIonExchange() const;

    [[nodiscard]] std::vector<std::string> GetSpeciesNames(int component_type,
                                                           const std::string& postfix = "") const;
    [[nodiscard]] std::vector<std::string> GetAdsorbedSpeciesNames() const;
    [[nodiscard]] std::vector<std::string> GetMineralNames() const;
    [[nodiscard]] std::vector<std::string> GetBasisSpeciesNames() const;
    [[nodiscard]] std::vector<std::string> GetIONames() const;
    [[nodiscard]] std::vector<std::string> GetSurfNames() const;

    [[nodiscard]] double GetSurfaceArea() const;
    [[nodiscard]] double GetDiffusionLayer() const;

    void SetSurfaceConc(double swat,
                        std::vector<double>& C_tot,
                        double& frac_DL,
                        InitChem* ICS_in = nullptr) const;

    InitChemData* getDataForSolution(const std::string& solutionName);
    [[nodiscard]] std::vector<std::vector<double>> GetInitialSolutions() const;

    void initialize_geochemical_state(std::size_t no_grid_blocks);
    void update_geochemical_state(std::size_t block_index);

    [[nodiscard]] bool basis_species_added_by_geochem();

    void SolveChem_I(double* ctot,
                     double* ctot_ads,
                     double* ctot_min,
                     double* log_a_min,
                     double Temp,
                     double Pres,
                     double porosity,
                     double dt,
                     double SA,
                     double frac_DL,
                     const std::array<double, 3>& mass_phase,
                     double& pH_in,
                     double& sigma_in,
                     double& psi_in);

  protected:
    void AddPhases(const PairOfPhaseNamesAndTypes& names_and_types,
                   const std::string& unique_identifier);

    /** Calls addPhases() assuming a single solution (i.e., no minerals,
     *  surface species, etc.). */
    void AddAqueousSolution(const std::string& solution_name);

  public:

    std::unique_ptr<InputReader> ptrInputReader_{nullptr};

    GeoChemState grid_state_for_netcdf_{};
    std::unique_ptr<NetCDFWriter> ptrNetCDF_{nullptr};
    GeoChemPhases geoChemPhases_{};

  protected:

    // These are the maximum numbers that can occur during the simulation
    int no_basis_species_{0};
    int no_aq_species_{0};
    int no_minerals_{0};
    int no_io_{0};

    std::vector<double> Cmin_;
    std::vector<double> Cchem_;
    std::vector<double> Ccads_;
    std::vector<double> Caq_;
    std::vector<double> Cdl_;
    std::vector<double> Charge_;

    std::array<double, 3> mass_phase_{ 1.0, 1.0, 1.0 };

    double Temp_{DefaultValues::Temperature};

    double sigma_{0.0};  // Surface charge
    double psi_{0.0};    // Surface potential
    double WHTO_{0.0};   // H2O

  protected:
    GeoLogger logger_{};

    std::unique_ptr<GCSolver> GCS_{nullptr};

    InputReader geochemicalDatabaseReader_;
    std::map<std::string, std::vector<std::string>> db_changes_made_by_user_;

    using InitChemPtr = std::unique_ptr<InitChem>;
    InitChemPtr ICS_full_{nullptr};
    static const inline std::string ICS_full_name_ = "ICS_FULL_";
    std::vector<InitChemPtr> ICS_sol_;
    std::vector<InitChemData> ICS_sol_data_;

    GeochemicalSolutionsManager solutions_manager_{};
    std::vector<double> val_log_;

    std::vector<double> c_old_;
    std::vector<double> dc_min_;
    std::vector<double> c_flux_;
    std::vector<double> c_min_old_;

  protected:
    void set_input_reader_keywords();

    void allocate_memory_for_solver_and_splay_tree_etc();
    void modify_selected_solver_options(GCSolver& GCS_in) const;
    void modify_selected_solver_options_json(GCSolver& GCS_in, nlohmann::json& json_parsed) const;

    void set_geochemical_database_modifications_from_user_input();
    void set_geochemical_database_modifications_from_json(const std::string& json_input);

    void create_ICS_full(const std::vector<std::string>& species_names_in_order={});

    void set_data_for_solution(const std::string& solutionName, double porosity);
    void set_data_for_solution(std::size_t solIdx, double porosity);

    void equilibrate_solution
        (
            // Requires a pointer to construct a BasVec object from ICS:
            InitChem* ICS,
            bool equilibrate,
            const std::string& case_name,
            int GEOCHEMDEBUG,
            std::vector<double>& Ct,
            std::vector<double>& Cads,
            double& pH,
            std::vector<double>& scharge,
            BasVecInfo* basVecInfo = nullptr,
            double frac_DL = -1.
        );
};

#endif  // GEO_CHEM_IF_IS_DEF_H
