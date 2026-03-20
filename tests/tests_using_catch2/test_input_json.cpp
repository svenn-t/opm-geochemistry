#if __has_include(<catch2/catch.hpp>)
#include <catch2/catch.hpp>
#else
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#endif

#include <sstream>
#include <string>
#include <iostream>
#include <map>

#include <opm/simulators/geochemistry/Core/GeoChemPhases.hpp>
#include <opm/simulators/geochemistry/Core/ChemInitChem.h>
#include <opm/simulators/geochemistry/GeoChemIF.h>
#include <opm/simulators/geochemistry/Common/Enums.hpp>

class TestGeoChemIF : public CGeoChemIF
{
    public:
        TestGeoChemIF(const std::string json_string) : CGeoChemIF()
        {
            set_geochemical_database_modifications_from_json(json_string);
        }

    const std::map<std::string, std::vector<std::string>>& getDbModifications()
    { return db_changes_made_by_user_; }
};

std::string json_input =
    R"({"SOLUTION": {
            "SOLUTION 0": {
                "pH": "7 charge",
                "Cl": "2e-2",
                "Ca": "1e-2"
            }
        },
        "EQUILIBRIUM_PHASES": {
            "EQUILIBRIUM_PHASES 0": {
                "CALCITE": "1 0",
                "CO2,g": "1 -3.5"
            }
        },
        "WRONG_PHASE": {
            "EQUILIBRIUM_PHASES XX": {
                "CALCITE_XX": "1 0",
                "CO2_XX": "1 -3.5"
            }
        }
    })";

double pH_to_match_json = 7.71221;

std::string json_input2 =
    R"({"SOLUTION": {
            "SOLUTION 0": {
                "pH": "7 charge",
                "Ca": "0.6e-3",
                "NO3": "1.2e-2",
                "K": "0.2e-3",
                "Na": "1e-3",
                "Cl": "1.2e-2"
            }
        },
        "IEXCHANGE": {
            "IEXCHANGE 0": {
                "X": " 1.1e-3"
            }
        }
    })";

std::string json_db =
    R"({"SOLUTION": {
            "SOLUTION 0": {
                "pH": "7 charge",
                "Ca": "1e-2",
                "Fee+2": "1e-3",
                "Hdg": "1e-4",
                "Hdg2": "1e-5",
                "NO3": "1e-6"
            }
        },
        "IEXCHANGE": {
            "IEXCHANGE 0": {
            "Z+": "1.1e-3"
            }
        },
        "BASIS_SPECIES": [
            "Fee+2 /HKF 6 55.847 -21870 -22050 -25.3 -0.7867 -9.6969 9.5479 -2.378 14.786 -4.6437 1.4382 /",
            "Hdg / HKF 4 0 4236 -1000 13.8 5.1427 4.7758 3.8729 -2.9764 27.6251 5.093 -0.209 /",
            "Hdg2 4 0 / 0 /"
        ],
        "SECONDARY_SPECIES": [
            "ZNO3   = Z+ + NO3- /ANA -0.7/"
        ],
        "EXCHANGE_SPECIES" : [
            "Z+"
        ],
        "MINERAL_PHASES": [
            "Hdg2(G) = Hdg2 /ANA 0 -9.3114    4.6473e-3   -49.335    1.4341    1.2815e5 /"
        ],
        "REMOVE_SPECIES": [
            "HNO3",
            "CaOH+"
        ]
    })";

TEST_CASE("Test json input format")
{
    GeoChemPhases phase_data;


     phase_data.resetFromJson(json_input);

     REQUIRE(phase_data.hasPhaseWithName(GeochemicalPhaseType::AQUEOUS_SOLUTION, "SOLUTION 0"));
     REQUIRE(phase_data.hasPhaseWithName(GeochemicalPhaseType::EQUILIBRIUM_MINERAL, "EQUILIBRIUM_PHASES 0"));


     REQUIRE(phase_data.hasElementInPhase(GeochemicalPhaseType::AQUEOUS_SOLUTION, "Ca"));
     REQUIRE(phase_data.hasElementInAnyPhase("Ca"));

     REQUIRE(phase_data.hasElementInPhase(GeochemicalPhaseType::AQUEOUS_SOLUTION, "Cl"));
     REQUIRE(phase_data.hasElementInAnyPhase("Cl"));

     REQUIRE(phase_data.hasElementInPhase(GeochemicalPhaseType::EQUILIBRIUM_MINERAL, "CALCITE"));
     REQUIRE(phase_data.hasElementInAnyPhase("CO2,g"));
     // NOT PRESENT
     REQUIRE(!phase_data.hasElementInAnyPhase("CALCITE_XX"));
     REQUIRE(!phase_data.hasElementInAnyPhase("CO2_XX"));
}

TEST_CASE("GeoChemPhases JSON reset from string and file should yield same result") {

    // Create object and reset from string
    GeoChemPhases phases_from_string;
    phases_from_string.resetFromJson(json_input);

    // Write JSON to temp file
    const std::string temp_filename = "temp_test_phases.json";
    {
        std::ofstream out(temp_filename);
        REQUIRE(out.is_open());
        out << json_input;
    }

    // Create another object and reset from file
    GeoChemPhases phases_from_file;
    phases_from_file.resetFromJson(temp_filename);

    // Cleanup the temp file
    std::remove(temp_filename.c_str());

    // Compare the two objects
    // You need to define a meaningful comparison function or operator== for GeoChemPhases
    REQUIRE(phases_from_string.hasPhaseWithName(GeochemicalPhaseType::AQUEOUS_SOLUTION, "SOLUTION 0"));
    REQUIRE(phases_from_file.hasPhaseWithName(GeochemicalPhaseType::AQUEOUS_SOLUTION, "SOLUTION 0"));

    REQUIRE(phases_from_string.hasPhaseWithName(GeochemicalPhaseType::EQUILIBRIUM_MINERAL, "EQUILIBRIUM_PHASES 0"));
    REQUIRE(phases_from_file.hasPhaseWithName(GeochemicalPhaseType::EQUILIBRIUM_MINERAL, "EQUILIBRIUM_PHASES 0"));

    REQUIRE(phases_from_string.hasElementInPhase(GeochemicalPhaseType::AQUEOUS_SOLUTION, "Ca"));
    REQUIRE(phases_from_file.hasElementInPhase(GeochemicalPhaseType::AQUEOUS_SOLUTION, "Ca"));

    REQUIRE(phases_from_string.hasElementInAnyPhase("Ca"));
    REQUIRE(phases_from_file.hasElementInAnyPhase("Ca"));

    REQUIRE(phases_from_string.hasElementInAnyPhase("CALCITE"));
    REQUIRE(phases_from_file.hasElementInAnyPhase("CALCITE"));
}

TEST_CASE("Test initialization of ICS_full and order determined by user")
{
    GeoChemPhases phase_data;
    const std::map<std::string, std::vector<std::string>> db_changes_made_by_user;
    std::map<GeochemicalPhaseType, GeoChemPhaseData> phases_to_be_used;
    std::vector<std::string> user_order = {"H", "Ca", "NO3", "K", "Na", "Cl", "X"};
    phase_data.resetFromJson(json_input2);

    const int max_size = phase_data.selectPhases(phase_data.getPhaseNamesAndTypes(), phases_to_be_used);

    REQUIRE(max_size == 7); // unique minerals and ions

    auto ICS_ptr = InitChem::create_from_input_data
        (
            nullptr,
            db_changes_made_by_user,
            phases_to_be_used,
            max_size,
            "TESTING",
            user_order
        );

    int i = 0;
    for( auto bb: user_order)
    {
        int pos_in_db = ICS_ptr->SM_basis_->get_row_index(bb);
        REQUIRE( pos_in_db > -1); //exists?
        REQUIRE( ICS_ptr->pos_[i] == pos_in_db); // gives relative position from input to database ordering
        i++;
    }
}

TEST_CASE("Test chemistry equal when using different orders")
{
    GeoChemPhases phase_data;
    const std::map<std::string, std::vector<std::string>> db_changes_made_by_user;
    std::map<GeochemicalPhaseType, GeoChemPhaseData> phases_to_be_used;
    std::vector<std::string> user_order = {"H", "Ca", "NO3", "K", "Na", "Cl", "X"};
    phase_data.resetFromJson(json_input2);

    const int max_size = phase_data.selectPhases(phase_data.getPhaseNamesAndTypes(), phases_to_be_used);

    REQUIRE(max_size == 7); // unique minerals and ions

    auto ICS_ptr_unordered = InitChem::create_from_input_data
        (
            nullptr,
            db_changes_made_by_user,
            phases_to_be_used,
            max_size,
            "TESTING"
        );

    auto ICS_ptr = InitChem::create_from_input_data
        (
            nullptr,
            db_changes_made_by_user,
            phases_to_be_used,
            max_size,
            "TESTING",
            user_order
        );

    auto vchem_u = ICS_ptr_unordered->c_vchem_;
    auto ads_u   = ICS_ptr_unordered->c_ads_;
    auto mineral_u = ICS_ptr_unordered->c_mineral_;
    auto log_a_u = ICS_ptr_unordered->log_a_mineral_;
    auto pos_u = ICS_ptr_unordered->pos_;

    auto vchem = ICS_ptr->c_vchem_;
    auto ads   = ICS_ptr->c_ads_;
    auto mineral = ICS_ptr->c_mineral_;
    auto log_a = ICS_ptr->log_a_mineral_;
    auto pos = ICS_ptr->pos_;

    ICS_ptr->PRINT_DEBUG_CHEM_ = 2;

    CGeoChemIF IF_unordered(std::move(ICS_ptr_unordered));
    CGeoChemIF IF(std::move(ICS_ptr));

    double Temp = 273.15 + 25, pHi=7, Pres = 10e5, porosity = 1, frac_DL = 0, specific_surface_area = 1, sigma =0, psi =0;

    std::array<double, 3> mass_phase_{ 1.0, 1.0, 1.0 }; // Sw So Sg
    double Sw = 1;
    double f_DL = 0; // Diffusive Layer Thickness
    IF.SetSurfaceConc(Sw, vchem, f_DL);
    IF.SolveChem_I(vchem.data(),
                ads.data(),
                mineral.data(),
                log_a.data(),
                Temp,
                Pres,
                porosity,
        /* dt= */ 0.0, // doing equilibrium calculations
                specific_surface_area,
                frac_DL,
                mass_phase_,
                pHi,
                sigma,
                psi);

    double pHi_u = 7;
    IF_unordered.SetSurfaceConc(Sw, vchem_u, f_DL);
    IF_unordered.SolveChem_I(vchem_u.data(),
                ads_u.data(),
                mineral_u.data(),
                log_a_u.data(),
                Temp,
                Pres,
                porosity,
        /* dt= */ 0.0, // doing equilibrium calculations
                specific_surface_area,
                frac_DL,
                mass_phase_,
                pHi_u,
                sigma,
                psi);

    double abs_tol = 1e-5;
    CHECK_THAT(pHi, Catch::Matchers::WithinAbs(pHi_u, abs_tol));

    auto max_size_it = std::max_element(pos.begin(), pos.end());
    std::size_t size_db = static_cast<std::size_t>(*max_size_it);

    std::vector<double> cout(size_db+1, 0.0), cout_u(size_db+1, 0.0), ads_out(size_db+1, 0.0), ads_out_u(size_db+1, 0.0);

    std::vector<double> result, result_ordered;
    for (int i = 0; i < pos.size(); ++i) // map to database ordering
    {
        if(pos[i] > -1)
        {
            cout[pos[i]] = vchem[i];
            ads_out[pos[i]] = ads[i];
        }
        if(pos_u[i] > -1)
        {
            cout_u[pos_u[i]] = vchem_u[i];
            ads_out_u[pos_u[i]] = ads_u[i];
        }
    }
    for (int i = 0; i < pos.size(); ++i) // check if equal
    {
        CHECK_THAT(cout[i], Catch::Matchers::WithinAbs(cout_u[i], abs_tol));
        CHECK_THAT(ads_out[i], Catch::Matchers::WithinAbs(ads_out_u[i], abs_tol));
    }
}

TEST_CASE("Test initialization of ICS_full with minerals")
{
    GeoChemPhases phase_data;
    const std::map<std::string, std::vector<std::string>> db_changes_made_by_user;
    std::map<GeochemicalPhaseType, GeoChemPhaseData> phases_to_be_used;

    phase_data.resetFromJson(json_input);

    const int max_size = phase_data.selectPhases(phase_data.getPhaseNamesAndTypes(), phases_to_be_used);

    REQUIRE(max_size == 5); // unique minerals and ions

    auto ICS_ptr = InitChem::create_from_input_data
        (
            nullptr,
            db_changes_made_by_user,
            phases_to_be_used,
            max_size,
            "TESTING"
        );

    const std::vector<std::string> basis_present = {"Ca+2", "H+", "Cl-", "H2O"};
    const std::vector<std::string> phases_present = {"CALCITE", "CO2,g"};

    for( auto bb: basis_present)
    {
        REQUIRE( ICS_ptr->SM_basis_->get_row_index(bb) > -1);
    }

    for( auto bb: phases_present)
    {
        REQUIRE( ICS_ptr->SM_mineral_->get_row_index(bb) > -1);
    }
}

TEST_CASE("Test solve from JSON and interpolate")
{
    GeoChemPhases phase_data;
    const std::map<std::string, std::vector<std::string>> db_changes_made_by_user;
    std::map<GeochemicalPhaseType, GeoChemPhaseData> phases_to_be_used;
    const std::vector<std::string> species_names_in_order;

    phase_data.resetFromJson(json_input);

    const int max_size = phase_data.selectPhases(phase_data.getPhaseNamesAndTypes(), phases_to_be_used);

    REQUIRE(max_size == 5); // unique minerals and ions

    auto ICS_ptr = InitChem::create_from_input_data
        (
            nullptr,
            db_changes_made_by_user,
            phases_to_be_used,
            max_size,
            "TESTING",
            species_names_in_order,
            /* INTERPOLATE = */ 1
        );
    auto vchem = ICS_ptr->c_vchem_;
    auto ads   = ICS_ptr->c_ads_;
    auto mineral = ICS_ptr->c_mineral_;
    auto log_a = ICS_ptr->log_a_mineral_;

    auto vchem2 = vchem;
    auto ads2   = ads;
    auto mineral2 = mineral;
    auto log_a2 = log_a;

    CGeoChemIF IF(std::move(ICS_ptr));

    double Temp = 273.15 + 25, pHi=7, Pres = 10e5, porosity = 1, frac_DL = 0, specific_surface_area = 1, sigma =0, psi =0;

    std::array<double, 3> mass_phase_{ 1.0, 1.0, 1.0 };
    IF.SolveChem_I(vchem.data(),
                ads.data(),
                mineral.data(),
                log_a.data(),
                Temp,
                Pres,
                porosity,
        /* dt= */ 0.0, // doing equilibrium calculations
                specific_surface_area,
                frac_DL,
                mass_phase_,
                pHi,
                sigma,
                psi);

    double abs_tol = 1e-5;
    CHECK_THAT(pHi, Catch::Matchers::WithinAbs(pH_to_match_json, abs_tol));

    // check if we call the solver once more using the same input, it does not invoke more calls, i.e. it is using splay tree look up
    pHi = 7; // reset pH
    IF.SolveChem_I(vchem2.data(),
                ads2.data(),
                mineral2.data(),
                log_a2.data(),
                Temp,
                Pres,
                porosity,
        /* dt= */ 0.0, // doing equilibrium calculations
                specific_surface_area,
                frac_DL,
                mass_phase_,
                pHi,
                sigma,
                psi);
    CHECK_THAT(pHi, Catch::Matchers::WithinAbs(pH_to_match_json, abs_tol)); // pH should be the same as before
    REQUIRE(IF.getSolverCalls() == 1);
}

TEST_CASE("Database modification with JSON")
{
    // Set up phases
    GeoChemPhases phase_data;
    phase_data.resetFromJson(json_db);
    std::map<GeochemicalPhaseType, GeoChemPhaseData> phases_to_be_used;
    const int max_size = phase_data.selectPhases(phase_data.getPhaseNamesAndTypes(), phases_to_be_used);

    // Use test class that inherits CGeoChemIF just to setup database modification from JSON
    TestGeoChemIF geoChemIF(json_db);
    auto db_mod = geoChemIF.getDbModifications();

    // Check modifications
    for (const auto& section : GeochemicalDatabaseKeyword::ALL_DATABASE_KEYWORDS) {
        // Skip SURFACE_SPECIES for now
        if (section == GeochemicalDatabaseKeyword::SURFACE_SPECIES) {
            continue;
        }
        REQUIRE(db_mod.find(section) != db_mod.end());
    }

    // Init. chemistry to internalize database modifications
    auto ICS_ptr = InitChem::create_from_input_data(nullptr,
                                                    db_mod,
                                                    phases_to_be_used,
                                                    max_size,
                                                    "TESTING");

    // Check internal database for modifications
    const std::vector<std::string> basis_species_expected = {"Fee+2", "Hdg", "Hdg2", "Z+"};
    for (const auto& elem : basis_species_expected) {
        REQUIRE(ICS_ptr->SM_basis_->get_row_index(elem) > -1);
    }

    const std::vector<std::string> minerals_expected = {"Hdg2(g)"};
    for (const auto& elem : minerals_expected) {
        REQUIRE(ICS_ptr->SM_mineral_->get_row_index(elem) > -1);
    }

    const std::vector<std::string> all_species_expected = {"ZNO3"};
    for (const auto& elem : all_species_expected) {
        REQUIRE(ICS_ptr->SM_all_->get_row_index(elem) > -1);
    }

    const std::vector<std::string> removed_species = {"HNO3", "CaOH+"};
    for (const auto& elem : removed_species) {
        REQUIRE(ICS_ptr->SM_basis_->get_row_index(elem) == -1);
    }
}
