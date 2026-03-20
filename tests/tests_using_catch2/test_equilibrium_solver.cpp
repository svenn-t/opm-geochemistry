#if __has_include(<catch2/catch.hpp>)
#include <catch2/catch.hpp>
#else
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#endif

#include <fstream>
#include <string>
#include <filesystem>
#include <opm/simulators/geochemistry/GeoChemIF.h>

#include "reference_results_eqsolver/EquilibriumTestCases/SetupAndRunEqsolver.hpp"
#include "reference_results_eqsolver/EquilibriumTestCases/RegisterEquilibriumTests.hpp"


/* Reads reference simulation data from a binary file. */
BasVecInfo getReferenceResults(GeoChemTestCase* testcase)
{
    const std::string archive_file = testcase->name() + ".cereal";

    BasVecInfo info_backup;
    std::ifstream cereal_infile(archive_file, std::ios::binary);
    cereal::BinaryInputArchive archive_infile(cereal_infile);
    archive_infile(info_backup);

    return info_backup;
}


TEST_CASE( "Test equilibrium solver") {

    static constexpr double abs_tolerance = 1.0e-5;
    // static constexpr double rel_tolerance = 1.0e-3;

    auto testFactory = EquilibriumTestCaseFactory();
    const auto& all_tests = testFactory.getAllEquilibriumTestCases();

    for(std::size_t i=0; i < all_tests.size(); ++i)
    {
        GeoChemTestCase* current_test = all_tests[i].get();

        std::cout << current_test->name() << "\n";
        std::filesystem::path cwd = std::filesystem::current_path();
        std::cout << "Current working directory: " << cwd << '\n';
        const auto reference_results = getReferenceResults(current_test);
        const auto eq_sol_results = gen_results_from_equilibrium_calculation(current_test);

        const auto reference_solution_props = reference_results.key_solution_properties_;
        for(auto const& [key, value]: eq_sol_results.key_solution_properties_){
            auto it = reference_solution_props.find(key);
            if(it != reference_solution_props.end()){
                const auto expected_value = it->second;
                CHECK_THAT(value, Catch::Matchers::WithinAbs(expected_value, abs_tolerance));
                // CHECK_THAT(value, Catch::Matchers::WithinRel(expected_value, rel_tolerance));
            }
        }

        const auto reference_species_data = reference_results.species_properties_;
        for(auto const& [key, value]: eq_sol_results.species_properties_){
            auto it = reference_species_data.find(key);
            if(it != reference_species_data.end()){
                const auto expected_value = it->second;
                if(std::fabs(value-expected_value) > abs_tolerance){
                    std::cout << "Key=(" << key.first << ", " << key.second << ")";
                    std::cout << ", Value=" << value << ", expected=" << expected_value << ".\n";
                }
                CHECK_THAT(value, Catch::Matchers::WithinAbs(expected_value, abs_tolerance));
                // CHECK_THAT(value, Catch::Matchers::WithinRel(expected_value, rel_tolerance));
            }
        }
    }

}
