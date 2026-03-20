#if __has_include(<catch2/catch.hpp>)
#include <catch2/catch.hpp>
#else
#include <catch2/catch_test_macros.hpp>
#endif

#include <sstream>
#include <string>
#include <iostream>
#include <opm/simulators/geochemistry/Core/GeoChemPhases.hpp>

#include "reference_results_eqsolver/EquilibriumTestCases/RegisterEquilibriumTests.hpp"

TEST_CASE("A few tests to check that the input for the reference test cases behave as expected...") 
{
    auto testFactory = EquilibriumTestCaseFactory();
    const auto& all_tests = testFactory.getAllEquilibriumTestCases();

    GeoChemPhases geochem_phases;
        
    for(auto& current_test: all_tests)
    {
        // Check that deleting all data from the object works...
        geochem_phases.clear();
        for(const auto& phaseType: GeoChemPhases::ALL_POSSIBLE_PHASES)
        {
            REQUIRE(geochem_phases.hasPhase(phaseType) == false);
            REQUIRE(geochem_phases.phaseSize(phaseType) == 0);
        }
        
        //...and finally read new data
        auto input = current_test->getInputFileAsString();
        std::istringstream inputStream(input);
        
        geochem_phases.updateFromInputStream(inputStream);
        
        REQUIRE(geochem_phases.isValid()); // there must be at least one solution
    }

}

