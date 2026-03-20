#include "ApprovalTests.hpp"
#if __has_include(<catch2/catch.hpp>)
#include <catch2/catch.hpp>
#else
#include <catch2/catch_test_macros.hpp>
#endif

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <opm/simulators/geochemistry/GeoChemIF.h>

#include "RunGeoChem/RegisterTests.hpp"
#include "RunGeoChem/SetupAndRunSolver.hpp"

void fOut_EffluentData(const std::vector<EffluentIonData>& effluentsForAllSteps, std::ostream& os)
{
    std::size_t timeIdx = 0;

    // Loop through each time step, and print the output
    for(const auto& effluents_t: effluentsForAllSteps)
    {
        os << "TimeIndex" << " = " << timeIdx << "\n";
        // Note: It is assumed that names_ & values_ have the same size.
        for(std::size_t i=0; i < effluents_t.names_.size(); ++i)
        {
            os << effluents_t.names_[i] << " = " << effluents_t.values_[i] << "\n";
        }
        ++timeIdx;
    }
}

TEST_CASE("Test 1D transport cases")
{
    GeoChemTestCaseFactory testFactory;

    SECTION("Calcite to Magnesite")
    {
        const GeoChemTestCase& test_case = *testFactory.getTestCase("Test_Calcite_to_Magnesite");
        auto results = gen_results_from_transport_calculation(test_case, /* serialize_flag= */ 1);
        ApprovalTests::Approvals::verify(results, fOut_EffluentData);
    }

    SECTION("test_rate_and_equilibrium_minerals")
    {
        const GeoChemTestCase& test_case = *testFactory.getTestCase("test_rate_and_equilibrium_minerals");
        auto results = gen_results_from_transport_calculation(test_case, /* serialize_flag= */ 1);
        ApprovalTests::Approvals::verify(results, fOut_EffluentData);
    }

    SECTION("Coreflood with interpolation")
    {
        const GeoChemTestCase& test_case = *testFactory.getTestCase("Coreflood_interpolation");
        auto results = gen_results_from_transport_calculation(test_case, /* serialize_flag= */ 1);
        ApprovalTests::Approvals::verify(results, fOut_EffluentData);
    }

    SECTION("InterpolationWithMineralsButNoSurfSpecies")
    {
        const GeoChemTestCase& test_case = *testFactory.getTestCase("InterpolationWithMineralsButNoSurfSpecies");
        auto results = gen_results_from_transport_calculation(test_case, /* serialize_flag= */ 1);
        ApprovalTests::Approvals::verify(results, fOut_EffluentData);
    }

    SECTION("Test_mineral_activation_energy")
    {
        const GeoChemTestCase& test_case = *testFactory.getTestCase("Test_mineral_activation_energy");
        auto results = gen_results_from_transport_calculation(test_case, /* serialize_flag= */ 1);
        ApprovalTests::Approvals::verify(results, fOut_EffluentData);
    }
}

