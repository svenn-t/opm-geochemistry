#include "RegisterTests.hpp"

#include <algorithm>

// Register any new tests here (and in the constructor)...
#include "cases_transport/Test_caco3_to_mgco3.hpp"
#include "cases_transport/Test_rate_and_equilibrium_minerals.hpp"
#include "cases_transport/Test_coreflood_interpolation.hpp"
#include "cases_transport/Test_InterpolationWithMineralsButNoSurfSpecies.hpp"
#include "cases_transport/Test_mineral_activation_energy.hpp"

GeoChemTestCaseFactory::GeoChemTestCaseFactory()
    : all_tests_()
{
    all_tests_.emplace_back(std::make_unique<Test_Calcite_to_Magnesite>());
    all_tests_.emplace_back(std::make_unique<Test_rate_and_equilibrium_minerals>());
    all_tests_.emplace_back(std::make_unique<Test_CorefloodInterpolation>());
    all_tests_.emplace_back(std::make_unique<Test_InterpolationWithMineralsButNoSurfSpecies>());
    all_tests_.emplace_back(std::make_unique<Test_MineralActivationEnergy>());
}

const GeoChemTestCase* GeoChemTestCaseFactory::getTestCase(std::string_view name) const
{
    auto result = std::find_if(
        std::begin(all_tests_),
        std::end(all_tests_),
        [name](const GeoChemTestCasePtr& test_case){ return test_case->name() == name; }
    );

    if(result != std::end(all_tests_))
    {
        return &*result->get();
    }
    return nullptr;
}

const std::vector<GeoChemTestCasePtr>& GeoChemTestCaseFactory::getAllTestCases() const
{
    return all_tests_;
}
