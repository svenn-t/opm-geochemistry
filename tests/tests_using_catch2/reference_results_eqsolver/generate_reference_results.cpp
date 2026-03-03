#include <memory>
#include <vector>

#include <opm/simulators/geochemistry/GeoChemIF.h>

#include "EquilibriumTestCases/SetupAndRunEqsolver.hpp"
#include "EquilibriumTestCases/RegisterEquilibriumTests.hpp"

int main(){

    auto testFactory = EquilibriumTestCaseFactory();
    const auto& all_tests = testFactory.getAllEquilibriumTestCases();

    for(std::size_t i=0; i < all_tests.size(); ++i)
    {
        GeoChemTestCase* current_test = all_tests[i].get();
        gen_results_from_equilibrium_calculation(current_test, /* serialize= */ true);
    }

    return 0;
}
