#ifndef EQSOLVER_REGISTER_TESTS_DEF
#define EQSOLVER_REGISTER_TESTS_DEF

#include <cstddef>
#include <memory>
#include <vector>

#include "GeoChemTestCase.hpp"

class EquilibriumTestCaseFactory
{

public:

    EquilibriumTestCaseFactory();

    const std::vector<GeoChemTestCasePtr>& getAllEquilibriumTestCases() const;

private:
    std::vector<GeoChemTestCasePtr> all_tests_;
};

#endif
