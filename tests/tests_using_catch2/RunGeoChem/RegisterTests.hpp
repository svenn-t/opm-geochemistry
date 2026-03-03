#ifndef GEO_CHEM_REGISTER_TESTS_DEF
#define GEO_CHEM_REGISTER_TESTS_DEF

#include <cstddef>
#include <memory>
#include <string_view>
#include <vector>

#include "GeoChemTestCase.hpp"

class GeoChemTestCaseFactory
{

public:

    GeoChemTestCaseFactory();

    const GeoChemTestCase* getTestCase(std::string_view name) const;

    const std::vector<GeoChemTestCasePtr>& getAllTestCases() const;

private:
    enum class Type{ EqSolver, OneDTransport };

    std::vector<GeoChemTestCasePtr> all_tests_;
};

#endif
