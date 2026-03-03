#include "RegisterEquilibriumTests.hpp"

// Regiser any new tests here...
#include "Test1.hpp"
#include "Test1_DL.hpp"
#include "Test1_IO_DL.hpp"
#include "Test1_IO_equil.hpp"
#include "Test1_IO_only.hpp"
#include "Test1_Redox.hpp"
#include "Test2_sc1.hpp"
#include "Test2_sc1_DL.hpp"
#include "Test_with_kinetics.hpp"
#include "Test_add_mineral_phase.hpp"

EquilibriumTestCaseFactory::EquilibriumTestCaseFactory()
: all_tests_()
{
    // Add new tests here...

    // Tests where we solve wrt. TOTH+
    all_tests_.emplace_back(std::make_unique<Test1>());
    all_tests_.emplace_back(std::make_unique<Test1_DL>());
    all_tests_.emplace_back(std::make_unique<Test1_IO_DL>());
    all_tests_.emplace_back(std::make_unique<Test1_IO_equil>());
    all_tests_.emplace_back(std::make_unique<Test1_IO_only>());

    // Tests where we compute pH from charge-balance
    all_tests_.emplace_back(std::make_unique<Test2_sc1>());
    all_tests_.emplace_back(std::make_unique<Test2_sc1_DL>());
    all_tests_.emplace_back(std::make_unique<TestWithKinetics>());
    all_tests_.emplace_back(std::make_unique<Test_add_mineral_phase>());

    // Redox
    all_tests_.emplace_back(std::make_unique<Test1_Redox>());
}

const std::vector<GeoChemTestCasePtr>& EquilibriumTestCaseFactory::getAllEquilibriumTestCases() const
{
    return all_tests_;
}
