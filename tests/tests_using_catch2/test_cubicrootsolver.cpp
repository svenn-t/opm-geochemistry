/*
Test CubicRootSolver on simple, known solutions to cubic equations
*/
#include <catch2/catch.hpp>

#include <opm/simulators/geochemistry/Numerical/RootSolvers.hpp>

#include <vector>

using CubicRoot = CubicRootSolver<double>;

TEST_CASE( "Cubic equation 1, three real roots" ) {
    // Setup equation parameters
    const double a = 1;
    const double b = -7;
    const double c = 14;
    const double d = -8;

    // Expected roots (in descending order)
    std::vector<double> expected = {4, 2, 1}; 

    // Calculated roots
    auto roots = CubicRoot::solve(a, b, c, d);

    // Loop over roots and check if they are correct
    for (std::size_t i = 0; i < roots.size(); ++i) {
        CHECK_THAT(roots[i], Catch::Matchers::WithinAbs(expected[i], 1.0e-10));
    }
}

TEST_CASE( "Cubic equation 2, one real root" ) {
    // Setup equation parameters
    const double a = 1;
    const double b = -1;
    const double c = -1;
    const double d = -2;

    // Expected roots (in descending order)
    std::vector<double> expected = {2}; 

    // Calculated roots
    auto roots = CubicRoot::solve(a, b, c, d);

    // Loop over roots and check if they are correct
    for (std::size_t i = 0; i < roots.size(); ++i) {
        CHECK_THAT(roots[i], Catch::Matchers::WithinAbs(expected[i], 1.0e-10));
    }
}

TEST_CASE( "Cubic equation 3, two equal and one distinct root" ) {
    // Setup equation parameters
    const double a = 1;
    const double b = 0;
    const double c = -3;
    const double d = 2;

    // Expected roots (in descending order)
    std::vector<double> expected = {1, -2}; 

    // Calculated roots
    auto roots = CubicRoot::solve(a, b, c, d);

    // Loop over roots and check if they are correct
    for (std::size_t i = 0; i < roots.size(); ++i) {
        CHECK_THAT(roots[i], Catch::Matchers::WithinAbs(expected[i], 1.0e-10));
    }
}