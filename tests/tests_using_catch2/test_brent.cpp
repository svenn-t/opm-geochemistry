/*
*   Test that Brent's method for finding the zero of a real-valued function
*   of one variable yields acceptable results on a selection of very simple
*   functions.
*/
#include <catch2/catch.hpp>

#include <cmath>
#include <vector>

#include <opm/simulators/geochemistry/Numerical/BrentMethod.hpp>

namespace test_brent_method
{  // to avoid potential name clashes

double x_squared(double x){ return x*x; }

}

using namespace test_brent_method;

using BrentMethod = BrentMethodReal<double>;

TEST_CASE( "Solve x^2 - 2 = 0") {

    auto f = [](double x) { return x_squared(x) - 2.0; };
    const auto expected = std::sqrt(2.0);

    int noIter = 0;
    const auto calculated = BrentMethod::solve(f, 0.0, 2.0, 30, 1.0e-8, noIter);

    CHECK_THAT(calculated, Catch::Matchers::WithinAbs(expected, 1.0e-10));

}
