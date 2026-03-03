/*
*   Test that Gauss-Legendre integration yields acceptable results for a
*   selection of very simple functions.
*/
#include <catch2/catch.hpp>

#include <cmath>
#include <vector>

#include <opm/simulators/geochemistry/Numerical/GaussLegendre.hpp>

static constexpr auto integrator2 = GaussLegendreIntegrator<double, 2>();
static constexpr auto integrator3 = GaussLegendreIntegrator<double, 3>();
static constexpr auto integrator4 = GaussLegendreIntegrator<double, 4>();
static constexpr auto integrator5 = GaussLegendreIntegrator<double, 5>();
static constexpr auto integrator6 = GaussLegendreIntegrator<double, 6>();
static constexpr auto integrator7 = GaussLegendreIntegrator<double, 7>();
static constexpr auto integrator8 = GaussLegendreIntegrator<double, 8>();
static constexpr auto integrator9 = GaussLegendreIntegrator<double, 9>();
static constexpr auto integrator10 = GaussLegendreIntegrator<double, 10>();


namespace test_gauss_legendre
{  // to avoid name clashes if defining the same functions elsewhere...

double x_squared(double x){
    return x*x;
}

double one_over_x(double x){
    return 1.0/x;
}

double exponential_function(double x){
    return std::exp(x);
}

double sine_function(double x){
    return std::sin(x);
}

}

using namespace test_gauss_legendre;


TEST_CASE( "Integrate x^2") {

    static constexpr double a = 1.0;
    static constexpr double b = 4.0;
    static const double analytical_integral = (1.0/3)*(std::pow(b, 3) - std::pow(a, 3));

    double numerical_integral;

    SECTION("2-point quadrature"){
        numerical_integral = integrator2.integrate(x_squared, a, b);
    }

    SECTION("3-point quadrature"){
        numerical_integral = integrator3.integrate(x_squared, a, b);
    }

    SECTION("4-point quadrature"){
        numerical_integral = integrator4.integrate(x_squared, a, b);
    }

    SECTION("5-point quadrature"){
        numerical_integral = integrator5.integrate(x_squared, a, b);
    }

    SECTION("6-point quadrature"){
        numerical_integral = integrator6.integrate(x_squared, a, b);
    }

    SECTION("7-point quadrature"){
        numerical_integral = integrator7.integrate(x_squared, a, b);
    }

    SECTION("8-point quadrature"){
        numerical_integral = integrator8.integrate(x_squared, a, b);
    }

    SECTION("9-point quadrature"){
        numerical_integral = integrator9.integrate(x_squared, a, b);
    }

    SECTION("10-point quadrature"){
        numerical_integral = integrator10.integrate(x_squared, a, b);
    }

    CHECK_THAT(numerical_integral, Catch::Matchers::WithinAbs(analytical_integral, 1.0e-10));
    CHECK_THAT(numerical_integral, Catch::Matchers::WithinRel(analytical_integral, 1.0e-10));
}

TEST_CASE( "Integrate 1/x") {

    static constexpr double a = 1.0;
    static constexpr double b = 4.0;
    static const double analytical_integral = std::log(b/a);

    double numerical_integral;
    // Use high tolerances by default
    double abs_tol = 0.1;
    double rel_tol = 0.1;

    SECTION("2-point quadrature"){
        numerical_integral = integrator2.integrate(one_over_x, a, b);
    }

    SECTION("3-point quadrature"){
        numerical_integral = integrator3.integrate(one_over_x, a, b);
    }

    SECTION("4-point quadrature"){
        numerical_integral = integrator4.integrate(one_over_x, a, b);
    }

    SECTION("5-point quadrature"){
        numerical_integral = integrator5.integrate(one_over_x, a, b);
    }

    SECTION("6-point quadrature"){
        numerical_integral = integrator6.integrate(one_over_x, a, b);
    }

    SECTION("7-point quadrature"){
        numerical_integral = integrator7.integrate(one_over_x, a, b);
    }

    SECTION("8-point quadrature"){
        numerical_integral = integrator8.integrate(one_over_x, a, b);
    }

    SECTION("9-point quadrature"){
        numerical_integral = integrator9.integrate(one_over_x, a, b);
    }

    SECTION("10-point quadrature"){
        numerical_integral = integrator10.integrate(one_over_x, a, b);
        // Require much more precision for 10-point quadrature
        abs_tol = 1.0e-8;
        rel_tol = 1.0e-8;
    }

    CHECK_THAT(numerical_integral, Catch::Matchers::WithinAbs(analytical_integral, abs_tol));
    CHECK_THAT(numerical_integral, Catch::Matchers::WithinRel(analytical_integral, rel_tol));
}

TEST_CASE( "Integrate e^x") {

    static constexpr double a = 1.0;
    static constexpr double b = 2.0;
    static const double analytical_integral = std::exp(b) - std::exp(a);

    double numerical_integral;
    // Use high tolerances by default
    double abs_tol = 0.1;
    double rel_tol = 0.1;

    SECTION("2-point quadrature"){
        numerical_integral = integrator2.integrate(exponential_function, a, b);
    }

    SECTION("3-point quadrature"){
        numerical_integral = integrator3.integrate(exponential_function, a, b);
    }

    SECTION("4-point quadrature"){
        numerical_integral = integrator4.integrate(exponential_function, a, b);
    }

    SECTION("5-point quadrature"){
        numerical_integral = integrator5.integrate(exponential_function, a, b);
    }

    SECTION("6-point quadrature"){
        numerical_integral = integrator6.integrate(exponential_function, a, b);
    }

    SECTION("7-point quadrature"){
        numerical_integral = integrator7.integrate(exponential_function, a, b);
    }

    SECTION("8-point quadrature"){
        numerical_integral = integrator8.integrate(exponential_function, a, b);
    }

    SECTION("9-point quadrature"){
        numerical_integral = integrator9.integrate(exponential_function, a, b);
    }

    SECTION("10-point quadrature"){
        numerical_integral = integrator10.integrate(exponential_function, a, b);
        // Require much more precision for 10-point quadrature
        abs_tol = 1.0e-10;
        rel_tol = 1.0e-10;
    }

    CHECK_THAT(numerical_integral, Catch::Matchers::WithinAbs(analytical_integral, abs_tol));
    CHECK_THAT(numerical_integral, Catch::Matchers::WithinRel(analytical_integral, rel_tol));
}

TEST_CASE( "Integrate sin(x)") {

    static constexpr double a = 1.0;
    static constexpr double b = 10.0;
    static const double analytical_integral = std::cos(a) - std::cos(b);

    double numerical_integral;
    // Use high tolerances by default
    double abs_tol = 0.1;
    double rel_tol = 0.1;

    // 2-point & 3-point quadrature is not expected to be good in this case..

    SECTION("4-point quadrature"){
        numerical_integral = integrator4.integrate(sine_function, a, b);
    }

    SECTION("5-point quadrature"){
        numerical_integral = integrator5.integrate(sine_function, a, b);
    }

    SECTION("6-point quadrature"){
        numerical_integral = integrator6.integrate(sine_function, a, b);
    }

    SECTION("7-point quadrature"){
        numerical_integral = integrator7.integrate(sine_function, a, b);
    }

    SECTION("8-point quadrature"){
        numerical_integral = integrator8.integrate(sine_function, a, b);
    }

    SECTION("9-point quadrature"){
        numerical_integral = integrator9.integrate(sine_function, a, b);
    }

    SECTION("10-point quadrature"){
        numerical_integral = integrator10.integrate(sine_function, a, b);
        // Require much more precision for 10-point quadrature
        abs_tol = 1.0e-10;
        rel_tol = 1.0e-10;
    }

    CHECK_THAT(numerical_integral, Catch::Matchers::WithinAbs(analytical_integral, abs_tol));
    CHECK_THAT(numerical_integral, Catch::Matchers::WithinRel(analytical_integral, rel_tol));
}
