#include <catch2/catch.hpp>

#include <cmath>

#include <opm/simulators/geochemistry/Extra/DiffuseLayerIntegrand.hpp>
#include <opm/simulators/geochemistry/Numerical/GaussLegendre.hpp>

static constexpr double abs_tolerance = 1.0e-10;
static constexpr double rel_tolerance = 1.0e-10;
static constexpr auto integrator = GaussLegendreIntegrator<double, 10>();
static constexpr double default_E0_ = 0.38826992498273194;

struct ExampleTestCase{

    ExampleTestCase(const std::string& name, double n1=0.0, double n2=0.0, double n_2=0.0, double E0=default_E0_)
    : name_(name)
    , n1_(n1)
    , n2_(n2)
    , n_2_(n_2)
    , E0_(E0)
    {
    }

    std::string name_;
    double n1_;
    double n2_;
    double n_2_;
    double E0_;

};

static auto get_all_test_cases() -> std::vector<ExampleTestCase>{

    std::vector<ExampleTestCase> all_test_cases;
    all_test_cases.emplace_back(ExampleTestCase("Pure NaCl", 0.219, 0.0, 0.0));  // only monovalents
    all_test_cases.emplace_back(ExampleTestCase("Pure Na2SO4", 0.438, 0.0, 0.219));  // no +2 ions
    all_test_cases.emplace_back(ExampleTestCase("Pure CaCl2", 0.0, 0.219, 0.0));  // no -2 ions
    all_test_cases.emplace_back(ExampleTestCase("Pure MgSO4", 0.0, 0.219, 0.219));  // no monovalents
    all_test_cases.emplace_back(ExampleTestCase("Negative discriminant example case", 0.00041550452510550571, 0.037628703995817941, 0.03812007099563542));

    return all_test_cases;
}

TEST_CASE("Compare numerical vs. analytical integration for surface excess integrals"){

    const auto all_test_cases = get_all_test_cases();

    // Note : It seems we cannot use variable name "case", name clash??
    for (const auto &tcase : all_test_cases)
    {
        const auto E0 = tcase.E0_;

        double numerical_integral = 0.0;
        double analytical_integral = 0.0;

        SECTION("Charge +1")
        {
            auto f = DiffuseLayerIntegrandOne(tcase.n1_, tcase.n2_, tcase.n_2_);
            analytical_integral = f.integrate(1.0, 1.0/E0);
            numerical_integral = integrator.integrate(f, 1.0, 1.0/E0);
        }

        SECTION("Charge -1")
        {
            auto f = DiffuseLayerIntegrandMinusOne(tcase.n1_, tcase.n2_, tcase.n_2_);
            analytical_integral = f.integrate(1.0, 1.0/E0);
            numerical_integral = integrator.integrate(f, 1.0, 1.0/E0);
        }

        SECTION("Charge +2")
        {
            auto f = DiffuseLayerIntegrandTwo(tcase.n1_, tcase.n2_, tcase.n_2_);
            analytical_integral = f.integrate(1.0, 1.0/E0);
            numerical_integral = integrator.integrate(f, 1.0, 1.0/E0);
        }

        SECTION("Charge -2")
        {
            auto f = DiffuseLayerIntegrandMinusTwo(tcase.n1_, tcase.n2_, tcase.n_2_);
            analytical_integral = f.integrate(1.0, 1.0/E0);
            numerical_integral = integrator.integrate(f, 1.0, 1.0/E0);
        }

        CHECK_THAT(numerical_integral, Catch::Matchers::WithinAbs(analytical_integral, abs_tolerance));
        CHECK_THAT(numerical_integral, Catch::Matchers::WithinRel(analytical_integral, rel_tolerance));
    }
}
