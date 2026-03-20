#if __has_include(<catch2/catch.hpp>)
#include <catch2/catch.hpp>
#else
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#endif

#include <map>
#include <string>
#include <fstream>

#include <opm/simulators/geochemistry/Thermo/hkf.h>
#include <opm/simulators/geochemistry/Thermo/eps_JN.h>
#include <opm/simulators/geochemistry/Thermo/ions.h>
#include <opm/simulators/geochemistry/Thermo/thermodata.h>
#include <opm/simulators/geochemistry/Thermo/water.h>

static constexpr double abs_tolerance_ = 1.0e-7;
static constexpr double abs_tolerance2_ = 1.0e-4;  // for most water props.

/*
* Stores reference cases 1, 2, and 3 from Table 5 of the IAPWS-97
* industrial thermodynamic formulation.
*
* See also:
*
*   - https://github.com/jjgomera/iapws.git
*   - http://twt.mpei.ac.ru/mcs/worksheets/iapws/IAPWS-IF97-Region1.xmcd
*   - http://www.casprod.eu/~leon/research/iapws-if97/listing.html
*/
struct WaterPropertiesIAPWS_97
{
    WaterPropertiesIAPWS_97(int reference_case)
    {
        switch (reference_case) {
            case 1:
            {
                pressure_ = 3.0e6;
                temperature_ = 300.0;
                v_ = 0.100215168e-2;
                h_ = 0.115331273e3;
                u_ = 0.112324818e3;
                s_ = 0.392294792;
                cp_ = 0.417301218e1;
                w_ = 0.150773921e4;
                break;
            }
            case 2:
            {
                pressure_ = 80.0e6;
                temperature_ = 300.0;
                v_ = 0.971180894e-3;
                h_ = 0.184142828e3;
                u_ = 0.106448356e3;
                s_ = 0.368563852;
                cp_ = 0.401008987e1;
                w_ = 0.163469054e4;
                break;
            }
            case 3:
            {
                pressure_ = 3.0e6;
                temperature_ = 500.0;
                v_ = 0.120241800e-2;
                h_ = 0.975542239e3;
                u_ = 0.971934985e3;
                s_ = 0.258041912e1;
                cp_ = 0.465580682e1;
                w_ = 0.124071337e4;
                break;
            }
        }
    }

    double pressure_;
    double temperature_;

    double v_;  // specific volume (m3/kg)
    double h_;  // specific enthalpy (kJ/kg)
    double u_;  // specific internal energy (kJ/kg)
    double s_;  // specific entropy (kJ/kg/K)
    double cp_;  // specific isobaric heat capacity (kJ/kg/K)
    double w_;  // speed of sound (m/s)

    //double cv_;  // specific isochoric heat capacity (kJ/kg/K)
};


TEST_CASE("Test calculation of thermodynamic properties for water")
{

    water water_props;

    static const std::array<WaterPropertiesIAPWS_97, 3> reference_cases =
    {
        WaterPropertiesIAPWS_97(1),
        WaterPropertiesIAPWS_97(2),
        WaterPropertiesIAPWS_97(3),
    };

    for(const auto& ref_props: reference_cases)
    {
        const double T = ref_props.temperature_;
        const double P = ref_props.pressure_;

        water_props.gibbsIAPWS(T, P);

        CHECK_THAT(water_props.v_, Catch::Matchers::WithinAbs(ref_props.v_, abs_tolerance_));

        // Note: Reference results are in kJ/kg or kJ/kg/K, we use Joules...
        CHECK_THAT(1.0e-3*water_props.h_, Catch::Matchers::WithinAbs(ref_props.h_, abs_tolerance2_));
        CHECK_THAT(1.0e-3*water_props.u_, Catch::Matchers::WithinAbs(ref_props.u_, abs_tolerance2_));
        CHECK_THAT(1.0e-3*water_props.s_, Catch::Matchers::WithinAbs(ref_props.s_, abs_tolerance2_));
        CHECK_THAT(1.0e-3*water_props.cp_, Catch::Matchers::WithinAbs(ref_props.cp_, abs_tolerance2_));
        CHECK_THAT(water_props.w_, Catch::Matchers::WithinAbs(ref_props.w_, abs_tolerance2_));

        water_props.printProperties();
    }

}

TEST_CASE("Test calculation of saturation pressure for water")
{
    water water_props;

    // Table 35 of IAPWS97 paper
    static constexpr std::array<std::pair<double, double>, 3> table35_tests =
    {
        std::make_pair<double, double>(300.0, 0.353658941e-2),
        std::make_pair<double, double>(500.0, 0.263889776e1),
        std::make_pair<double, double>(600.0, 0.123443146e2)
    };

    for(const auto& [T, Psat] : table35_tests)
    {
        const double satPress = 1.0e-6*water_props.PsatIAPWS(T);
        CHECK_THAT(satPress, Catch::Matchers::WithinAbs(Psat, abs_tolerance_));
    }
}
