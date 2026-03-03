#ifndef EQSOLVER_TEST_ADD_MINERAL_PHASE_DEF
#define EQSOLVER_TEST_ADD_MINERAL_PHASE_DEF

#include "GeoChemTestCase.hpp"

struct Test_add_mineral_phase: public GeoChemTestCase
{
    std::string name() const override
    {
        return "test_add_mineral_phase";
    }

    std::string getInputFileAsString() const override
    {
        return R"""(Temp 70.0
Pres 8.0e5

equilibrate 1

# Because the reference simulation used a different set of HKF parameters for the GCO3H^0 complex
SECONDARY_SPECIES
GCO3H = GCO3- - H2O + H+ / HKF  -92250  -98900  28.1  15.2964  -10.0382  -55.4493  7.4092  36.8069  3.5851  -0.3107 /
#GCO3H = GCO3- + H+ / HKF  -146940.8956 -98900  28.1    6.2466  7.4711  2.8136  -3.0879 38.4529 5.3534  -0.1934 /
/end

MINERAL_PHASES
WITHERITE = Ba+2 + HCO3- - H+ / HKF 45.81 -278400  -297500  26.8  21.5  11.06  -3.91 /
bCALCITE = Ca+2 + HCO3- - H+ / ANA 36.93 1.84867 /
/end

geochem

solution 0
pH 7 charge
Ba 0.005
K 0.01
Cl 0.02
Ca 1.0e-8
HCO3 1.0e-8
/end

equilibrium_phases
WITHERITE 1 0
bCALCITE 1 0
CO2(G) 1 -3.5
/end

iexchange
X 1.0e-2
/ end

/end)""";

    }
};

#endif
