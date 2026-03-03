#ifndef EQSOLVER_TEST_WITH_KINETICS_DEF
#define EQSOLVER_TEST_WITH_KINETICS_DEF

#include "GeoChemTestCase.hpp"

struct TestWithKinetics: public GeoChemTestCase
{
    std::string name() const override
    {
        return "test_with_kinetics";
    }

    std::string getInputFileAsString() const override
    {
        return R"""(Temp 50.0
Pres 8.0e7
equilibrate 1

# Because the reference simulation used a different set of HKF parameters for the GCO3H^0 complex
SECONDARY_SPECIES
GCO3H = GCO3- - H2O + H+ / HKF  -92250  -98900  28.1  15.2964  -10.0382  -55.4493  7.4092  36.8069  3.5851  -0.3107 /
/end

geochem

solution 0
pH 7 charge
Ca 0.042
Cl 0.084
Mg 0.1
SO4 0.1
/end

# calcite 0.01 1 0 37.8 1.2e-8 8.4 0.12 0.5 1
rate
# rate = [(k_1*exp(-Ea/Rg)(1/T-1/298.15)+k_2*a_H*exp(-Ea/Rg)(1/T-1/298.15)]*(1-SI^n)^m
#mineral wt-fraction Sg log_af logEa_1 k_1 logEa_2 k_2 n   m
calcite 0.999 1 0 37.8 3.43E-03 8.4 1.11E+03 1 1
magnesite 0.0 1 0 60 .7E-06 0.0 0.1 1 1
anhydrite	0.0	1 0 60 .7E-08 0.0 0.0 1 1
/ end

equilibrium_phases
CO2(G) 1 -3.5
/end

/end)""";

    }
};

#endif
