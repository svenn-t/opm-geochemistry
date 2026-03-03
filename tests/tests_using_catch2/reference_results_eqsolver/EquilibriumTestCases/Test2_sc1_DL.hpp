#ifndef EQSOLVER_TEST2_SC1_DL_DEF
#define EQSOLVER_TEST2_SC1_DL_DEF

#include "GeoChemTestCase.hpp"

struct Test2_sc1_DL: public GeoChemTestCase
{
    std::string name() const override
    {
        return "test2_sc1_DL";
    }

    std::string getInputFileAsString() const override
    {
        return R"""(Temp 50.0
Pres 10.0e5
equilibrate 1

# Because the reference simulation used a different set of HKF parameters for the GCO3H^0 complex
SECONDARY_SPECIES
GCO3H = GCO3- - H2O + H+ / HKF  -92250  -98900  28.1  15.2964  -10.0382  -55.4493  7.4092  36.8069  3.5851  -0.3107 /
/end

geochem

solution 0
pH 7 charge
Na 0.01
K 0.01
Cl 0.02
SO4 0.005
Ca 0.005
HCO3 1.0e-8
/end

equilibrium_phases
calcite 1 0
CO2(G) 1 -3.5
/end

complex
method 1
#specific surface area m^2/L pore volume, size of diffusion layer
s_area 1000 10
#Name sites/nm^2
GCa 2
GCO3 2
/ end

iexchange
X 1.0e-2
/ end

/end)""";

    }
};

#endif
