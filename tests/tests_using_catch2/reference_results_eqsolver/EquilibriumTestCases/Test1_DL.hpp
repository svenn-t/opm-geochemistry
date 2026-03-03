#ifndef EQSOLVER_TEST1_DL_DEF
#define EQSOLVER_TEST1_DL_DEF

#include "GeoChemTestCase.hpp"

struct Test1_DL: public GeoChemTestCase
{
    std::string name() const override
    {
        return "test1_DL";
    }

    std::string getInputFileAsString() const override
    {
        return R"""(Temp 130.0
Pres 1.0e7
equilibrate 1

# Because the reference simulation used a different set of HKF parameters for the GCO3H^0 complex
SECONDARY_SPECIES
GCO3H = GCO3- + H+ - H2O / HKF  -92250  -98900  28.1  15.2964  -10.0382  -55.4493  7.4092  36.8069  3.5851  -0.3107 /
/end

geochem

solution 0
pH 7
Na 48e-3
SO4 24e-3
K 24e-3
Cl 24e-3
Ca 1.0e-8
HCO3 1.0e-8
/end

complex
method 1
#specific surface area m^2/L pore volume, size of diffusion layer
s_area 1000 10
#Name sites/nm^2
GCa 2
GCO3 2
/ end

equilibrium_phases
calcite 1 0
/end

/end)""";

    }
};

#endif
