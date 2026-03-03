#ifndef EQSOLVER_TEST1_IO_ONLY_DEF
#define EQSOLVER_TEST1_IO_ONLY_DEF

#include "GeoChemTestCase.hpp"

struct Test1_IO_only: public GeoChemTestCase
{
    std::string name() const override
    {
        return "test1_IO_only";
    }

    std::string getInputFileAsString() const override
    {
        return R"""(Temp 130.0
Pres 1.0e7
equilibrate 1

# Because the reference simulation used a different set of HKF parameters for the GCO3H^0 complex
SECONDARY_SPECIES
GCO3H = GCO3- - H2O + H+ / HKF  -92250  -98900  28.1  15.2964  -10.0382  -55.4493  7.4092  36.8069  3.5851  -0.3107 /
/end

geochem

solution 0
pH 7
Na 48e-3
SO4 24e-3
K 24e-3
Cl 24e-3
/end

iexchange
X 1.0e-3
/ end

/end)""";

    }
};

#endif
