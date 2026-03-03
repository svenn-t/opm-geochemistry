#ifndef EQSOLVER_TEST1_REDOX
#define EQSOLVER_TEST1_REDOX

#include "GeoChemTestCase.hpp"

struct Test1_Redox: public GeoChemTestCase
{
    std::string name() const override
    {
        return "test1_redox";
    }

    std::string getInputFileAsString() const override
    {
        return R"""(Temp 25.0
Pres 1.010325e5
equilibrate 1

geochem

solution 0
pH 7 charge
e- 1e-4
Ca 1e-3
HCO3 1e-3
Na 0.1
Cl 0.2
SO4 0.05
/end

equilibrium_phases
calcite 1 0.0
CO2(G) 1 -3.5
N2(G) 1 -3.0
/end

/end)""";

    }
};

#endif
