#ifndef TRANSPORT_TEST_CALCITE_TO_MAGNESITE_DEF
#define TRANSPORT_TEST_CALCITE_TO_MAGNESITE_DEF

#include "GeoChemTestCase.hpp"

struct Test_Calcite_to_Magnesite: public GeoChemTestCase
{
    std::string name() const override
    {
        return "Test_Calcite_to_Magnesite";
    }

    GeoChemTestType type() const override
    {
        return GeoChemTestType::OneDTransport;
    }

    std::string getInputFileAsString() const override
    {
        return R"""(Temp 130
Pres 8.0e5

VolRate 0.02
Volume 33.3
Porosity  0.42

NoBlocks 5
Imp 0
Flush 1.0

MINERAL_PHASES
ncalcite = Ca+2 + HCO3- - H+ / HKF 36.934000    -269880.000000  -288552.000000  22.150000   24.980000   5.240000    -6.200000   1200.000000 /
/end

chemtol
1e-7 1e-8
/end

inject epoch1
1 10.0
/end

inject epoch2
2 20.0
/end

inject epoch3
3 5.0
/end

inject epoch4
2 20.0
/end

inject epoch5
2 10.0 0.08
/end

inject epoch6
2 10 0.02
/end

inject epoch7
2 10.0 0.08
/end

----------------------------------------------------------------
geochem
----------------------------------------------------------------

solution 0
-- distilled water
pH 7 charge
Na 1.0e-8
Cl 1.0e-8
/ end

solution 1
pH 7 charge
Cl 0.657
Na 0.657
/ end

solution 2
pH 7 charge
Cl 0.438
Mg 0.219
/ end

solution 3
-- distilled water
pH 7 charge
Na 1.0e-8
Cl 1.0e-8
/ end

rate
# rate = (k_1*exp(-Ea/Rg)(1/T-1/298.15)+k_2*a_H*exp(-Ea/Rg)(1/T-1/298.15))*(1-SI^n)^m
#mineral   wt-fraction	Sg   log_af   logEa_1	  k_1       logEa_2	k_2	        n	m
ncalcite 	1.0			1		0		37.8	3.43E-02	8.4		1.11E+03	1	1
magnesite	0.0			1		0		60		1.5E-08	    0.0		0.0			1	1
/ end

iexchange
X 0.02
/ end

/end)""";
    }
};

#endif
