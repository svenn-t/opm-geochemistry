#ifndef TRANSPORT_INTERPOLATION_MINERALS_BUT_NO_SURF_DEF
#define TRANSPORT_INTERPOLATION_MINERALS_BUT_NO_SURF_DEF

#include "GeoChemTestCase.hpp"

struct Test_InterpolationWithMineralsButNoSurfSpecies: public GeoChemTestCase
{
    std::string name() const override
    {
        return "InterpolationWithMineralsButNoSurfSpecies";
    }

    GeoChemTestType type() const override
    {
        return GeoChemTestType::OneDTransport;
    }

    std::string getInputFileAsString() const override
    {
        return R"""(Temp 80
Pres 800000.0

Tf  10

VolRate 0.2
Volume 40
Porosity 0.4

NoBlocks 9

Imp 0
Flush 0.75
Interpolate 100

BASIS_SPECIES
# Name a0 Mw
SCN- 5 58.08 / HKF 909090 9 /
/end

chemtol
1e-7 1e-8
/end

geochem

--Initial solution in column
solution 0
pH 7 charge
HCO3	0.002
Cl	0.623
Mg	0.0445
Ca	0.013
Na	0.5
K	0.01
/end

rate
# rate = (k_1*exp(-Ea/Rg)(1/T-1/298.15)+k_2*exp(-Ea/Rg)(1/T-1/298.15)*(1-SI^n)^m
#mineral   wt-fraction	Sg   log_af   logEa_1	  k_1         logEa_2	k_2	        n	m
calcite 	1.0	1	0	37.8	3.43E-03	8.4	1.11E+01	1	1
dolomite	0	1	0	60	.7E-09	     	0.0	0.0		1	1
anhydrite	0.0	1	0	60	.3E-07		0.0	0.0		1	1
/ end

solution 1
pH 7 charge
HCO3	0.002
Cl	0.525
SO4	0.024
SCN-	0.024
Mg	0.0445
Ca	0.013
Na	0.45
K	0.034
/end

solution 2
pH 7 charge
-- sea_water
Ca		0.0130
Cl		0.1251
H   1.0
HCO3    0.002
K		0.01
Mg		0.0445
Na		0.05
SO4		0.024
/end

inject Epoch1
1 5
/end

inject Epoch2
2 5
/end

inject Epoch3
1 5
/end

/end

)""";}
};

#endif
