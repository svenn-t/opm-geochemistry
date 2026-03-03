#ifndef TRANSPORT_TEST_RATE_AND_EQUILIBRIUM_MINERALS_DEF
#define TRANSPORT_TEST_RATE_AND_EQUILIBRIUM_MINERALS_DEF

#include "GeoChemTestCase.hpp"

struct Test_rate_and_equilibrium_minerals: public GeoChemTestCase
{
    std::string name() const override
    {
        return "test_rate_and_equilibrium_minerals";
    }

    GeoChemTestType type() const override
    {
        return GeoChemTestType::OneDTransport;
    }

    std::string getInputFileAsString() const override
    {
        return R"""(Temp 130
Pres 800000.0

Tf  10

VolRate 0.2
Volume 38.31
Porosity 0.4174

NoBlocks 10

Imp 0
Flush 0.5

WriteBlock 0
WriteNetCDF 0
Interpolate 0

NumericalJacobian 0
debug 0

BASIS_SPECIES
# Name a0 Mw
SCN- 5 58.08 / HKF 909090 9 /
/end

SECONDARY_SPECIES
# DeltaG	DeltaH	S	a1	a2	a3	a4	c1	c2	omega
NaX = X- + Na+ / ANA 0. /
CaX2 = 2X- + Ca+2 / ANA -0.8 /
MgX2 = 2X- + Mg+2 / ANA -0.6 /
KX = X- + K+ / ANA -0.7 /
BaX2 = 2X- + Ba+2 / ANA -0.91 /
SrX2 = 2X- + Sr+2 / ANA -0.91 /
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

complex
method 1
s_area 3000.0 9.096887879206474
GCa 4.95
GCO3 4.95
/ end

iexchange
X 0.05
/ end

rate
# rate = (k_1*exp(-Ea/Rg)(1/T-1/298.15)+k_2*exp(-Ea/Rg)(1/T-1/298.15)*(1-SI^n)^m
#mineral   wt-fraction	Sg   log_af   logEa_1	  k_1         logEa_2	k_2	        n	m
#calcite 	38.454	1	0	37.8	3.43E-03	8.4	1.11E+01	1	1
calcite 	1.0	1	0	37.8	3.43E-03	8.4	1.11E+01	1	1
magnesite	   0	1	0	60	.7E-09	     	0.0	0.0		1	1
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
Ca 0.05
Cl 0.1
/end

inject Epoch1
1 5
/end

inject Epoch2
2 5
/end

/end)""";
    }
};

#endif
