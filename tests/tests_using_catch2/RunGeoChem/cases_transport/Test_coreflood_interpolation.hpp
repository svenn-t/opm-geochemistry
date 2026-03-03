#ifndef TRANSPORT_COREFLOOD_INTERPOLATION_DEF
#define TRANSPORT_COREFLOOD_INTERPOLATION_DEF

#include "GeoChemTestCase.hpp"

struct Test_CorefloodInterpolation: public GeoChemTestCase
{
    std::string name() const override
    {
        return "Coreflood_interpolation";
    }

    GeoChemTestType type() const override
    {
        return GeoChemTestType::OneDTransport;
    }

    std::string getInputFileAsString() const override
    {
        return R"""(Temp 130
Pres 8.0e5

Tf  96
VolRate 0.03
Volume 33.3
Porosity  0.42
NoBlocks 50

WriteNetCDF 0
Imp 0
Flush 1.0

debug 0
Interpolate 1000

chemtol
1e-7 1e-8
/end

----------------------------------------------------------------
geochem
----------------------------------------------------------------

solution 0
--ekofisk_water
Ca		0.0000001
Cl		0.0000001
H   1.0
HCO3		0.0000001
K		0.0000001
Mg		0.0000001
Na		0.0000001
SO4		0.0000001
/ end

solution 1
-- sea_water
Ca		0.0130
Cl		0.1251
H   1.0
HCO3		0.002
K		0.01
Mg		0.0445
Na		0.05
SO4		0.024
/ end

solution 2
--torfisk_water (not currently used)
Ca		0.013
Cl		0.1251
H   1.0
HCO3		0.002
K		0.01
Mg		0.0445
Na		0.05
SO4		0.024
/ end

rate
# rate = (k_1*exp(-Ea/Rg)(1/T-1/298.15)+k_2*exp(-Ea/Rg)(1/T-1/298.15)*(1-SI^n)^m
#mineral   wt-fraction    Sg   log_af   logEa_1          k_1         logEa_2            k_2                n    m
calcite   		1.0   		1             0             40e+3     	3.43e-3               8.4e+3       1.11E+03             1             1
#calcite   		38.454   		1             0             40e+3     	3.43e-3               8.4e+3       1.11E+03             1             1
magnesite          0.0          1             0             60e+3   	0.00000004              0.0          0.0                  1             1
anhydrite          0.0          1             0             60e+3   	0.00000009                0.0          0.0                  1             1
#barite            0.0          2             0             30.8e+3     1.11E-03               0.0          1.11E-03             1              1
#talc              0.01         2             0             42          8.84E-09               0            0.00E+00             1             1
#strontianite      0.01        1             0             30.8         1.11E-10              30.8       0.00E+00                1             1
#muscovite         0.01        1             0             22           8.84E-15              30.8       0.00E+00                1             1
#kaolinite         0.01        1             0             22.2         1.47E-09               0             9.26E-08            1             1
#witherite         0.01        1             0             30.8         1.11E-04               48           1.11E-03             1             1
#paragonite        0.01        1             0             0             8.84E-17              23.6       8.84E+03               1             1
#quartz            0.01        1             0             0             2.80E-07               0         0.00E+00              1             1
/ end

iexchange
X 0.005
/ end

#complex
#method 1
#s_area 2000 10.0
#GCa 4.95
#GCO3 4.95
#/ end

/end)""";}
};

#endif
