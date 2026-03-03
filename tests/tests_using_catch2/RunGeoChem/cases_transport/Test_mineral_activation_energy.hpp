#ifndef TRANSPORT_TEST_MINERAL_ACTIVATION_ENERGY_DEF
#define TRANSPORT_TEST_MINERAL_ACTIVATION_ENERGY_DEF

// We set calcite rate with a combination of k1 and Ea1.
// Results should agree with previous sim. using only k1.
// In addition, includes anhydrite and magnesite.
struct Test_MineralActivationEnergy: public GeoChemTestCase
{
    std::string name() const override
    {
        return "Test_mineral_activation_energy";
    }

    GeoChemTestType type() const override
    {
        return GeoChemTestType::OneDTransport;
    }

    std::string getInputFileAsString() const override
    {
        return R"""(Temp 130
Pres 8.0e5

Tf  25
VolRate 0.027
Volume 30.0
Porosity 0.5
NoBlocks 100

# Corresponds to a timestep of ~0.1 days: F=Q[mL/min]*dt[min]/Vp_block[mL].
Imp 1
Flush 6.5

WriteNetCDF 0
debug 0
Interpolate 0

chemtol
1e-7 1e-8
/end

----------------------------------------------------------------
geochem
----------------------------------------------------------------

solution 0
--distilled water
Ca		0.0000001
Cl		0.0000001
H               1.0
#pH 7 charge
HCO3	        0.0000001
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
#pH 7 charge
HCO3		0.002
K		0.01
Mg		0.0445
Na		0.05
SO4		0.024
/ end

rate
# rate = (k_1*exp(-Ea/Rg)(1/T-1/298.15)+k_2*exp(-Ea/Rg)(1/T-1/298.15)*(1-SI^n)^m
#mineral   wt-fraction    Sg   log_af   logEa_1      k_1         logEa_2        k_2                n             m
#
calcite          1        1       0      1.0e3      9.4528E-07      0.0          0.0               1             1
# Corresponds to:
#calcite          1        1       0      0.0        1.05E-06      .0            0.0               1             1
#
magnesite        0        1       0      0.0        2.01E-09       0.0           0.0               2             1
anhydrite        0        1       0      0.0        5.42E-07       0.0           0.0               1             1
/ end

debug 0
Interpolate 0

chemtol
1e-7 1e-8
/end

/end)""";
    }
};

# endif // TRANSPORT_TEST_MINERAL_ACTIVATION_ENERGY_DEF