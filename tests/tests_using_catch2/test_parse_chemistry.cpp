//  TODO:
//
//  There are (at least) two problems with how species charges and coefficients
//  are read:
//
//  1) If a specie name includes '-' or '+' in the name, the charge will be
//     incorrect (e.g., "K-FELDSPAR").
//  2) Species with negative coefficients are not read correctly when they
//     are read on their own. For example, when "-H+" is parsed as one term
//     in the equation for calcite, the code does the right thing.
//     However, the code does not work if we read "-H+" in isolation...
//     [Luckily, we have likely never encountered the error in practice]
//
#if __has_include(<catch2/catch.hpp>)
#include <catch2/catch.hpp>
#else
#include <catch2/catch_test_macros.hpp>
#endif

#include <opm/simulators/geochemistry/IO/ParseChemistry.hpp>

TEST_CASE("Parse neutral species")
{
    auto data = ParsedSpecieData("3Tra");
    REQUIRE(data.specie_name_ == "Tra");
    REQUIRE(data.specie_name_without_charge_ == "Tra");
    REQUIRE(data.coefficient_ == 3.0);
    REQUIRE(data.charge_ == 0.0);
}

TEST_CASE("Parse monovalent cations with negative coefficient")
{
    auto data = ParsedSpecieData("-4Na+");
    REQUIRE(data.specie_name_ == "Na+");
    REQUIRE(data.specie_name_without_charge_ == "Na");
    REQUIRE(data.coefficient_ == -4.0);
    REQUIRE(data.charge_ == 1.0);

    auto data_hplus = ParsedSpecieData("-H+");
    REQUIRE(data_hplus.specie_name_ == "H+");
    REQUIRE(data_hplus.specie_name_without_charge_ == "H");
    REQUIRE(data_hplus.coefficient_ == -1.0);
    REQUIRE(data_hplus.charge_ == 1.0);
    
}

TEST_CASE("Parse monovalent cation with positive coefficient")
{
    auto data = ParsedSpecieData("17K+");
    REQUIRE(data.specie_name_ == "K+");
    REQUIRE(data.specie_name_without_charge_ == "K");
    REQUIRE(data.coefficient_ == 17.0);
    REQUIRE(data.charge_ == 1.0);
}

TEST_CASE("Parse monovalent anion with positive coefficient")
{
    auto data = ParsedSpecieData("2Cl-");
    REQUIRE(data.specie_name_ == "Cl-");
    REQUIRE(data.specie_name_without_charge_ == "Cl");
    REQUIRE(data.coefficient_ == 2.0);
    REQUIRE(data.charge_ == -1.0);
}

TEST_CASE("Parse monovalent anion with negative coefficient")
{
    auto data = ParsedSpecieData("-6SCN-");
    REQUIRE(data.specie_name_ == "SCN-");
    REQUIRE(data.specie_name_without_charge_ == "SCN");
    REQUIRE(data.coefficient_ == -6.0);
    REQUIRE(data.charge_ == -1.0);
}

TEST_CASE("Parse divalent cation")
{
    auto data = ParsedSpecieData("2Ca+2");
    REQUIRE(data.specie_name_ == "Ca+2");
    REQUIRE(data.specie_name_without_charge_ == "Ca");
    REQUIRE(data.coefficient_ == 2.0);
    REQUIRE(data.charge_ == +2.0);
}

TEST_CASE("Parse divalent anion")
{
    auto data = ParsedSpecieData("3SO4-2");
    REQUIRE(data.specie_name_ == "SO4-2");
    REQUIRE(data.specie_name_without_charge_ == "SO4");
    REQUIRE(data.coefficient_ == 3.0);
    REQUIRE(data.charge_ == -2.0);
}

TEST_CASE("Parse an ion exchange secondary specie (analytical logK)")
{
    const auto&[data_vector, data_logK] = parseSecondarySpecie("CaX2 = 2X- + Ca+2 / ANA -0.8 /");
    
    REQUIRE(data_vector.size() == 3);
    REQUIRE(data_vector[0].specie_name_ == "CaX2");
    REQUIRE(data_vector[0].charge_ == 0.0);
    REQUIRE(data_vector[1].specie_name_ == "X-");
    REQUIRE(data_vector[1].charge_ == -1.0);
    REQUIRE(data_vector[2].specie_name_ == "Ca+2");
    REQUIRE(data_vector[2].charge_ == 2.0);

    REQUIRE(data_logK.model_ == LogKModel::ANA);
    REQUIRE(data_logK.mol_volume_ == 0.0); // Only used for minerals.

    // For analytical secondary species, we always add 3 zeroes first...
    REQUIRE(data_logK.values_.size() == 1+3);
    for(std::size_t i=0; i < 3; ++i)    REQUIRE(data_logK.values_[i] == 0.0);
    REQUIRE(data_logK.values_[3] == -0.8);
}


TEST_CASE("Parse HKF exchange secondary specie")
{
    const auto&[data_vector, data_logK] = parseSecondarySpecie("ACEH = ACE- + H+ /HKF -94760 -116100 42.7 11.6198 5.218 2.5088 -2.9946 42.076 -1.5417 -0.15 /");
    
    REQUIRE(data_vector.size() == 3);
    REQUIRE(data_vector[0].specie_name_ == "ACEH");
    REQUIRE(data_vector[0].charge_ == 0.0);
    REQUIRE(data_vector[0].coefficient_ == 1.0);
    REQUIRE(data_vector[1].specie_name_ == "ACE-");
    REQUIRE(data_vector[1].charge_ == -1.0);
    REQUIRE(data_vector[1].coefficient_ == 1.0);
    REQUIRE(data_vector[2].specie_name_ == "H+");
    REQUIRE(data_vector[2].charge_ == 1.0);
    REQUIRE(data_vector[2].coefficient_ == 1.0);

    REQUIRE(data_logK.model_ == LogKModel::HKF);
    REQUIRE(data_logK.mol_volume_ == 0.0); // Only used for minerals.

    REQUIRE(data_logK.values_.size() == 10);
    REQUIRE(data_logK.values_[0] == -94760);
    REQUIRE(data_logK.values_[1] == -116100);
    REQUIRE(data_logK.values_[2] == 42.7);
    REQUIRE(data_logK.values_[3] == 11.6198);
    REQUIRE(data_logK.values_[4] == 5.218);
    REQUIRE(data_logK.values_[5] == 2.5088);
    REQUIRE(data_logK.values_[6] == -2.9946);
    REQUIRE(data_logK.values_[7] == 42.076);
    REQUIRE(data_logK.values_[8] == -1.5417);
    REQUIRE(data_logK.values_[9] == -0.15);
}

TEST_CASE("Parse ANA mineral with a single, constant logK value")
{
    // Note: T=130 C.
    const auto&[data_vector, data_logK] = parseMineralPhase("bCALCITE = Ca+2 + HCO3- - H+ /ANA  36.918 0.373304");
    
    REQUIRE(data_vector.size() == 4);
    REQUIRE(data_vector[0].specie_name_ == "bCALCITE");
    REQUIRE(data_vector[0].charge_ == 0.0);
    REQUIRE(data_vector[0].coefficient_ == 1.0);
    REQUIRE(data_vector[1].specie_name_ == "Ca+2");
    REQUIRE(data_vector[1].charge_ == +2.0);
    REQUIRE(data_vector[1].coefficient_ == 1.0);
    REQUIRE(data_vector[2].specie_name_ == "HCO3-");
    REQUIRE(data_vector[2].charge_ == -1.0);
    REQUIRE(data_vector[2].coefficient_ == 1.0);
    REQUIRE(data_vector[3].specie_name_ == "H+");
    REQUIRE(data_vector[3].charge_ == 1.0);
    REQUIRE(data_vector[3].coefficient_ == -1.0);
    
    REQUIRE(data_logK.model_ == LogKModel::ANA);
    REQUIRE(data_logK.mol_volume_ == 36.918);
    
    // For analytical mineral phases, always add one zero first...
    REQUIRE(data_logK.values_.size() == 1 + 1);
    for(std::size_t i=0; i < 1; ++i)    REQUIRE(data_logK.values_[i] == 0.0);
    REQUIRE(data_logK.values_[1] == 0.373304);
}

TEST_CASE("Test parsing of ANA minerals with only one non-zero coefficient, but more than 7 coefficients in total")
{
    const auto&[data_vector, data_logK] = parseMineralPhase("MyMineral = Ca+2 + SO4-2 /ANA  42.0 0.314159 0 0 0 0 0 0 0 0 0");
    
    REQUIRE(data_vector.size() == 3);
    REQUIRE(data_vector[0].specie_name_ == "MyMineral");
    REQUIRE(data_vector[0].charge_ == 0.0);
    REQUIRE(data_vector[0].coefficient_ == 1.0);
    REQUIRE(data_vector[1].specie_name_ == "Ca+2");
    REQUIRE(data_vector[1].charge_ == +2.0);
    REQUIRE(data_vector[1].coefficient_ == 1.0);
    REQUIRE(data_vector[2].specie_name_ == "SO4-2");
    REQUIRE(data_vector[2].charge_ == -2.0);
    REQUIRE(data_vector[2].coefficient_ == 1.0);
    
    REQUIRE(data_logK.model_ == LogKModel::ANA);
    REQUIRE(data_logK.mol_volume_ == 42.0);
    
    // If more than 7 numbers are input (besides mol weight), they are removed.
    REQUIRE(data_logK.values_.size() == MAX_SIZE_ANA_MINERAL);  // The maximum size (!!)
    
    for(std::size_t i=0; i < MAX_SIZE_ANA_MINERAL; ++i)
    {
        if(i != 1)      REQUIRE(data_logK.values_[i] == 0.0);
    }
    REQUIRE(data_logK.values_[1] == 0.314159);
}

TEST_CASE("Attempt to parse ANA mineral with no logK value entered")
{
    const auto&[data_vector, data_logK] = parseMineralPhase("MyMineral = Ca+2 + SO4-2 /ANA  17.0");
    // Two zeroes are always added + one more since we input no logK value...
    REQUIRE(data_logK.values_.size() == 2 + 1);
    for(const auto& element: data_logK.values_) REQUIRE(element == 0.0);
}

TEST_CASE("Parse HKF mineral")
{
    const auto&[data_vector, data_logK] = parseMineralPhase("K-FELDSPAR = Al+3 + 2H2O - 4H+ + K+ + 3SiO2 /HKF 108.87 -895374 -949188 51.13 76.617 4.311 -29.945 0.0 -99 -99 -99");
    
    REQUIRE(data_vector.size() == 6);
    
    REQUIRE(data_vector[0].specie_name_ == "K-FELDSPAR");
    REQUIRE(data_vector[0].charge_ == -1.0);  // TODO: B/c of "-" in name...
    REQUIRE(data_vector[0].coefficient_ == 1.0);
    REQUIRE(data_vector[1].specie_name_ == "Al+3");
    REQUIRE(data_vector[1].charge_ == 3.0);
    REQUIRE(data_vector[1].coefficient_ == 1.0);
    REQUIRE(data_vector[2].specie_name_ == "H2O");
    REQUIRE(data_vector[2].charge_ == 0.0);
    REQUIRE(data_vector[2].coefficient_ == 2.0);
    REQUIRE(data_vector[3].specie_name_ == "H+");
    REQUIRE(data_vector[3].charge_ == 1.0);
    REQUIRE(data_vector[3].coefficient_ == -4.0);
    REQUIRE(data_vector[4].specie_name_ == "K+");
    REQUIRE(data_vector[4].charge_ == 1.0);
    REQUIRE(data_vector[4].coefficient_ == 1.0);
    REQUIRE(data_vector[5].specie_name_ == "SiO2");
    REQUIRE(data_vector[5].charge_ == 0.0);
    REQUIRE(data_vector[5].coefficient_ == 3.0);
    
    REQUIRE(data_logK.model_ == LogKModel::HKF);
    REQUIRE(data_logK.mol_volume_ == 108.87);

    REQUIRE(data_logK.values_.size() == 7);  // More than > 7 are removed (!!)
    REQUIRE(data_logK.values_[0] == -895374);
    REQUIRE(data_logK.values_[1] == -949188);
    REQUIRE(data_logK.values_[2] == 51.13);
    REQUIRE(data_logK.values_[3] == 76.617);
    REQUIRE(data_logK.values_[4] == 4.311);
    REQUIRE(data_logK.values_[5] == -29.945);
    REQUIRE(data_logK.values_[6] == 0.0);
}
