/*
  Copyright 2025 Equinor ASA.

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/
#include "config.h"

#define BOOST_TEST_MODULE OpmGeochemInterfaceTests

#include <boost/test/unit_test.hpp>

#include <opm/simulators/geochemistry/OpmGeoChemInterface.hpp>

#include <algorithm>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>


BOOST_AUTO_TEST_CASE(InitializeFromDeckTest)
{
    // Setup
    std::string empty_file_name = "";
    bool charge_balance = true;
    const std::pair<double, double> tol = std::make_pair<double, double>(1e-5, 1e-6);
    int splay_tree_resolution = 10;

    // Initialize interface
    std::shared_ptr<OpmGeoChemInterface> geoChemInterface = std::make_shared<OpmGeoChemInterface>();
    const std::vector<std::string> species = { "H", "NA", "SO4", "K", "CL", "CA", "HCO3" };
    const std::vector<std::string> minerals = { "CAL" };
    const std::vector<std::string> ion_ex = { "X" };
    geoChemInterface->initialize_from_opm_deck(empty_file_name,
                                               species,
                                               minerals,
                                               ion_ex,
                                               charge_balance,
                                               tol,
                                               splay_tree_resolution);

    // Check number of species, minerals, and ion exchangers
    BOOST_CHECK_EQUAL(species.size(), geoChemInterface->numberOfAqueousBasisSpecies());
    BOOST_CHECK_EQUAL(minerals.size(), geoChemInterface->numberOfMinerals());
    BOOST_CHECK_EQUAL(ion_ex.size(), geoChemInterface->numberOfIonExchange());
    BOOST_CHECK_EQUAL(species.size() + ion_ex.size(), geoChemInterface->numberOfBasisSpecies());

    // Check if species names are located in interface
    const auto& basis_species = geoChemInterface->GetBasisSpeciesNames();
    for (const auto& elem : species) {
        BOOST_CHECK_MESSAGE(
            std::find(basis_species.begin(), basis_species.end(), elem) != basis_species.end(),
            "Species missing in GetBasisSpeciesNames() = " << elem
        );
    }

    // Check mineral name
    const auto& mineral_species = geoChemInterface->GetMineralNames();
    BOOST_CHECK_MESSAGE(
        std::find(mineral_species.begin(), mineral_species.end(), minerals[0]) != mineral_species.end(),
        "Mineral missing in GetMineralNames() = " << minerals[0]
    );

    // Check ion exchange name
    const auto& ionex_species = geoChemInterface->GetIONames();
    BOOST_CHECK_MESSAGE(
        std::find(basis_species.begin(), basis_species.end(), ion_ex[0]) != basis_species.end(),
        "Ion exchange missing in GetBasisSpeciesNames() = " << ion_ex[0]
    );
    BOOST_CHECK_MESSAGE(
        std::find(ionex_species.begin(), ionex_species.end(), ion_ex[0]) != ionex_species.end(),
        "Ion exchange missing in GetIONames() = " << ion_ex[0]
    );

    // Calculate initial mineral concentration test
    // Formula = total_rock_density / mole_weight_mineral * (1-poro) / poro * weight_fraction_mineral
    std::unordered_map<std::string, double> weight_mineral{ { minerals[0], 0.001 } };
    std::vector<double> Cmin;
    double porosity = 0.25;
    geoChemInterface->calculate_initial_mineral_concentration(Cmin, porosity, weight_mineral);
    double Cmin_expected = 0.08093;
    BOOST_CHECK_CLOSE(Cmin_expected, Cmin[0], 1e-2);

    // Calculate surface concentrations
    // NOTE: ion exchangers are inserted in Ctot _after_ aqueous species
    std::unordered_map<std::string, double> Cion { { ion_ex[0], 1.1e-3 } };
    double swat = 0.5;
    std::size_t nBasis = geoChemInterface->numberOfBasisSpecies();
    std::vector<double> Ctot(nBasis, 0.0);
    double frac_dl = 0.0;
    geoChemInterface->set_surface_concentrations(swat, Ctot, frac_dl, Cion);
    std::size_t nAq = nBasis - 1;
    BOOST_CHECK_CLOSE(Ctot[nAq], (1.0 / swat) * Cion.at(ion_ex[0]), 1e-6);

    // Check get functions
    for (const auto& elem : geoChemInterface->get_log_a_mineral()) {
        BOOST_CHECK_EQUAL(elem, 0.0);
    }
}
