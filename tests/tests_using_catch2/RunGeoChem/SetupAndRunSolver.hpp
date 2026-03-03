#ifndef SETUP_AND_RUN_GEOCHEM_SOLVER_DEF
#define SETUP_AND_RUN_GEOCHEM_SOLVER_DEF

#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "GeoChemTestCase.hpp"

#include <opm/simulators/geochemistry/StandaloneSolvers.hpp>

std::vector<EffluentIonData> gen_results_from_transport_calculation(const GeoChemTestCase& testcase,
                                                                    int serialize=0);

#endif
