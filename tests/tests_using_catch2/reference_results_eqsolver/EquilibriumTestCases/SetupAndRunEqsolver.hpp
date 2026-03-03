#ifndef SETUP_AND_RUN_EQSOLVER_DEF
#define SETUP_AND_RUN_EQSOLVER_DEF

#include <sstream>
#include <string>

#include "GeoChemTestCase.hpp"

#include <opm/simulators/geochemistry/StandaloneSolvers.hpp>

BasVecInfo gen_results_from_equilibrium_calculation(GeoChemTestCase* testcase,
                                                    bool serialize=false);

#endif
