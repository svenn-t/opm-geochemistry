#include "SetupAndRunEqsolver.hpp"

BasVecInfo gen_results_from_equilibrium_calculation(GeoChemTestCase* testcase,
                                                    bool serialize)
{
    auto string_input = testcase->getInputFileAsString();
    // Possibly add keyword "SERIALIZE" at the top, to generate a binary file
    // with sim. results
    string_input = serialize ? "SERIALIZE 1\n" + string_input : string_input;
    std::stringstream input(string_input);

    EquilibriumSolver solver;
    return solver.solve(testcase->name(), input);
}
