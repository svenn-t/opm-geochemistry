#include "SetupAndRunSolver.hpp"

std::vector<EffluentIonData> gen_results_from_transport_calculation(const GeoChemTestCase& testcase,
                                                                    int serialize)
{
    if(testcase.type() != GeoChemTestType::OneDTransport) return {};

    auto string_input = testcase.getInputFileAsString();
    // Add "SERIALIZE 1" at the top to return a non-empty vector of results
    // Add "SERIALIZE 2" to also generate a binary file with sim. results
    if(serialize)
        string_input = "SERIALIZE " + std::to_string(serialize) + "\n" + string_input;
    std::stringstream input(string_input);

    OneDimensionalTransportSolver solver;
    return solver.solve(testcase.name(), input);
}
