/*
 * MIT License
 *
 * Copyright (C) 2025 Aksel Hiorth
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
*/
#ifndef GEOCHEM_PHASES_IS_DEF_H
#define GEOCHEM_PHASES_IS_DEF_H

#include <nlohmann/json.hpp>

#include <array>
#include <istream>
#include <iterator>
#include <map>
#include <string>
#include <tuple>
#include <set>
#include <vector>

enum class GeochemicalPhaseType
{
    AQUEOUS_SOLUTION = 0,
    SURFACE_COMPLEX = 1,
    EXCHANGE_SITES = 2,
    EQUILIBRIUM_MINERAL = 3,
    RATE_MINERAL = 4,
    GAS_PHASE = 5

};

/* Converts enum GeochemicalPhaseType to const char*, which can be further converted to std::string. */
constexpr const char* PhaseKeyword(GeochemicalPhaseType phase_type) noexcept
{
    switch(phase_type)
    {
    case GeochemicalPhaseType::AQUEOUS_SOLUTION: return "SOLUTION";
    case GeochemicalPhaseType::SURFACE_COMPLEX: return "COMPLEX";
    case GeochemicalPhaseType::EXCHANGE_SITES: return "IEXCHANGE";
    case GeochemicalPhaseType::EQUILIBRIUM_MINERAL: return "EQUILIBRIUM_PHASES";
    case GeochemicalPhaseType::GAS_PHASE: return "GAS_PHASE";
    case GeochemicalPhaseType::RATE_MINERAL: return "RATE";
    }
    return "";
}

// Define some useful aliases
using GeoChemBlockContent = std::map<std::string, std::string>;  //  parsed input data
using GeoChemPhaseData = std::map<std::string, GeoChemBlockContent>;
using PairOfPhaseNamesAndTypes = std::pair<std::vector<std::string>, std::vector<GeochemicalPhaseType>>;

/**
 * Class to keep track of all the aqueous solutions, minerals, ion exchangers,
 * etc. that can possibly exist during a simulation.
 *
 * Note that a keyword can be followed by a name, making the keywords unique,
 * e.g., "solution 0" and "solution 1" will be stored individually.
 */

// TODO: What if the same phase is entered twice, e.g., "SOLUTION X"...?
//       [If one keyword is present several times, only the last will be kept?]
//
class GeoChemPhases
{

  public:

    static constexpr std::array<GeochemicalPhaseType, 6> ALL_POSSIBLE_PHASES =
        {{
             GeochemicalPhaseType::AQUEOUS_SOLUTION,
             GeochemicalPhaseType::EXCHANGE_SITES,
             GeochemicalPhaseType::SURFACE_COMPLEX,
             GeochemicalPhaseType::EQUILIBRIUM_MINERAL,
             GeochemicalPhaseType::GAS_PHASE,
             GeochemicalPhaseType::RATE_MINERAL
         }};

    GeoChemPhases();
    GeoChemPhases(const GeoChemPhases&) = delete;
    GeoChemPhases(const GeoChemPhases&&) = delete;
    GeoChemPhases& operator=(const GeoChemPhases&) = delete;
    GeoChemPhases& operator=(const GeoChemPhases&&) = delete;

    void resetFromInputStream(std::istream& inputStream);

    /** accepts input of type: R"({"SOLUTION": {"SOLUTION 0": {"pH": "7 charge", "Ca": "1e-2", "Cl": "2e-2"}}})";*/
    // void resetFromJson(const std::string& json_input);
    nlohmann::json resetFromJson(const std::string& json_input);

    /** Updates an existing instance of GeoChemPhases based on the input.*/
    void updateFromInputStream(std::istream& inputStream);

    void clear();
    [[nodiscard]] bool isValid() const;

    /** @return True if there is at least one phase of the given type. */
    [[nodiscard]] bool hasPhase(GeochemicalPhaseType type) const;

    /** @return True if there is at least one element with name in the given type. */
    [[nodiscard]] bool hasElementInPhase(GeochemicalPhaseType type, const std::string& name) const;

    /** @return True if there is at least one element with name in any phase */
    [[nodiscard]] bool hasElementInAnyPhase(const std::string& name) const;

    /** @return True if there is a phase of the given type with the given name
     *          (case-sensitive).
     */
    [[nodiscard]] bool hasPhaseWithName(GeochemicalPhaseType type,
                                        const std::string& name) const;

    [[nodiscard]] int phaseSize(GeochemicalPhaseType type) const;

    /** @return: A vector with the names of all phases of the input type. */
    [[nodiscard]] std::vector<std::string> getAllNamesOfType(GeochemicalPhaseType type) const;

    /** Gain access to the underlying std::map which holds information about
     *  types, names, and values. */
    [[nodiscard]] const GeoChemPhaseData& getData(GeochemicalPhaseType type) const;

    /** @return: A pair of vectors holding 1) all the unique names and
    *            2) the corresponding types. */
    [[nodiscard]] PairOfPhaseNamesAndTypes getPhaseNamesAndTypes() const;

    void add_to_phase(GeochemicalPhaseType type,
                      const std::string& solution_name,
                      const std::string& specie_name,
                      const std::string& specie_values);

    // **********************************************************************
    //              HELPER FUNCTIONS USED TO INITIALIZE GEOCHEMICAL MODELS
    // **********************************************************************

    /**
     * Lazy routine:
     *
     *  Assumes that the first aqueous solution is for the "reservoir",
     *  and that the rest are to be injected.
     *  All surfaces, minerals, etc. are assumed to be part of reservoir.
     *
     * @return A tuple consisting of 1) the reservoir phases names & types,
     *        and 2) the names of the solutions to be injected.
     */
    [[nodiscard]] std::tuple<PairOfPhaseNamesAndTypes, std::vector<std::string>> GetInjectReservoirSimple() const;

    /** Similar to  GetInjectReservoirSimple(), but instead of choosing the
     *  first solution to be in the reservoir, look for it by name.
     */
    [[nodiscard]] std::tuple<PairOfPhaseNamesAndTypes, std::vector<std::string>> GetReservoirAndInjected(const std::string& reservoir_solution_name) const;

    /**
     * @return The maximum size to be passed to an InitChem "constructor".
     */
    int selectPhases(const PairOfPhaseNamesAndTypes& names_and_types,
                     std::map<GeochemicalPhaseType,
                         GeoChemPhaseData>& phasesToBeUsed) const;

  private:
    std::map<GeochemicalPhaseType, GeoChemPhaseData> phaseData_;

};

#endif
