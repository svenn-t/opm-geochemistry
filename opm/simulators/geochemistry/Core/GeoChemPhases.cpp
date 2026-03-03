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
#include <opm/simulators/geochemistry/Core/GeoChemPhases.hpp>
#include <opm/simulators/geochemistry/Common/CustomExceptions.hpp>
#include <opm/simulators/geochemistry/IO/InputReader.hpp>
#include <opm/simulators/geochemistry/Utility/HelpfulStringMethods.hpp>
#include <opm/simulators/geochemistry/Utility/MiscUtilityFunctions.hpp>

#include <fstream>
#include <cassert>

GeoChemPhases::GeoChemPhases()
    : phaseData_()
{
    for(const auto& phaseType: GeoChemPhases::ALL_POSSIBLE_PHASES)
    {
        phaseData_[phaseType] = GeoChemPhaseData();
    }
}

void GeoChemPhases::resetFromInputStream(std::istream& inputStream)
{
    clear();
    updateFromInputStream(inputStream);
    assert(isValid());
    if (!isValid()) throw MissingKeywordException("Geochemical model could not be initialized, there must be at least one solution!");
}

// void GeoChemPhases::resetFromJson(const std::string& json_input)
nlohmann::json GeoChemPhases::resetFromJson(const std::string& json_input)
{
    nlohmann::json json_parsed;
    try {
        // Check if json_input is a file path
        std::ifstream file(json_input);
        if (file.good()) {
            std::ifstream in(json_input);
            if (!in) {
                throw std::runtime_error("Could not open JSON file: " + json_input);
            }
            in >> json_parsed;
        } else {
            // Fallback: parse as raw JSON string
            json_parsed = nlohmann::json::parse(json_input);
        }

        for (const auto& phaseType : GeoChemPhases::ALL_POSSIBLE_PHASES)
        {
            const std::string key = PhaseKeyword(phaseType);
            if (json_parsed.contains(key))
            {
                const auto& section = json_parsed[key];
                for (const auto& [phase_name, values] : section.items())
                {
                    for (const auto& [element_name, element_val] : values.items())
                    {
                        add_to_phase(phaseType, phase_name, element_name, element_val);
                    }
                }
            }
        }

    } catch (const nlohmann::json::parse_error& e) {
        std::cerr << "JSON parse error in GeoChemPhases: " << e.what() << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Error in GeoChemPhases::resetFromJson: " << e.what() << std::endl;
    }

    return json_parsed;
}
void GeoChemPhases::updateFromInputStream(std::istream& inputStream)
{
    InputReader chemFileReader({}, /*case_sensitive=*/ false);
    std::vector<std::string> end_keywords = { R"(/end)", "end", R"(/ end)" };
    chemFileReader.define_main_block("geochem", end_keywords);

    for(const auto& phaseType: GeoChemPhases::ALL_POSSIBLE_PHASES)
    {
        chemFileReader.define_block_keyword(PhaseKeyword(phaseType), end_keywords);
    }

    chemFileReader.read(inputStream);

    for (const auto& phaseType : GeoChemPhases::ALL_POSSIBLE_PHASES)
    {
        const auto phaseName = static_cast<std::string>(PhaseKeyword(phaseType));

        for (const auto& uniqueName : chemFileReader.unique_names(phaseName))
        {
            const auto& content_as_vector = chemFileReader.get_block_keyword_content(uniqueName);

            for (const auto& content : content_as_vector)
            {
                std::pair<std::string, std::string> key = get_first_word(content);

                add_to_phase(phaseType, uniqueName, key.first, key.second);
            }
        }
    }
}

void GeoChemPhases::clear() {
    phaseData_.clear();
}

bool GeoChemPhases::isValid() const
{
    // At minimum, we require there to be one aqueous solution!
    if (phaseSize(GeochemicalPhaseType::AQUEOUS_SOLUTION) == 0) return false;
    return true;
}

bool GeoChemPhases::hasPhase(GeochemicalPhaseType type) const
{
    return phaseData_.count(type) > 0;  // Note: count can only return 0 or 1
}

bool GeoChemPhases::hasElementInPhase(GeochemicalPhaseType type, const std::string& name) const
{
    if (hasPhase(type))
    {
        for (const auto&  [elem, phase] : phaseData_.at(type))
        {
            for (const auto&  [name_i, val_i] : phase)
                if (name_i == name)
                    return true;
        }
    }

    return false;
}

bool GeoChemPhases::hasElementInAnyPhase(const std::string& name) const
{
    for(const auto& phaseType: GeoChemPhases::ALL_POSSIBLE_PHASES)
    {
        if (hasElementInPhase(phaseType, name))
            return true;
    }
    return false;
}

bool GeoChemPhases::hasPhaseWithName(GeochemicalPhaseType type, const std::string& name) const
{
    const auto all_names_of_type = getAllNamesOfType(type);
    return element_is_in_container(name, all_names_of_type);
}

int GeoChemPhases::phaseSize(GeochemicalPhaseType type) const
{
    if(phaseData_.find(type) == phaseData_.end()){
        return 0;
    }
    return static_cast<int>(phaseData_.at(type).size());
}

std::vector<std::string> GeoChemPhases::getAllNamesOfType(GeochemicalPhaseType type) const
{
    if(phaseData_.find(type) == phaseData_.end())
    {
        return {};  // Not found...
    }

    const auto& mapForType = phaseData_.at(type);

    std::vector<std::string> vec;
    for (const auto& [name, value_not_used] : mapForType)
    {
        vec.push_back(name);
    }
    return vec;
}

/* Assumes phaseSize(type) >= 1. */
const GeoChemPhaseData& GeoChemPhases::getData(GeochemicalPhaseType type) const
{
    return phaseData_.at(type);
}

PairOfPhaseNamesAndTypes GeoChemPhases::getPhaseNamesAndTypes() const
{
    std::vector<std::string> names;
    std::vector<GeochemicalPhaseType> types;

    for (const auto& phaseType : GeoChemPhases::ALL_POSSIBLE_PHASES)
    {
        if (phaseData_.find(phaseType) != phaseData_.end())
        {
            for (const auto& [name, not_used] : phaseData_.at(phaseType))
            {
                names.push_back(name);
                types.push_back(phaseType);
            }
        }
    }
    return std::make_pair(names, types);
}

void GeoChemPhases::add_to_phase(GeochemicalPhaseType type,
                                 const std::string& solution_name,
                                 const std::string& specie_name,
                                 const std::string& specie_values)
{
    phaseData_[type][solution_name].insert_or_assign(specie_name, specie_values);
}

std::tuple<PairOfPhaseNamesAndTypes, std::vector<std::string>> GeoChemPhases::GetInjectReservoirSimple() const
{
    PairOfPhaseNamesAndTypes reservoir_phases;
    std::vector<std::string> names_of_injected_solutions;

    const auto solution_phases = getPhaseNamesAndTypes();

    bool read_first_solution = false;

    for (std::size_t i = 0; i < solution_phases.first.size(); ++i)
    {
        if (solution_phases.second[i] == GeochemicalPhaseType::AQUEOUS_SOLUTION && !read_first_solution)
        {
            reservoir_phases.first.push_back(solution_phases.first[i]);
            reservoir_phases.second.push_back(solution_phases.second[i]);
            read_first_solution = true;
        }
        else if (solution_phases.second[i] == GeochemicalPhaseType::AQUEOUS_SOLUTION)
        {
            names_of_injected_solutions.push_back(solution_phases.first[i]);
        }
        else // the rest of the phases are for the reservoir
        {
            reservoir_phases.first.push_back(solution_phases.first[i]);
            reservoir_phases.second.push_back(solution_phases.second[i]);
        }
    }

    return std::make_tuple(reservoir_phases, names_of_injected_solutions);
}

std::tuple<PairOfPhaseNamesAndTypes, std::vector<std::string>> GeoChemPhases::GetReservoirAndInjected(const std::string& reservoir_solution_name) const
{
    PairOfPhaseNamesAndTypes reservoir_phases;
    std::vector<std::string> names_of_injected_solutions;

    const auto solution_phases = getPhaseNamesAndTypes();

    bool found_reservoir_solution = false;
    for (std::size_t i = 0; i < getPhaseNamesAndTypes().first.size(); ++i)
    {
        const auto& name = solution_phases.first[i];
        const auto& type = solution_phases.second[i];

        if (!found_reservoir_solution
            && type == GeochemicalPhaseType::AQUEOUS_SOLUTION
            && name == reservoir_solution_name)
        {
            reservoir_phases.first.push_back(name);
            reservoir_phases.second.push_back(type);
            found_reservoir_solution = true;
        }
        else if (type == GeochemicalPhaseType::AQUEOUS_SOLUTION)
        {
            names_of_injected_solutions.push_back(name);
        }
        else // the rest of the phases are for the reservoir
        {
            reservoir_phases.first.push_back(name);
            reservoir_phases.second.push_back(type);
        }
    }

    return std::make_tuple(reservoir_phases, names_of_injected_solutions);
}

int GeoChemPhases::selectPhases(
    const PairOfPhaseNamesAndTypes& names_and_types,
    std::map<GeochemicalPhaseType,
        GeoChemPhaseData>& phasesToBeUsed
) const
{
    // TODO: This function could be improved...(or removed?)
    //
    phasesToBeUsed.clear();

    const auto& phases_names = names_and_types.first;
    const auto& phases_types = names_and_types.second;

    std::vector<std::string> added_keys;

    // Go through the species in the input phases and add those that exist
    // [inside this GeoChemPhases object]
    for (std::size_t i = 0; i < phases_types.size(); ++i)
    {
        const auto& phase_type = phases_types[i];
        const auto& phase_name = phases_names[i];

        if (hasPhaseWithName(phase_type, phase_name))
        {
            phasesToBeUsed[phase_type][phase_name] = getData(phase_type).at(phase_name);

            const auto& data = getData(phase_type).at(phase_name);
            for (const auto& [key, content_not_used] : data)
            {
                added_keys.push_back(key);
            }
        }
    }

    // Note: Since some keys are not basis species, max_size will generally
    //       be an overestimate.
    std::set<std::string> unique_keys(added_keys.begin(), added_keys.end());
    return static_cast<int>(unique_keys.size());  // max_size
}
