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
#include <opm/simulators/geochemistry/IO/ParseChemistry.hpp>

#include <opm/simulators/geochemistry/Common/CustomExceptions.hpp>
#include <opm/simulators/geochemistry/Utility/HelpfulStringMethods.hpp>

/**
* Construct ParsedSpeciesData instance from a string.
* Example: "2Ca+2" --> species_name="Ca+2", without_charge="Ca", coefficient=2.0, charge=+2.
*/
ParsedSpecieData::ParsedSpecieData(const std::string& input_data)
{
    auto data = trim_left(input_data);

    coefficient_ = std::atof(data.c_str());
    if(coefficient_ == 0.0 && starts_with(data, "-"))
    {
        coefficient_ = -1.0;
    }
    else if (coefficient_ == 0.0)
    {
        coefficient_ = 1.0;
    }

    const std::size_t offset = (coefficient_ < 0.0) ? 1 : 0;
    const std::size_t pos = data.find_first_of("+-", offset);

    charge_ = 0.0;
    if (pos != std::string::npos)
    {
        charge_ = std::atof(data.substr(pos, std::string::npos).c_str());
        if (charge_ == 0.0)
        {
            if (data[pos] == '-')   charge_ = -1.0;
            else                    charge_ = 1.0;
        }
    }

    std::string specie_name_without_coefficient;
    for (char c : data)
    {
        // Is the second condition always true? ('9' comes after '.' in ascii?)
        if ((c > '9' && c != '.') || !specie_name_without_coefficient.empty())
            specie_name_without_coefficient.push_back(c);
    }

    specie_name_ = specie_name_without_coefficient;
    specie_name_without_charge_ = specie_name_.substr(0, specie_name_.find_first_of("+-"));
}

/** If the input can not be read, returns sensible defaults.*/
std::tuple<double, double> parse_chemtol(const std::vector<std::string>& block_content)
{
    double mbal_crit = 1.0e-5;
    double pH_crit = 1.0e-5;
    
    if(!block_content.empty())
    {
        const std::vector<std::string> lineTab = split_string(block_content[0]);
        if (!lineTab.empty()) mbal_crit = std::stod(lineTab[0]);
        if (lineTab.size() > 1) pH_crit = std::stod(lineTab[1]);
    }
    return std::make_tuple(mbal_crit, pH_crit);
}

std::tuple<std::vector<ParsedSpecieData>, LogKModel, std::vector<double>> parseReaction(const std::string& str_reaction)
{
    std::vector<ParsedSpecieData> vector_of_species_data;

    // Parse left-hand side
    std::stringstream iss(substring_between(str_reaction, "", "/"));
    std::string tmp_data;
    iss >> tmp_data;
    vector_of_species_data.emplace_back(tmp_data);

    // Next, read the dependent (right-hand side) species
    char p;
    while (iss >> p >> tmp_data)
    {
        vector_of_species_data.emplace_back(tmp_data);
        if (p == '-') vector_of_species_data.back().coefficient_ *= -1.0;
        else assert(p == '+' || p == '=');
    }

    // Check if "=" is part of first specie name 
    if (vector_of_species_data[0].specie_name_.find("=") != std::string::npos)
    {
        error_and_exit("Possible error detected \"=\" should not be part of specie name {:s}, include a space in the reaction.", vector_of_species_data[0].specie_name_);
    }
    // Finally, read model type and values (which will be inspected later)
    const std::size_t indexHKF = str_reaction.find("HKF");
    const std::size_t indexANA = str_reaction.find("ANA");
    if (indexHKF != std::string::npos)
    {
        assert(indexANA == std::string::npos); // Cannot have both.
        auto values = parse_doubles(substring_between(str_reaction, "HKF", "/"));
        return std::make_tuple(std::move(vector_of_species_data), LogKModel::HKF, std::move(values));
    }
    else if (indexANA != std::string::npos)
    {
        auto values = parse_doubles(substring_between(str_reaction, "ANA", "/"));
        return std::make_tuple(std::move(vector_of_species_data), LogKModel::ANA, std::move(values));
    }
    else
    {
        throw ModelDoesNotExistException("Secondary specie must be specified by either ANA or HKF keyword.");
    }
}

/** Parse a single line of the SECONDARY_SPECIES keyword. */
std::tuple<std::vector<ParsedSpecieData>, LogKData> parseSecondarySpecie(const std::string& str_secondary_specie)
{
    auto [vector_of_species_data, model, values] = parseReaction(str_secondary_specie);

    LogKData data;
    data.model_ = model;
    data.mol_volume_ = 0.0;
    if (model == LogKModel::ANA)
    {
        data.values_.resize(values.size() + 3);
        for(std::size_t i=0; i < 3; ++i)    data.values_[i] = 0.0;
        std::copy_n(std::make_move_iterator(begin(values)), values.size(), data.values_.begin()+3);
    }
    else data.values_ = std::move(values);

    return std::make_tuple(std::move(vector_of_species_data), std::move(data));
}

/** Parse a single line of the MINERAL_PHASES keyword. */
std::tuple<std::vector<ParsedSpecieData>, LogKData> parseMineralPhase(const std::string& str_mineral_phase)
{
    auto [vector_of_species_data, model, values] = parseReaction(str_mineral_phase);

    LogKData data;
    data.model_ = model;

    // Remove molar volume from vector
    if (values.size() > 0)
    {
        if (values[0] < 0)
        {
            error_and_exit("Molar volume (first entry) should not be negative for reaction and data {:s}.", str_mineral_phase);
        }
        data.mol_volume_ = values[0];
        values.erase(values.begin());
    }
    else
    {
        error_and_exit("No molar volume entered for reaction and data {:s}.", str_mineral_phase);
    }

    int max_size = (model == LogKModel::ANA) ? MAX_SIZE_ANA_MINERAL : MAX_SIZE_HKF_MINERAL;

    if (model == LogKModel::ANA && values.empty())
    {
        data.values_ = std::vector<double>(3, 0.0);
    }
    else if(model == LogKModel::ANA)
    {
        data.values_.resize(values.size()+1);
        for(std::size_t i=0; i < 1; ++i)    data.values_[i] = 0.0;
        std::copy_n(std::make_move_iterator(begin(values)), values.size(), data.values_.begin()+1);
    }
    else data.values_ = std::move(values);

    if(data.values_.size() > max_size) data.values_.resize(max_size);

    return std::make_tuple(std::move(vector_of_species_data), std::move(data));
}