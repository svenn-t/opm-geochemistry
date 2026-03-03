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
#ifndef PARSE_CHEMISTRY_DEF_H
#define PARSE_CHEMISTRY_DEF_H

#include <algorithm>
#include <cassert>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#include <opm/simulators/geochemistry/IO/ParseUtil.hpp>
#include <opm/simulators/geochemistry/IO/ErrorMsg.hpp>

static constexpr int MAX_SIZE_ANA_MINERAL = 8;  // deltaG, deltaH, a0, a1, a2, a3, a4, a5 [NOTE: S is not included]
static constexpr int MAX_SIZE_HKF_MINERAL = 7;  // deltaG, deltaH, S, a1, a2, a3, a7
static constexpr double DEFAULT_ION_SIZE = 4.0;

enum class LogKModel{ HKF=0, ANA=1 };

/**
* Helper struct used when parsing chemical species.
*/
struct ParsedSpecieData 
{
    ParsedSpecieData(const std::string& input_data);
    
    std::string specie_name_;
    std::string specie_name_without_charge_;
    double coefficient_;
    double charge_;
};

/** Helper struct for parsing either ANA or HKF parameters. */
struct LogKData
{
    LogKModel model_{LogKModel::HKF};
    double mol_volume_{ 0.0 };  // only used for minerals
    std::vector<double> values_;

};

std::tuple<double, double> parse_chemtol(const std::vector<std::string>& block_content);

std::tuple<std::vector<ParsedSpecieData>, LogKModel, std::vector<double>> parseReaction(const std::string& str_reaction);
std::tuple<std::vector<ParsedSpecieData>, LogKData> parseSecondarySpecie(const std::string& str_secondary_specie);
std::tuple<std::vector<ParsedSpecieData>, LogKData> parseMineralPhase(const std::string& str_mineral_phase);

#endif
