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
/**
 * @file
 */
#ifndef CHEMGLOBAL_H
#define CHEMGLOBAL_H

#include <cmath>
#include <cstdlib>
#include <numeric>
#include <type_traits>
#include <string_view>
#include <vector>

#include <opm/simulators/geochemistry/Common/Constants.hpp>

using RealVec = std::vector<double>;

class ChemTable;
void calculate_mol_weight_mineral(ChemTable& Min, ChemTable& Bas);

static inline double POW10(double x)
{
    return std::exp(NumericalConstants::LNTEN*x);
}

bool is_gas_buffer(std::string_view buffer_name);

int number_of_combinations(int size);

void convert_into_binary(int number_in_base10, int* storage_array, int array_size);

/**
 * Creates a unique key for each combination of positive array entries.
 * Note that the actual values does not matter, as long as they are positive.
 * Specifically, the effect of the function consist of two steps:
 *  1) Convert the input array into a binary number by replacing all positive
 *     numbers with 1 and all other numbers by 0.
 *  2) Return the resulting binary number in base 10.
 *
 * Examples:
 *  [1, 2, 0]           ->     (110)_2     ==     6
 *  [42, 99, 0]         ->     (110)_2     ==     6
 *  [0.003, 0, 1]       ->     (101)_2     ==     5
 *  [0, 1e-4, 0.1]      ->     (011)_2     ==     3
 */
template<class T, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
int create_unique_key_from_positive_elements(const T* array_of_T, int array_size)
 {
     int key = 0;
     for(int i=0; i < array_size; ++i){
         key = 2*key + (array_of_T[i] > 0 ? 1 : 0 );
     }
     return key;
 }

int negative_floor_function(double x);

#endif
