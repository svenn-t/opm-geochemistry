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
#include <opm/simulators/geochemistry/Common/ChemGlobal.h>

#include <opm/simulators/geochemistry/Core/ChemTable.h>

void calculate_mol_weight_mineral(ChemTable& Min,  ChemTable& Bas)
{
    Min.mol_weight_.resize(Min.noRows_);

    for(int i=0; i < Min.noRows_; ++i)
    {
        Min.mol_weight_[i] = 0.0;
        for(int j=0; j < Min.noColumns_; ++j){
            Min.mol_weight_[i] += Min.M_[i][j]*Bas.mol_weight_[j];
        }
    }
}
bool is_gas_buffer(std::string_view buffer_name)
{
    if(buffer_name.find(",g") != std::string::npos)     return true;
    if(buffer_name.find(",G") != std::string::npos)     return true;
    if(buffer_name.find("(g)") != std::string::npos)    return true;
    if(buffer_name.find("(G)") != std::string::npos)    return true;
    return false;
}

int number_of_combinations(int size)
{
    return static_cast<int>(pow(2.0, size));
}

/**
 * The binary number is stored in the input dynamic array.
 * Thus, it is assumed that the input number is >= 0 and < 2^array_size.
 */
void convert_into_binary(int number_in_base10, int* storage_array, int array_size)
{
    int quotient = number_in_base10;
    for (int i = 0; i < array_size; ++i) {
        auto [quot, rem] = std::div(quotient, 2);
        storage_array[array_size - i - 1] = rem;
        quotient = quot;
    }
}

/**
 * Note that this function rounds differently than the function std::floor in
 * the standard library. Example: If the input is -3, the return value is -4.
 */
int negative_floor_function(double x)
{
    if(x < 0)       return static_cast<int>(x) - 1;
    else            return static_cast<int>(x);
}
