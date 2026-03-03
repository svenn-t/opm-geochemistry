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
*
* A collection of utility functions for file-related stuff.
*/
#ifndef FILE_FUNCTIONS_DEF_H
#define FILE_FUNCTIONS_DEF_H

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

FILE* my_fopen(const char* filename, const char* mode);

std::string remove_quotation_marks_from_filename(const std::string& str);

/**
 * @param [in,out] out: Output stream to write to.
 * @param [in] values: Vector of values to write to the stream.
 * @param [in] no_values: Number of values to print.
 * @param [in] delim: Character to write in-between subsequent values.
 */
template<typename T>
void write_output_to_stream(std::ostream& out, const std::vector<T>& values, std::size_t no_values, char delim='\t')
{
    std::size_t no_written_values = 0;

    for(const auto& val: values)
    {
        if(no_written_values == 0) out << val;
        else out << delim << val;

        ++no_written_values;
        if(no_written_values == no_values) break;
    }
    out << "\n";
}

/**
 * Overloaded function that prints all of the values in the input vector.
 *
 * @param [in,out] out: Output stream to write to.
 * @param [in] values: Vector of values to write to the stream.
 * @param [in] delim: Character to write in-between subsequent values.
 */
template<typename T>
void write_output_to_stream(std::ostream& out, const std::vector<T>& values, char delim='\t')
{
    write_output_to_stream(out, values, values.size(), delim);
}

#endif
