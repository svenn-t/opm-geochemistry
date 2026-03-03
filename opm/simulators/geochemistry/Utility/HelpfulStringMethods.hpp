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
* A selection of functions used to manipulate std::strings.
*/
#ifndef HelpfulStringMethods_hpp
#define HelpfulStringMethods_hpp

#include <algorithm>
#include <sstream>
#include <string>
#include <vector>

template <typename T>
std::string convert_to_string(const T& value)
{
  std::ostringstream ss;
  ss << value;
  return ss.str();
}

/* Case-insensitive. Returns an iterator to the first element of the second group.*/
std::vector<std::string>::iterator partition_vector_by_key(std::vector<std::string>& unsorted_vector, const std::string& key);

std::vector<std::string> read_lines_from_stream(std::istream& inputStream, const std::string& endMarker);

int in_list(const std::string& name, const std::vector<std::string>& string_list);

std::pair<std::string, std::string> get_first_word(const std::string& line);

std::string remove_charge_from_species_name(const std::string& species_name);

/* Note: Finds the first occurrences of each delimiter. */
std::string substring_between(const std::string& str, const std::string& left_delimiter, const std::string& right_delimiter);

std::vector<std::string> split_by_delimiter(const std::string& string_to_split, char delimiter='\t');
std::vector<std::string> split_by_delimiter(const std::string& string_to_split, const std::string& delimiter="\t");
std::vector<std::string> split_string(const std::string& string_to_split);

bool string_contains(std::string_view str, std::string_view sub_str);

bool starts_with(const std::string& str, const std::string& start);

bool ends_with(const std::string & str, const std::string& ending);

bool has_content(const std::string& str);

std::string to_lower_case(const std::string& string_to_convert);
std::vector<std::string> to_lower_case(const std::vector<std::string>& strings_to_convert);

std::string to_upper_case(const std::string& string_to_convert);
std::vector<std::string> to_upper_case(const std::vector<std::string>& strings_to_convert);

std::string trim_left(const std::string& str, const std::string& t=" \t\n\r\f\v");
std::string trim_right(const std::string& str, const std::string& t=" \t\n\r\f\v");
std::string trim_left_right(const std::string& str, const std::string& t=" \t\n\r\f\v");




/* Suffixes are everything after ".". */
std::string remove_suffixes(const std::string& str);

/** Allows a mix of lower and upper case letters. */
int index_of_string_in_vector(const std::string& input_string, const std::vector<std::string>& vec);

/** Same as index_of_string_in_vector, except we do all comparisons in upper case. */
int index_of_string_in_vector_upper_case(const std::string& input_string, const std::vector<std::string>& vec);

#endif
