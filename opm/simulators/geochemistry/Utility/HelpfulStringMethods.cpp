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
#include <opm/simulators/geochemistry/Utility/HelpfulStringMethods.hpp>

std::vector<std::string>::iterator partition_vector_by_key(std::vector<std::string>& unsorted_vector, const std::string& key) {

    const std::string KEY = to_upper_case(key);

    return std::stable_partition(unsorted_vector.begin(), unsorted_vector.end(), 
        [KEY](const std::string& str) {
            return to_upper_case(str).find(KEY) != std::string::npos;
        });
}

std::vector<std::string> read_lines_from_stream(std::istream& inputStream, const std::string& endMarker)
{
    std::vector<std::string> lines_read;
    std::string buf;
    while (std::getline(inputStream, buf))
    {
        if(buf.find(endMarker) != std::string::npos){
            break; // found end marker
        }
        else if(buf.find("#") == std::string::npos){
            lines_read.push_back(buf);
        }
    }
    
    std::getline(inputStream, buf);
    
    return lines_read;
}

std::string remove_charge_from_species_name(const std::string& species_name){
    std::string s1 = split_by_delimiter(species_name, "+")[0];
    return split_by_delimiter(s1, "-")[0];
}

std::string substring_between(const std::string& str, const std::string& left_delimiter, const std::string& right_delimiter)
{
    std::size_t idx_left = str.find(left_delimiter);
    if (idx_left == std::string::npos) idx_left = 0;
    else idx_left += left_delimiter.size();

    std::size_t idx_right = str.find(right_delimiter);
    const std::size_t noChars = (idx_right == std::string::npos) ? idx_right : idx_right - idx_left;

    return str.substr(idx_left, noChars);
}

std::pair<std::string,std::string> get_first_word(const std::string& line)
{
    // returns the first word, and the rest of the line
    std::pair<std::string, std::string> key_val={"",""};
    std::vector<std::string> Tab = split_string(line);
    std::vector<std::string>::iterator it;
    if (Tab.empty())
        return key_val;
    else if (Tab.size() == 1)
    {
        key_val.first = trim_left_right(Tab[0]);
        return key_val;
    }
    else
    {
        key_val.first = trim_left_right(Tab[0]);
        for (it = Tab.begin()+1; it != Tab.end(); ++it)
            key_val.second += (*it) + " ";
        key_val.second = trim_left_right(key_val.second);
        return key_val;
    }

}

int in_list(const std::string& name, const std::vector<std::string>& string_list)
{
    // checks if an element is in the list case insensitive, returns -1 if not found
    if (string_list.empty())
        return -1;
    else
    {
        std::string NAME = to_upper_case(name);
        std::vector<std::string>::const_iterator it;
        for (it = string_list.begin(); it != string_list.end(); ++it)
        {
            std::string ELEMENT = to_upper_case(*it);
            if (ELEMENT == NAME)
                return std::distance(string_list.begin(), it);
        }
        return -1;
    }
}

std::vector<std::string> split_by_delimiter(const std::string& string_to_split,
                                            char delimiter)
{
    std::vector<std::string> ret;

    std::stringstream ss(string_to_split);
    std::string str;
    while (std::getline(ss, str, delimiter))
    {
        ret.push_back(str);
    }

    return ret;
}

std::vector<std::string> split_by_delimiter(const std::string& string_to_split,
                                            const std::string& delimiter)
{
    std::vector<std::string> list_of_strings = std::vector<std::string>();
    
    size_t pos = 0;
    std::string token;
    
    std::string str(string_to_split);
    while ((pos = str.find(delimiter)) != std::string::npos)
    {
        token = str.substr(0, pos);
        str.erase(0, pos + delimiter.length());
        list_of_strings.push_back(token);
    }
    
    if(!str.empty()) list_of_strings.push_back(str);
    
    if(list_of_strings.empty())
    {
        list_of_strings.push_back(string_to_split);
    }
    
    return list_of_strings;
}

std::vector<std::string> split_string(const std::string& string_to_split)
{
    // Splits a string by white space, and stores the individual tokens (if any) in a
    // vector, which is subsequently returned.
    std::vector<std::string> list_of_strings = std::vector<std::string>();
    
    std::istringstream iss(string_to_split);
    do
    {
        std::string s;
        iss >> s;
        if(!s.empty())
            list_of_strings.push_back(s);
        
    }while(iss);
    
    return list_of_strings;
}

bool string_contains(std::string_view str, std::string_view sub_str)
{
    return str.find(sub_str) != std::string::npos;
}

bool starts_with(const std::string& str, const std::string& start)
{
    if (start.size() > str.size())
    {
        return false;
    }
    return std::equal(start.begin(), start.end(), str.begin());
}

bool ends_with(const std::string & str, const std::string& ending)
{
    if (ending.size() > str.size()){
        return false;
    }
    return std::equal(ending.rbegin(), ending.rend(), str.rbegin());
}

/** @return True, if the string is empty after trimming it. */
bool has_content(const std::string& str)
{
    if(trim_left_right(str).empty())    return false;
    return true;
}

std::string to_lower_case(const std::string& string_to_convert)
{
    // Not perfect, will not work with multi-byte characters...
    std::string new_str;
    new_str.resize(string_to_convert.size());
    std::transform(string_to_convert.begin(), string_to_convert.end(), new_str.begin(), ::tolower);
    return new_str;
}



std::string to_upper_case(const std::string& string_to_convert)
{
    // Not perfect, will not work with multi-byte characters...
    std::string new_str;
    new_str.resize(string_to_convert.size());
    std::transform(string_to_convert.begin(), string_to_convert.end(), new_str.begin(), ::toupper);
    return new_str;
}

std::string trim_left(const std::string& str, const std::string& t)
{
    std::string s(str);
    s.erase(0, s.find_first_not_of(t));
    return s;
}

std::string trim_right(const std::string& str, const std::string& t)
{
    std::string s(str);
    s.erase(s.find_last_not_of(t) + 1);
    return s;
}

std::string trim_left_right(const std::string& str, const std::string& t)
{
    return trim_left(trim_right(str, t), t);
}

std::vector<std::string> to_lower_case(const std::vector<std::string>& strings_to_convert){
    auto copy_of_vector = strings_to_convert;
    for(std::size_t i=0; i < copy_of_vector.size(); ++i)
    {
        copy_of_vector[i] = to_lower_case(copy_of_vector[i]);
    }
    return copy_of_vector;
}

std::vector<std::string> to_upper_case(const std::vector<std::string>& strings_to_convert)
{
    auto copy_of_vector = strings_to_convert;
    for(std::size_t i=0; i < copy_of_vector.size(); ++i)
    {
        copy_of_vector[i] = to_upper_case(copy_of_vector[i]);
    }
    return copy_of_vector;
}

std::string remove_suffixes(const std::string& str) {
    return str.substr(0, str.find("."));
}

int index_of_string_in_vector(const std::string& input_string, const std::vector<std::string>& vec)
{
    for (auto it = vec.cbegin(); it != vec.cend(); ++it)
    {
        if (*it == input_string) return it - vec.begin();
    }
    return -1;
}

int index_of_string_in_vector_upper_case(const std::string& input_string, const std::vector<std::string>& vec)
{
    for (auto it = vec.cbegin(); it != vec.cend(); ++it)
    {
        if ( to_upper_case(*it) == to_upper_case(input_string) ) return it - vec.begin();
    }
    return -1;
}
