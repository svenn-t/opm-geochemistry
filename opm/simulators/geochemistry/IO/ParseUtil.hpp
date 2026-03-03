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
#ifndef IO_PARSE_UTIL_IS_INCLUDED
#define IO_PARSE_UTIL_IS_INCLUDED

#include <iostream>
#include <istream>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>

#include <optional>

/* Parses doubles until reaching something that cannot be converted to a double. */
std::vector<double> parse_doubles(const std::string& str);

template<typename T>
using base_type = typename std::remove_cv<typename std::remove_reference<T>::type>::type;

// Note: We have only actually tested this function on ints and doubles (which are the types it is supposed to use work for).
template<class T, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
std::optional<T> parse_value_or_error(std::istream& inputStream, const std::string& err_msg)
{
    base_type<T> value{};
    if(! (inputStream >> value) )
    {
        std::cout << err_msg << "\n";
        return std::nullopt;
    }
    return std::make_optional(value);
}

#endif
