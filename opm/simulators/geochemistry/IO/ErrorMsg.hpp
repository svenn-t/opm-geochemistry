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
#ifndef ERROR_MSG_IS_DEF_HPP
#define ERROR_MSG_IS_DEF_HPP

#include <opm/simulators/geochemistry/IO/fmt_include.hpp>

template<typename... Args>
static void error(fmt::format_string<Args...> message, Args&&... args)
{
    fmt::print("ERROR: ");
    fmt::print(message, std::forward<Args>(args)...);
    fmt::print("\n");
}

template<typename... Args>
static void error_and_exit(fmt::format_string<Args...> message, Args&&... args)
{
    error(message, std::forward<Args>(args)...);
    std::exit(0);
}

template<typename... Args>
static void issue_warning(const char* message, Args... args)
{
    fmt::print("WARNING: ");
    fmt::print(message, args...);
    fmt::print("\n");
}

#endif