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
#ifndef GEO_LOGGER_IS_DEF_HPP
#define GEO_LOGGER_IS_DEF_HPP

#include <cassert>
#include <cstdio>
#include <iostream>
#include <memory>
#include <utility>
#include <vector>

#include <opm/simulators/geochemistry/IO/fmt_include.hpp>

struct file_closer
{
    void operator()(std::FILE* fp){
        std::fclose(fp);
    }
};
using uniqueFilePtr = std::unique_ptr<std::FILE, file_closer>;

enum class LogLevel
{
    // Small letters to avoid potential conflicts with macros...
    debug = 0,
    info = 1,
    warning = 2,
    error = 3,
    critical = 4,
    none = 5 // must be the last
};


class GeoLogger
{
public:

    GeoLogger(LogLevel level=LogLevel::error);

    void setLogLevel(LogLevel level);
    void enableFileOutput(const char* fileName="tmp_logfile.txt");
    void disableFileOutput();

    template<typename... Args>
    void info(const char* message, Args... args)
    {
        if(logLevel_ <= LogLevel::info) log("INFO", message, args...);
    }

    template<typename... Args>
    void warning(const char* message, Args... args)
    {
        if(logLevel_ <= LogLevel::warning) log("WARNING", message, args...);
    }

    template<typename... Args>
    void error(const char* message, Args... args)
    {
        if(logLevel_ <= LogLevel::error) log("ERROR", message, args...);
    }

    template<typename... Args>
    void critical(const char* message, Args... args)
    {
        if(logLevel_ <= LogLevel::critical) log("CRITICAL ERROR", message, args...);
    }

    template<typename... Args>
    void debug(const char* message, Args... args)
    {
        if(logLevel_ <= LogLevel::debug) log("DEBUG", message, args...);
    }

private:
    LogLevel logLevel_{LogLevel::none};

    bool writeOutputToFile_{false};
    const char* filePath_{nullptr};
    uniqueFilePtr filePtr_{nullptr};

private:

    template<typename... Args>
    void log(const char* header, const char* message, Args... args)
    {
        // Note: fmt::runtime needed for C++20...
        //
        FILE* f = filePtr_.get();
        if(f)
        {
            assert(writeOutputToFile_);
            fmt::print(f, fmt::runtime("{}: "), header);
            fmt::print(f, fmt::runtime(message), args...);
            fmt::print(f, "\n");
        }
        else
        {
            fmt::print(fmt::runtime("{}: "), header);
            fmt::print(fmt::runtime(message), args...);
            fmt::print("\n");
        }
    }
};

#endif  // GEO_LOGGER_IS_DEF_HPP