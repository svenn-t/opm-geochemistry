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
#ifndef INPUT_FILE_HANDLER_HPP
#define INPUT_FILE_HANDLER_HPP
#include <opm/simulators/geochemistry/Utility/HelpfulStringMethods.hpp>

#include <fstream>
#include <string>

/** Helper class to create stream from file, and quit if file not opened.*/
class InputFileHandler
{
public:
    explicit InputFileHandler(std::string input_file_name);

    std::string getName() { return case_name_; }
    std::ifstream& getStream() { return input_stream_; }

private:

    std::string input_file_name_{};
    std::string case_name_{};
    std::ifstream input_stream_{};
    void fail() const;
};

#endif //INPUT_FILE_HANDLER_HPP
