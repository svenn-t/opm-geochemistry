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
#include <opm/simulators/geochemistry/IO/FileFunctions.hpp>

/**
 * Opens a file and checks if it is open correctly. If not, exits the program.
 *
 * @param [in] filename Name of file.
 * @param [in] mode The file mode (read, write, etc.)
 * @return Pointer to the opened file.
*/
FILE* my_fopen(const char* filename, const char* mode)
{
    FILE *fptr; // file pointer

    fptr = fopen(filename, mode);
    if (fptr == (FILE*)nullptr)
    {
        printf("error opening %s\n",filename);
        std::exit(0);
    }
    return fptr;
}

/** @return New string after removing all instances of the escape character '\"' in the original string.*/
std::string remove_quotation_marks_from_filename(const std::string& str)
{
    std::string new_str(str);

    std::size_t n = new_str.find("\"");
    while (n != std::string::npos)
    {
        new_str.erase(n, 1);
        n = new_str.find("\"");
    }
    return new_str;
}
