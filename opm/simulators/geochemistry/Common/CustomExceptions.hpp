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
#ifndef CUSTOM_EXCEPTIONS_IS_DEF
#define CUSTOM_EXCEPTIONS_IS_DEF

#include <exception>
#include <string>

class IntegrandIsNotWellDefined : public std::runtime_error
{

public:
    IntegrandIsNotWellDefined(const std::string& error_msg) : std::runtime_error(error_msg) { };
};

class MissingKeywordException : public std::runtime_error
{

public:
    MissingKeywordException(const std::string& error_msg) : std::runtime_error(error_msg) { };
};

class ModelDoesNotExistException : public std::runtime_error
{

public:
    ModelDoesNotExistException(const std::string& error_msg) : std::runtime_error(error_msg) { };
};

class NotInitializedException : public std::runtime_error
{

public:
    NotInitializedException(const std::string& error_msg) : std::runtime_error(error_msg) { };
};

class InvalidInputException : public std::logic_error
{

public:
    InvalidInputException(const std::string& error_msg): std::logic_error(error_msg) { };
};

class MultipleInitializationException : public std::logic_error
{

public:
    MultipleInitializationException(const std::string& error_msg): std::logic_error(error_msg) { };
};

class NotImplementedException : public std::logic_error
{

public:
    NotImplementedException(const std::string& error_msg): std::logic_error(error_msg) { };
};


#endif
