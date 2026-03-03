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
#ifndef GEOCHEM_INPUT_T_IS_DEF_H
#define GEOCHEM_INPUT_T_IS_DEF_H

#include <memory>

#include <opm/simulators/geochemistry/Utility/HelpfulStringMethods.hpp>
#include <opm/simulators/geochemistry/IO/InputReader.hpp>

class input_t
{
    
public:
    
    input_t();
    
    std::unique_ptr<InputReader> header_;
    
    void read_header(std::string fname);
    std::string get_file_name() const;
    
private:
    
    bool has_defaults_;
    std::string fname_;
    
    std::vector<std::string> allowed_types_;
    std::vector<std::string> default_names_;
    std::vector<std::string> default_type_;
    std::vector<std::string> end_keywords_;
    std::vector<std::string> block_keywords_;
    std::map<std::string, std::string> simple_key_value_pairs_ = std::map<std::string, std::string>();
    
    // Helper functions
    void initialize_defaults();
    int set_default_key_values(const std::string& name, const std::string& value, const std::string& type);
    
};

#endif  // GEOCHEM_INPUT_T_IS_DEF_H
