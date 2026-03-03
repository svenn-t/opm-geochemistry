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
#include <opm/simulators/geochemistry/IO/input_m.h>

input_t::input_t()
: header_(nullptr)
, has_defaults_(false)
, fname_()
, allowed_types_()
, default_names_()
, default_type_()
, end_keywords_()
, block_keywords_()
, simple_key_value_pairs_()
{
    initialize_defaults();
}

void input_t::initialize_defaults(){
    
    // First, reset everything...
    fname_.clear();
    allowed_types_.clear();
    default_names_.clear();
    default_type_.clear();
    end_keywords_.clear();
    block_keywords_.clear();
    simple_key_value_pairs_.clear();
    
    // Default values and parameters for equilibrium solver and 1D solver
    allowed_types_ = { "FLOAT", "INT" };

    set_default_key_values("equilibrate","0","INT");
    
    set_default_key_values("Temp", "25.0", "FLOAT");
    set_default_key_values("Pres", "10e5","FLOAT");

    set_default_key_values("NoBlocks","10", "INT");
    set_default_key_values("VolRate","0.1", "FLOAT"); // Flow rate ml/min
    set_default_key_values("Volume","25.0", "FLOAT"); // Total pore volume (mL)
    set_default_key_values("Porosity", "1.0", "FLOAT");
    set_default_key_values("Tf","24.0", "FLOAT");     // Simulation time (hours)

    set_default_key_values("WriteBlock","0", "INT");
    set_default_key_values("WriteEndBlock","0", "INT");
    set_default_key_values("WriteNetCDF","0",  "INT");

    set_default_key_values("Imp", "1", "INT"); // Implicit or explicit integration?
    set_default_key_values("Flush", "0.5", "FLOAT");   // Volume flushed out of each cell

    set_default_key_values("Debug","0", "INT");
    set_default_key_values("Interpolate","0",  "INT");
    set_default_key_values("Serialize","0",  "INT");

    end_keywords_ = { R"(/end)", "end", R"(/ end)" };
    block_keywords_ = { "CHEMTOL","INJECT" };
    has_defaults_ = true;
}

void input_t::read_header(const std::string fname)
{
    if(!has_defaults_) initialize_defaults();
    has_defaults_ = false;
    
    fname_ = fname;
    header_ = std::make_unique<InputReader>(simple_key_value_pairs_, false);  // case-insensitive keywords
    
    for (auto it = block_keywords_.cbegin(); it != block_keywords_.cend(); ++it)
    {
        header_->define_block_keyword(*it, end_keywords_);
    }
    std::ifstream infile_stream(fname, std::ios::binary);
    header_->read_key_value_pairs(infile_stream);
    
    for (auto it = default_names_.cbegin(); it != default_names_.cend(); ++it)
    {
        // default values are overwritten if new values are read
        simple_key_value_pairs_[*it] = header_->get_simple_keyword_value(*it);
    }
}

std::string input_t::get_file_name() const{
    return fname_;
}

int input_t::set_default_key_values(const std::string &name, const std::string &value, const std::string &type)
{
    const std::string TYPE = to_upper_case(type);
    
    if (std::count(allowed_types_.begin(), allowed_types_.end(), TYPE))
    {
        default_names_.push_back(name);
        default_type_.push_back(TYPE);
        simple_key_value_pairs_[name] = value;
        return 1;
    }
    else
    {
        std::cout << " input: no type named " << type << std::endl;
        return 0;
    }
}
