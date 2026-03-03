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
#ifndef SERIALIZE_FOR_TESTING_IS_DEF_H
#define SERIALIZE_FOR_TESTING_IS_DEF_H

#include <cereal/types/map.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/utility.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/archives/binary.hpp>

/**
 * For storing numerical output in a simple data structure, mainly for use in automatic regression testing.
 */
struct BasVecInfo
{
    std::map<std::string, double> key_solution_properties_;
    std::map<std::pair<std::string, std::string>, double> species_properties_;
    
    template <class Archive>
    void save(Archive& ar) const
    {
        ar(key_solution_properties_, species_properties_);
    }
        
    template <class Archive>
    void load(Archive& ar)
    {
        ar(key_solution_properties_, species_properties_);
    }
    
};

/* Helper struct that can be used to store effluent ion data at a given moment in time. */
struct EffluentIonData
{
    std::vector<std::string> names_;
    std::vector<double> values_;
    
    template <class Archive>
    void save(Archive& ar) const
    {
        ar(names_, values_);
    }
        
    template <class Archive>
    void load(Archive& ar)
    {
        ar(names_, values_);
    }
    
};


#endif
