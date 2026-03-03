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
/**
* @file
*
* A collection of utility functions, many for working with templated containers.
*/
#ifndef MISC_UTILITY_FUNC_IS_DEF_H
#define MISC_UTILITY_FUNC_IS_DEF_H

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <numeric>
#include <set>
#include <vector>

#include <optional>

template<class Real>
// Remark: It might make more sense to take the input arguments by value here? (they will usually be doubles)
static constexpr Real restrictByMaximumAbsoluteValue(const Real& proposedValue, const Real& maxAbsoluteValue)
{
    assert(maxAbsoluteValue > 0.0);
    if (proposedValue > maxAbsoluteValue) return maxAbsoluteValue;
    else if (proposedValue < -maxAbsoluteValue) return -maxAbsoluteValue;
    else return proposedValue;
}

// ***************************************************************************
//                   FUNCTIONS FOR WORKING WITH VARIADIC TEMPLATES
// ***************************************************************************

/**
 * @return True iff the first input argument is also among the rest.
 * */
template<typename T, typename ... Values>
bool anyOf(T elem, Values ...value)
{
    return (... || (elem == value));
}

template<typename T>
static void printParameterPack(T last)
{
    std::cout << last << "\n";
}

template <typename T, typename... args>
static void printParameterPack(T first, args... rest)
{
    std::cout << first << "\n";
    printParameterPack(rest...);
}

// ***************************************************************************
// FUNCTIONS OPERATING ON STANDARD LIBRARY CONTAINERS OF A KNOWN TYPE
// ***************************************************************************

template<typename dummy=int> // Just so we can compile it without a cpp file...
bool vector_has_infinite_numbers(const std::vector<double>& vec)
{
    return std::any_of(vec.begin(),
                       vec.end(),
                       [](double x){return std::isinf(x); });
}

// ***************************************************************************
// FUNCTIONS OPERATING ON STANDARD LIBRARY CONTAINERS WITH TEMPLATED TYPE
// ***************************************************************************

template<typename T>
void append_to_vector(std::vector<T>& vec_before, const std::vector<T>& vec_to_add)
{
    if (!vec_to_add.empty())
    {
        vec_before.insert(vec_before.end(), vec_to_add.begin(), vec_to_add.end());
    }
}

enum class Comparison{ EQUAL, LEFT_LARGER, LEFT_SMALLER };

template<typename T>  // Objects of type T are comparable
Comparison compare_vectors(const std::vector<T>& vec1, const std::vector<T>& vec2)
{
    assert(vec1.size() == vec2.size());

    for (std::size_t i = 0; i < vec1.size(); ++i)
    {
        if (vec1[i] > vec2[i]) return Comparison::LEFT_LARGER;
        if (vec2[i] > vec1[i]) return Comparison::LEFT_SMALLER;
    }
    return Comparison::EQUAL;
}

template<typename T>
void fill_beginning_of_vector_from_vector(const std::vector<T>& from, std::vector<T>& to, std::size_t max_no_elements)
{
    std::size_t no_elements = std::min(to.size(), max_no_elements);
    no_elements = std::min(from.size(), no_elements);  // cannot copy more elements than are actually there

    for(std::size_t i=0; i < no_elements; ++i)
    {
        to[i] = from[i];
    }
}

template<typename T>
void fill_beginning_of_vector_from_vector(const std::vector<T>& from, std::vector<T>& to)
{
    fill_beginning_of_vector_from_vector(from, to, to.size());
}

template<typename SequenceContainer, typename ValueType = typename SequenceContainer::value_type>
void fill_value(SequenceContainer& container, ValueType value)
{
    std::fill(container.begin(), container.end(), value);
}

template<typename SequenceContainer, typename  ValueType = typename SequenceContainer::value_type>
void fill_zero(SequenceContainer& container)
{
    fill_value(container, static_cast<ValueType>(0.0));
}

/** The caller should obviously make sure that the size of the array is big enough (>= max_no_values). */
template<typename SequenceContainer, typename ValueType = typename SequenceContainer::value_type>
void fill_beginning_of_array_from_container(ValueType* array, const SequenceContainer& container, std::size_t max_no_values)
{
    const std::size_t no_values_to_fill = std::min(max_no_values, container.size());
    for(std::size_t i=0; i < no_values_to_fill; ++i)
    {
        array[i] = container[i];
    }
}

/**
 * Does NOT alter the size of the container. If the size is smaller than
 * no_values, only insert as many values as there are room for. */
template<typename SequenceContainer, typename ValueType = typename SequenceContainer::value_type>
void fill_beginning_of_container_from_array(SequenceContainer& container, const ValueType* array, std::size_t max_no_values)
{
    const std::size_t no_values_to_fill = std::min(max_no_values, container.size());
    for(std::size_t i=0; i < no_values_to_fill; ++i)
    {
        container[i] = array[i];
    }
}

/**
 * Moves elements from a source vector into a destination vector (at the end).
 * The source vector is left in undefined but safe-to-destruct state.
 */
template<typename T>
void move_elements_from_vector_to_vector(std::vector<T>& src, std::vector<T>& dest)
{
    const auto minimum_capacity = dest.size() + src.size();
    if(dest.capacity() < minimum_capacity) dest.reserve(minimum_capacity);

    dest.insert(dest.end(),
                std::make_move_iterator(src.begin()),
                std::make_move_iterator(src.end())
    );
}

template<typename T>
void multiply_vector_by_value(std::vector<T>& vec, T value)
{
    std::transform(vec.begin(), vec.end(), vec.begin(), [&value](auto& v){ return v*value; });
}

template<typename T>
std::vector<T> copy_vector_and_multiply_by_value(const std::vector<T>& vec, T value)
{
    std::vector<T> new_vec = vec; // Make copy
    multiply_vector_by_value(new_vec, value);
    return new_vec;
}

template<typename T>
void copy_from_vector_to_vector(const std::vector<T>& vec_from, std::vector<T>& vec_to)
{
    assert(vec_to.size() >= vec_from.size());
    std::copy_n(vec_from.begin(), vec_from.size(), vec_to.begin());
}

template<typename T>
std::size_t number_of_non_null_values(const std::vector<std::optional<T>>& vec)
{
    return std::count_if(vec.cbegin(), vec.cend(), [](const std::optional<T>& element){ return element.has_value(); });
}

template<typename T>
void print_vector_of_vectors(const std::vector<std::vector<T>>& vec, std::ostream& out=std::cout)
{
    for(std::size_t rowIndex=0; rowIndex < vec.size(); ++ rowIndex)
    {
        const auto noColumns = vec[rowIndex].size();

        if(noColumns > 0) out << vec[rowIndex][0];
        for(std::size_t columnIndex=1; columnIndex < noColumns; ++ columnIndex)
        {
            out << "\t" << vec[rowIndex][columnIndex];
        }
        out << "\n";
    }
}

template <typename T>
void remove_element_from_vector(std::vector<T>& vec, std::size_t pos)
{
    auto it = vec.begin();
    std::advance(it, pos);
    vec.erase(it);
}

template <typename T>
std::vector<std::size_t> sort_vector_indices(const std::vector<T>& vec,bool reverse=true)
{
    // Original index locations
    std::vector<std::size_t> indices(vec.size());
    std::iota(indices.begin(), indices.end(), 0);

    // Sort indices based on values in vec
    if(reverse)
        std::sort(indices.begin(), indices.end(),
                  [&vec](std::size_t i1, std::size_t i2) {return vec[i1] > vec[i2]; });
    else
        std::sort(indices.begin(), indices.end(),
            [&vec](std::size_t i1, std::size_t i2) {return vec[i1] < vec[i2]; });

    return indices;
}

template<typename T>
bool vector_has_duplicates(const std::vector<T>& vec)
{
    std::set<T> set_of_vec(vec.begin(), vec.end());
    return set_of_vec.size() != vec.size();
}

/** First creates a sorted copy of the input vector, then returns the copy after removing duplicates. */
template<typename T>
std::vector<T> remove_duplicates(const std::vector<T>& vec)
{
    std::vector<T> ret = vec;
    std::sort(ret.begin(), ret.end());
    ret.erase(std::unique(ret.begin(), ret.end()), ret.end());  // Note: Unique only works for consecutive duplicates
    return ret;
}

template<typename T>
std::vector<T> find_duplicates(const std::vector<T>& vec)
{
    std::vector<T> duplicates;
    if(vec.size() <= 1){
        return duplicates; // there cannot be duplicates in this case
    }

    auto sorted_vec = vec;
    std::sort(sorted_vec.begin(), sorted_vec.end());

    // Note: I am sure this could be implemented much more elegantly...
    std::size_t no_repeats = 0;
    for(std::size_t i=1; i < sorted_vec.size(); ++i){

        if(sorted_vec[i] == sorted_vec[i-1]){
            ++no_repeats;
        }

        if(no_repeats >= 1 && i == sorted_vec.size() - 1){
            duplicates.push_back(sorted_vec[i]);
        }
        else if(no_repeats >= 1){
            duplicates.push_back(sorted_vec[i-1]);
            no_repeats = 0;
        }

    }

    return duplicates;
}

/** Returns the index of the first occurrence of the element in the container,
 *  and -1 if the element does not exist.*/
template<typename SequenceContainer, typename T = typename SequenceContainer::value_type>
int index_of_element_in_container(const T& elem, const SequenceContainer& container)
{
    auto it = std::find(container.cbegin(), container.cend(), elem);
    if(it != container.cend()) return static_cast<int>(it - container.cbegin());
    else return -1;
}

template<typename T>
int index_of_element_in_set(const T& elem, const std::set<T>& set)
{
    const auto it = set.find(elem);
    if(it != set.cend())
    {
        return static_cast<int>(std::distance(set.cbegin(), it));
    }
    return -1;
}

template<typename SequenceContainer, typename T = typename SequenceContainer::value_type>
bool element_is_in_container(const T& elem, const SequenceContainer& container)
{
    const auto it = std::find(container.cbegin(), container.cend(), elem);
    return (it != container.cend());
}

// ****************************************************************************************
//                   FUNCTIONS OPERATING (ONLY) ON DYNAMIC ARRAYS
// ****************************************************************************************

template<typename T>
std::vector<T> array_as_vector(const T* array, int array_size)
{
    std::vector<T> vec;
    vec.resize(array_size);
    for(int i=0; i < array_size; ++i) vec[i] = array[i];
    return vec;
}

/* Returns -1 if the element is not found (exactly). */
template<typename T>
int index_of_element_in_array(T elem, const T* array, int array_size)
{
    for(int i=0; i < array_size; ++i)
    {
        if(elem == array[i])
        {
            return i;
        }
    }
    return -1;
}

template<typename T>
void insert_into_first_element_of_dynamic_array(T* array, const T elem, int array_size)
{
    for (int i = 1; i < array_size; ++i){
        array[array_size - i] = array[array_size - 1 - i];
    }
    array[0] = elem;
}

/** Overloaded version of the corresponding function for std::vectors. */
template<typename T>
std::size_t number_of_non_null_values(const T* array, int array_size)
{
    std::size_t count = 0;
    for(int i=0; i < array_size; ++i)
    {
        if(array[i]) ++count;
    }
    return count;
}

#endif
