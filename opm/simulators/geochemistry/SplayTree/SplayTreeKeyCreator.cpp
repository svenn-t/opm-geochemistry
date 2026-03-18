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
#include <opm/simulators/geochemistry/SplayTree/SplayTreeKeyCreator.hpp>

#include <cmath>
#include <opm/simulators/geochemistry/Common/ChemGlobal.h>

SplayTreeKeyCreator::SplayTreeKeyCreator() = default;

SplayTreeKeyCreator::SplayTreeKeyCreator(
    int splay_tree_resolution,
    int max_no_basis_species,
    int max_no_minerals,
    bool has_ion_exchange,
    bool has_surface_complexes
)
: no_basis_species_(max_no_basis_species)
, no_minerals_(max_no_minerals)
, num_key_log_(max_no_basis_species + 1)
, num_key_lin_(max_no_minerals + 1)
, num_key_(max_no_basis_species + max_no_minerals + 2)
, splay_tree_resolution_(splay_tree_resolution)
, log_bin_width_(std::pow(10.0, 1.0 / splay_tree_resolution))
, lin_interval_length_(1.0)
{
    inv_log10_of_log_bin_width_ = 1.0 / std::log10(log_bin_width_);

    val_scale_log_.resize(num_key_log_, 1.0);
    val_scale_lin_.resize(num_key_lin_, 1.0);

    val_scale_log_inv_.resize(num_key_log_);
    val_scale_lin_inv_.resize(num_key_lin_);

    // Set scaling values to find key
    for (int i = 0; i < num_key_log_; ++i)
    {
        val_scale_log_inv_[i] = 2.0 * log_bin_width_ / (val_scale_log_[i] * (log_bin_width_ + 1.0));
    }

    // Calculate the lin base
    for (int i = 0; i < num_key_lin_; ++i)
    {
        val_scale_lin_[i] *= lin_interval_length_;
        val_scale_lin_inv_[i] = 1. / val_scale_lin_[i];
    }

    // Finally, compute the number of values to store
    no_stored_values_ = max_no_basis_species + 1;  // Include pH.
    if(has_surface_complexes)
    {   // Two extra for surface charge and surface potential
        no_stored_values_ =  2 * max_no_basis_species + 3;
    }
    else if(has_ion_exchange)
    {
        no_stored_values_ = 2 * max_no_basis_species + 1;
    }
}

int SplayTreeKeyCreator::splayTreeResolution() const
{
    return splay_tree_resolution_;
}

int SplayTreeKeyCreator::keySize() const
{
    return num_key_;
}

int SplayTreeKeyCreator::numberOfValuesToStore() const
{
    return no_stored_values_;
}

std::vector<int> SplayTreeKeyCreator::createKey(const double* c_bas, double dt, double temperature, const double* log_a_min) const
{
    std::vector<int> new_key;
    new_key.resize(num_key_);

    // Note: Zero values are evaluated to -inf (e.g., relevant when dt==0).
    for (int i = 0; i < no_basis_species_; ++i)
    {
        new_key[i] = negative_floor_function(inv_log10_of_log_bin_width_ * log10(val_scale_log_inv_[i] * c_bas[i]));
    }
    new_key[no_basis_species_] = negative_floor_function(inv_log10_of_log_bin_width_ * log10(val_scale_log_inv_[num_key_log_-1] * dt));

    new_key[no_basis_species_ + 1] = negative_floor_function(val_scale_lin_inv_[0] * temperature + 0.5);

    for (int i = 0; i < no_minerals_; ++i)
    {
        new_key[no_basis_species_ + 2 + i] = negative_floor_function(val_scale_lin_inv_[i+1] * log_a_min[0] + 0.5);
    }

    return new_key;
}
