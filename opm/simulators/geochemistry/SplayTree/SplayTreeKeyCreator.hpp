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
#ifndef SPLAY_TREE_KEY_CREATOR_IS_DEF_HPP
#define SPLAY_TREE_KEY_CREATOR_IS_DEF_HPP

#include <vector>

/** Also calculates the number of values to store in each splay tree node. */
class SplayTreeKeyCreator
{
public:
    SplayTreeKeyCreator();
    SplayTreeKeyCreator(
        int splay_tree_resolution,
        int max_no_basis_species,
        int max_no_minerals,
        bool has_ion_exchange,
        bool has_surface_complexes
    );

    int splayTreeResolution() const;
    int keySize() const;
    int numberOfValuesToStore() const;

    std::vector<int> createKey(const double* c_bas, double dt, double temperature, const double* log_a_min) const;

private:
    int no_basis_species_{0};
    int no_minerals_{0};

    int no_stored_values_{0};

    int num_key_log_{0};    // Number of key values with logarithmic binning  == no_basis_species_ + 1
    int num_key_lin_{0};    // Number of key values with linear binning       == no_minerals_ + 1
    int num_key_{0};        // == num_key_log_ + num_key_lin_

    int splay_tree_resolution_{0};  // == 0 means no interpolation
    double log_bin_width_{10.0};
    double inv_log10_of_log_bin_width_{1.0};
    double lin_interval_length_{1.0};

    std::vector<double> val_scale_log_;  // Scaling factor logarithmically binned values (default: 1.0)
    std::vector<double> val_scale_log_inv_;
    std::vector<double> val_scale_lin_;  // Scaling factor for linearly binned values (default: 1.0)
    std::vector<double> val_scale_lin_inv_;
};

#endif  // SPLAY_TREE_KEY_CREATOR_IS_DEF_HPP