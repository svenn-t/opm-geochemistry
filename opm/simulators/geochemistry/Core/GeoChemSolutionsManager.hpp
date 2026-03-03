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
#ifndef GEO_CHEMICAL_SOLUTIONS_MANAGER_HPP
#define GEO_CHEMICAL_SOLUTIONS_MANAGER_HPP

#include <vector>

#include <opm/simulators/geochemistry/Common/ChemGlobal.h>
#include <opm/simulators/geochemistry/SplayTree/SplayTree.h>
#include <opm/simulators/geochemistry/SplayTree/SplayTreeKeyCreator.hpp>

class BasVec;
class InitChem;

class GeochemicalSolutionsManager
{
  public:

    GeochemicalSolutionsManager();
    ~GeochemicalSolutionsManager();

    /**
     * @return The number of values that the solver needs to store
     *         (regardless of whether the splay tree look-up is used or not).
     * */
    std::size_t init(int splay_tree_resolution,
                     const std::vector<std::string>& basis_species,
                     const std::vector<std::string>& minerals,
                     bool has_ion_exchange,
                     bool has_surface_complexes);

    int splayTreeResolution() const;

    std::vector<int> createKey(const double* c_bas,
                               double dt,
                               double temperature,
                               const double* log_a_min) const;

    void print_information() const;

    SplayTree* get(int key);

    BasVec** get_vchem(const double* c_bas,
                       const double* c_min,
                       InitChem* ICS_full);

  private:

    std::vector<std::string> specie_names_;
    std::vector<std::string> mineral_names_;

    int no_basis_combinations_{0};
    int no_mineral_combinations_{0};

    // For retrieving the geochemical database and the set of basis species
    InitChem** ICS_key_{nullptr};
    BasVec*** Vchem_key_{nullptr};

    // Anonymous struct for workspace arrays
    struct
    {
        std::vector<int> local_binary_key_;
        std::vector<int> global_binary_key_;
    } workspace;

    SplayTreeKeyCreator sp_key_creator_{};
    std::vector<SplayTreeUniquePtr> sp_;  // Indexed by mineral combination

  private:

    void construct_BasVec_objects_for_InitChem(int max_number_of_mineral_combinations,
                                               BasVec** Vchem_in,
                                               InitChem* ICS_in
    );
};

#endif // GEO_CHEMICAL_SOLUTIONS_MANAGER_HPP
