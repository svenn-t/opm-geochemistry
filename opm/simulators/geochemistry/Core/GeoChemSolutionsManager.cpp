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
#include <opm/simulators/geochemistry/Core/GeoChemSolutionsManager.hpp>

#include <utility>

#include <opm/simulators/geochemistry/Core/ChemBasVec.h>
#include <opm/simulators/geochemistry/Core/ChemInitChem.h>

GeochemicalSolutionsManager::GeochemicalSolutionsManager() = default;

GeochemicalSolutionsManager::~GeochemicalSolutionsManager()
{
    //  We must delete the BasVec objects first.
    if(Vchem_key_)
    {
        for (int i = 0; i < no_basis_combinations_; ++i)
        {
            if (Vchem_key_[i])
            {
                for (int j = 0; j < no_basis_combinations_; ++j)
                {
                    delete Vchem_key_[i][j];
                }
            }
            delete[] Vchem_key_[i];
        }
        delete[] Vchem_key_;
    }

    if (ICS_key_)
    {
        for (int i = 0; i < no_basis_combinations_; ++i)
        {
            delete ICS_key_[i];
        }
        delete[] ICS_key_;
    }
}

std::size_t GeochemicalSolutionsManager::init(int splay_tree_resolution,
                                   const std::vector<std::string>& basis_species,
                                   const std::vector<std::string>& minerals,
                                   bool has_ion_exchange,
                                   bool has_surface_complexes)
{
    const auto max_no_basis_species = static_cast<int>(basis_species.size());
    const auto max_no_minerals      = static_cast<int>(minerals.size());
    specie_names_                   = basis_species;
    mineral_names_                  = minerals;
    no_basis_combinations_          = number_of_combinations(max_no_basis_species);
    no_mineral_combinations_        = number_of_combinations(max_no_minerals);

    sp_key_creator_ = SplayTreeKeyCreator(
        splay_tree_resolution,
        max_no_basis_species,
        max_no_minerals,
        has_ion_exchange,
        has_surface_complexes
    );

    sp_.clear();
    if (splay_tree_resolution > 0)
    {   // Allocate memory for splay trees
        sp_.reserve(no_mineral_combinations_);
        for (int i_min = 0; i_min < no_mineral_combinations_; ++i_min)
        {
            sp_.emplace_back(std::make_unique<SplayTree>(
                i_min,
                sp_key_creator_.keySize(),
                sp_key_creator_.numberOfValuesToStore())
            );
        }
    }

    ICS_key_ = new InitChem*[no_basis_combinations_];
    Vchem_key_ = new BasVec**[no_basis_combinations_];

    for(int i = 0; i < no_basis_combinations_; ++i)
    {
        ICS_key_[i]     = nullptr;
        Vchem_key_[i]   = nullptr;
    }

    workspace.local_binary_key_.resize(max_no_basis_species, 0);
    workspace.global_binary_key_.resize(max_no_basis_species, 0);

    return sp_key_creator_.numberOfValuesToStore();
}

int GeochemicalSolutionsManager::splayTreeResolution() const
{
    return sp_key_creator_.splayTreeResolution();
}

std::vector<int> GeochemicalSolutionsManager::createKey(const double* c_bas, double dt, double temperature, const double* log_a_min) const
{
    return sp_key_creator_.createKey(c_bas, dt, temperature, log_a_min);
}

void GeochemicalSolutionsManager::print_information() const
{
    if (sp_key_creator_.splayTreeResolution() > 0)
    {
        const auto max_no_minerals = static_cast<int>(mineral_names_.size());

        auto tmp_binary_key = std::vector<int>(no_mineral_combinations_);
        for (int i = 0; i < no_mineral_combinations_; ++i)
        {
            if(sp_[i]->numberOfLookups() == 0)  continue;  // Only print out the mineral comb. that were used.

            std::cout << "SplayTree size: " << sp_[i]->numberOfTreeNodes() << ", number of successful lookups: " << sp_[i]->numberOfLookups() << "\n";

            convert_into_binary(i, tmp_binary_key.data(), max_no_minerals);
            std::cout << "SplayTree mineral combination:\n";
            for (int j = 0; j < max_no_minerals; ++j)
            {
                if (tmp_binary_key[j]) std::cout << mineral_names_[j] << "\t";
            }
            std::cout << "\n";
        }
    }
}

SplayTree* GeochemicalSolutionsManager::get(int key)
{
    SplayTree* pSplay_tree = sp_[key].get();
    assert(pSplay_tree); // Must exist!
    return pSplay_tree;
}

/**
 * In addition to constructing InitChem and BasVec objects, this function generates a unique key
 * for the current basis combination.
 *
 * Example:
 *  If Ca, Cl, H, HCO3, and Mg are possible ions, we also check if minerals are present and add basis sources at minerals.
 *  Since ctot contains the total concentration including surface species, we do not need to check for species at the surface.
 *
 *  Note: Could save space by not allocating full Vchem vector...
 *
 * @param [in] c_bas When called, contains aqueous concentrations + surface concentrations.
 * @param [in] c_min  Concentration of minerals in rock.
 * @param [in] ICS_full Contains the full list of possible species and minerals.
 * @return The full Vchem vector for a given basis combination and all possible mineral combinations in the rock.
 */
BasVec** GeochemicalSolutionsManager::get_vchem(const double* c_bas,
                                     const double* c_min,
                                     InitChem* ICS_full)
{
    auto cb_cp = array_as_vector(c_bas, ICS_full->size_basis_);

    for (int i = 0; i < ICS_full->size_basis_; ++i)
    {
        for (int imin = 0; imin < ICS_full->size_min_; ++imin)
        {
            // For each combination of basis (ion) and mineral species, need to
            // know how many moles of the ion is in one mole of the mineral...
            const auto stoichiometric_coeff = ICS_full->SM_mineral_->M_[ICS_full->pos_min_[imin]][ICS_full->pos_[i]];
            cb_cp[i] += std::fabs(stoichiometric_coeff * c_min[imin]);
        }
    }

    const int key_bas = create_unique_key_from_positive_elements(cb_cp.data(), ICS_full->size_basis_);

    if (!Vchem_key_[key_bas])
    {
        if(!ICS_key_[key_bas])  // Replaces the old (deleted) function get_ICS()
        {
            ICS_key_[key_bas] = InitChem::create_from_ICS_full(*ICS_full,
                                                               c_bas,
                                                               c_min,
                                                               cb_cp);
        }
        InitChem* ICSp = ICS_key_[key_bas];

        // Here, we allocate memory for each possible combination of minerals.
        Vchem_key_[key_bas] = new BasVec*[no_basis_combinations_];
        for (int i = 0; i < no_basis_combinations_; ++i)
        {
            Vchem_key_[key_bas][i] = nullptr;
        }
        construct_BasVec_objects_for_InitChem(ICS_full->size_min_, Vchem_key_[key_bas], ICSp);
    }

    return Vchem_key_[key_bas];
}

/**
 * The input InitChem object contains information about all minerals that are
 * allowed to form when using the given basis.
 *
 * As an example, suppose 3 minerals can form: anhydrite, calcite, and dolomite.
 * Then, there are 2^3=8 possible mineral combinations, each of which can be
 * given a binary representation:
 *
 *  (0 0 0) = no minerals.
 *  (0 0 1) = (0, 0 , dolomite) i.e. only dolomite.
 *  (0 1 0) = (0, calcite , 0) i.e. only calcite.
 *  (0 1 1) = (0, calcite , dolomite) i.e., calcite & dolomite, but no anhydrite.
 *
 *  ... etc ....
 *
 * Note however that before computing a unique integer key for each BasVec object,
 * these binary representations must be translated to a binary representation
 * with N bits, where N >= 3 is the maximum number of minerals that can occur
 * during the entire simulation (i.e., N=ICS_full->size_min >= ICS_in->size_min).
 * This is because the solver always computes keys relative to the full set of minerals.
 *
 * A detailed example to illustrate the problem:
 *
 *  Suppose the following 11 ions can potentially exist: Ba, Br, Ca, Cl, H, HCO3, K, Mg, Na, SO4, and Sr.
 *  Therefore, there are 2^11 = 2048 possible combinations (key_bas=0, 1, ..., 2047)
 *  Suppose also that the following 4 minerals can exist: Anhydrite, Barite, Calcite, and Dolomite,
 *  i.e., there are 2^4=16 combinations (key=0, 1, ..., 15).
 *
 *  Suppose we want to construct all BasVec objects for the specific basis that contains all ions except barium.
 *  If the ions are sorted alphabetically, this corresponds to key_bas=1023 (binary representation: 01111111111).
 *  Without barium present, the only possible minerals are anhydrite, calcite, and dolomite, hence we get
 *  the 8 binary representations listed above. However, we still need to account for the fact that barite is the
 *  second mineral in the full list, thus we must change the binary representations by adding a column of zeros
 *  after the first column:
 *
 *  (0 0 0 0)
 *  (0 0 0 1)
 *  (0 0 1 0)
 *  (0 0 1 1)
 *  ... etc ....
 *
 *  In general: For each mineral in the full set that is missing from the basis under consideration,
 *  we have to insert a column of zeros at the appropriate column index.
 *
 * @param [in] max_number_of_mineral_combinations
 * @param [in,out] Vchem_in: Array of BasVec* objects in which to insert new ones, one for each possible mineral combination.
 * @param [in] ICS_in: Holds information about basis species and possible minerals.
 */
void GeochemicalSolutionsManager::construct_BasVec_objects_for_InitChem(
int max_number_of_mineral_combinations,
BasVec** Vchem_in,
InitChem* ICS_in
)
{
    const int numberOfLocalMineralCombinations = number_of_combinations(ICS_in->size_min_);

    for (int i = 0; i < numberOfLocalMineralCombinations; ++i)
    {
        fill_zero(workspace.local_binary_key_);
        fill_zero(workspace.global_binary_key_);

        convert_into_binary(i, workspace.local_binary_key_.data(), ICS_in->size_min_);

        for (int j = 0; j < ICS_in->size_min_; ++j)
        {
            // When setting mineral concentrations in ICS_in, the actual values are immaterial,
            // the constructor of the new BasVec object only needs to know whether they are
            // positive or not. The actual values will be provided later by the solver.
            workspace.global_binary_key_[ICS_in->krmap_[j]] = workspace.local_binary_key_[j];

            const auto bkey_j = static_cast<double>(workspace.local_binary_key_[j]);
            ICS_in->c_mineral_[j] = bkey_j;
            if (j < ICS_in->size_sup_min_) ICS_in->c_sup_min_[j] = bkey_j;
            else ICS_in->c_buffer_[j - ICS_in->size_sup_min_] =  bkey_j;
        }

        // Ditto: We only need to know whether the values are positive or zero:
        const int key = create_unique_key_from_positive_elements(workspace.global_binary_key_.data(),
                                                                 max_number_of_mineral_combinations);

        Vchem_in[key] = BasVec::createFromInitChem(ICS_in);

        if (ICS_in->PRINT_DEBUG_CHEM_) Vchem_in[key]->write("Vchem_key_" + std::to_string(key) + ".out");
    }
}