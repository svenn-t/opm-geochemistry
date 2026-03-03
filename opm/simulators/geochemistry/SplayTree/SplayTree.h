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
#ifndef SPLAY_TREE_H
#define SPLAY_TREE_H

#include <iostream>
#include <memory>

#include <opm/simulators/geochemistry/SplayTree/TreeNode.h>
#include <opm/simulators/geochemistry/Common/ChemGlobal.h>
#include <opm/simulators/geochemistry/Core/ChemBasVec.h>
#include <opm/simulators/geochemistry/Utility/MiscUtilityFunctions.hpp>

class SplayTree
{
public:
    SplayTree(int mineral_key, int num_key, int num_values);
    ~SplayTree();

    TreeNode* root();

    int mineralCombinationKey() const;

    std::size_t numberOfTreeNodes() const;
    std::size_t numberOfLookups() const;

    bool SplayNewKey(std::vector<int> key);
    void insertValues(int new_mineral_combination_key, double* new_values);

private:
    int mineral_combination_key_{-1};
    std::size_t num_success_lookup_{0};

    int num_key_{0};
    std::vector<int> new_key_;

    int num_values_{0};
    TreeNode* new_node_{nullptr};

    std::size_t numberOfTreeNodes_{0};
    TreeNode* root_{nullptr};
    TreeNode* nil_node_{nullptr};
    TreeNode* header_{nullptr};

 private:

    void rotateLeftChildOfRoot();
    void rotateRightChildOfRoot();
    void destroyTreeNode(TreeNode* node);

};

using SplayTreeUniquePtr = std::unique_ptr<SplayTree>;

#endif
