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
#include <opm/simulators/geochemistry/SplayTree/SplayTree.h>

#include <algorithm>

SplayTree::SplayTree(int mineral_key, int num_key, int num_values)
: mineral_combination_key_(mineral_key)
, num_key_(num_key)
, num_values_(num_values)
{
    new_key_.resize(num_key_);

    nil_node_ = new TreeNode(num_key_);
    nil_node_->left_ = nil_node_;
    nil_node_->right_ = nil_node_;
    header_ = new TreeNode(); // No need to store any values, only children are set during splaying.
    root_ = nil_node_;

    new_node_ = new TreeNode(num_key_);
    numberOfTreeNodes_ = 1;
}

SplayTree::~SplayTree()
{
    destroyTreeNode(root_);
    if(nil_node_)   delete nil_node_;
    if(header_)     delete header_;
    if(new_node_)   delete new_node_;
}

TreeNode* SplayTree::root()
{
    return root_;
}

int SplayTree::mineralCombinationKey() const
{
    return mineral_combination_key_;
}

std::size_t SplayTree::numberOfTreeNodes() const
{
    return numberOfTreeNodes_;
}

std::size_t SplayTree::numberOfLookups() const
{
    return num_success_lookup_;
}

/** 
* @param [in] key: New key that was generated prior to calling this function.
* @return: True iff the input key did not already exist in the splay tree,
*          meaning that a new solution has to be computed.
*/
bool SplayTree::SplayNewKey(std::vector<int> key)
{
    new_key_ = std::move(key);

    copy_from_vector_to_vector(new_key_, new_node_->key_);  // Copy key.

    if (root_ == nil_node_)
    {
        new_node_->left_ = nil_node_;
        new_node_->right_ = nil_node_;
        root_ = new_node_;

        return true;  // The tree was empty, so the key is certainly new.
    }

    // Otherwise, if tree is non-empty, continue...
    ++num_success_lookup_;

    copy_from_vector_to_vector(new_key_, nil_node_->key_);  // Copy key.

    header_->left_ = nil_node_;
    header_->right_ = nil_node_;

    TreeNode* left_tree_max = header_;
    TreeNode* right_tree_min = header_;
    while (compare_vectors(new_key_, root_->key_) != Comparison::EQUAL)
    {
        if (compare_vectors(new_key_, root_->key_) == Comparison::LEFT_SMALLER)
        {
            if (compare_vectors(new_key_, root_->left_->key_) == Comparison::LEFT_SMALLER)
            {
                rotateLeftChildOfRoot();
            }

            if (root_->left_ == nil_node_) break;

            // Link right
            right_tree_min->left_ = root_;
            right_tree_min = root_;
            root_ = root_->left_;
        }
        else
        {
            if (compare_vectors(new_key_, root_->right_->key_) == Comparison::LEFT_LARGER)
            {
                rotateRightChildOfRoot();
            }

            if (root_->right_ == nil_node_) break;

            // Link left
            left_tree_max->right_ = root_;
            left_tree_max = root_;
            root_ = root_->right_;
        }
    }

    // Reassemble
    //  - Right child of header_ points to left tree.
    //  - Left child of header_ points to right tree.
    left_tree_max->right_ = root_->left_;
    right_tree_min->left_ = root_->right_;
    root_->left_ = header_->right_;
    root_->right_ = header_->left_;

    if (compare_vectors(new_key_, root_->key_) == Comparison::LEFT_SMALLER)
    {
        new_node_->left_ = root_->left_;
        new_node_->right_ = root_;
        root_->left_ = nil_node_;
        root_ = new_node_;
    }
    else if (compare_vectors(new_key_, root_->key_) == Comparison::LEFT_LARGER)
    {
        new_node_->right_ = root_->right_;
        new_node_->left_ = root_;
        root_->right_ = nil_node_;
        root_ = new_node_;
    }

    return (root_ == new_node_);  // If false, it means the key already existed in the tree.
}

/** Must only be called when the root is the same as new_node. */
void SplayTree::insertValues(int new_mineral_combination_key, double* new_values)
{
    assert(root_ == new_node_);
    assert(num_values_ > 0);

    root_->val_.resize(num_values_); // Note: This is needed.

    if (new_mineral_combination_key == mineral_combination_key_)
    {   // The only situation in which this condition can fail to hold
        // is when there is a new super-saturated mineral. If the solver
        // diverges, it simply returns the old (i.e., input) key.
        for(int i=0; i < num_values_; ++i)
        {
            root_->val_[i] = new_values[i];
        }
    }

    // TODO: Here, the splay tree and its root can end up having different
    //       mineral keys, due to the appearance of a new mineral... Ok?
    root_->mineral_combination_key_ = new_mineral_combination_key;

    new_node_ = new TreeNode(num_key_);
    ++numberOfTreeNodes_;
}

void SplayTree::rotateLeftChildOfRoot()
{
    TreeNode* tmp = root_;
    root_ = tmp->left_;
    tmp->left_ = root_->right_;
    root_->right_ = tmp;
}

void SplayTree::rotateRightChildOfRoot()
{
    TreeNode* tmp = root_;
    root_ = tmp->right_;
    tmp->right_ = root_->left_;
    root_->left_ = tmp;
}

void SplayTree::destroyTreeNode(TreeNode* node)
{
    if(!node)   return;
    else if(node == nil_node_)
    {
        // We cannot delete the nil_node_ just yet...
        return;
    }
    else if (node == header_ || node == new_node_)
    {
        return;  // We should never enter here...
    }
    destroyTreeNode(node->left_);
    destroyTreeNode(node->right_);
    delete node;
}
