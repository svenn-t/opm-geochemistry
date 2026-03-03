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
#ifndef THERMO_DATA_H
#define THERMO_DATA_H

#include <array>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iterator>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include <opm/simulators/geochemistry/LinearAlgebra/LinearSolverRoutines.hpp>
#include <opm/simulators/geochemistry/Core/ChemTable.h>
#include <opm/simulators/geochemistry/Utility/HelpfulStringMethods.hpp>

class thnode
{

  public:

    int read_data(std::ifstream& dbf, int N, std::set<std::string>& basis_set_);

    static std::string remove_charge(std::string);

  public:

    bool in_list_{false};
    bool in_basis_{false};

    std::string name_{};
    std::string nick_name_{};
    std::string original_nick_name_{};
    std::string type_{};
    std::string stoi_{};
    std::string stoi_b_{};
    std::string ref_{};
    std::string date_{};

    std::vector<double> a_;
    std::vector<double> c_;
    std::vector<std::string> basis_;
    std::vector<double> dstoi_;
    std::vector<double> dstoi_full_;
    std::vector<double> dstoi_full_new_;

    double charge_{0.0};
    double charge_res_{0.0};  // Difference in charge between basis species and node-species.
    double G_{0.0};
    double S_{0.0};
    double H_{0.0};
    double Mw_{0.0};
};

class thermodata
{
  public:
    thermodata();
    ~thermodata();

    void create_files_from_slop07_ah();

  private:

    void supcrt2tables();
    void updateNode(std::vector<thnode*>& node, std::ifstream& dbf, int N);
    bool read_new_basis();

    void sub_basis_set();
    void abs_basis_set(std::vector<thnode*>& nodes);

    static thnode* get_thnode(const std::string& name, std::vector<thnode*>& nodes);

    void write_basis(const std::string& file_name, std::vector<thnode*>& nodes);
    void write_basis_new(const std::string& file_name, std::vector<thnode*>& nodes);
    void write_mineral(const std::string& file_name, std::vector<thnode*>& nodes);

    void db2cpp_basis(const std::string& data_file="tmp/supcrt_basis_basis.dat");
    void db2cpp_aq(const std::string& data_file="tmp/supcrt_basis_aq_complexes.dat");
    void db2cpp_minerals(const std::string& data_file_mineral="tmp/supcrt_min.dat",
                                const std::string& data_file_gas="tmp/supcrt_gas.dat", int skip_gas_header=1);

    std::size_t pos_electron() const;


    void extract_species_tokens_and_write_to_stream(std::ostream& os, const std::string& buffer, char delim='\t');
    int find_specie(const std::string name, std::vector<thnode*>& nodes);
    void remove_specie(const std::string name, std::vector<thnode*>& nodes);

  private:

    static const inline std::string separator_{ "\t" };
    static constexpr int output_precision_{8};

    int OXIDATION_STATES_{1};  // 1=Keep, 0=Skip.

    std::set<std::string> full_basis_set_;
    std::set<std::string> sub_basis_set_;
    std::set<std::string> new_basis_set_;

    double** Basis_Trans_{nullptr};
    double** Basis_Trans_inv_{nullptr};

    std::vector<double> new_basis_set_charge_;
    std::vector<double> new_basis_set_Mw_;
    std::vector<thnode*> min_;
    std::vector<thnode*> min1_;
    std::vector<thnode*> min2_;
    std::vector<thnode*> min3_;
    std::vector<thnode*> aq_;
    std::vector<thnode*> gas_;

};

#endif
