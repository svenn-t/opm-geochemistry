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
#ifndef CHEM_TABLE_H
#define CHEM_TABLE_H

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <list>
#include <map>
#include <numeric>
#include <random>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>
#include <iomanip>

#include <opm/simulators/geochemistry/Common/ChemGlobal.h>
#include <opm/simulators/geochemistry/Common/Enums.hpp>
#include <opm/simulators/geochemistry/Database/ChemDatabaseHKF.hpp>
#include <opm/simulators/geochemistry/IO/ErrorMsg.hpp>
#include <opm/simulators/geochemistry/IO/fmt_include.hpp>
#include <opm/simulators/geochemistry/IO/FileFunctions.hpp>
#include <opm/simulators/geochemistry/IO/ParseChemistry.hpp>
#include <opm/simulators/geochemistry/Utility/HelperMacros.hpp>
#include <opm/simulators/geochemistry/Utility/HelpfulStringMethods.hpp>
#include <opm/simulators/geochemistry/Utility/MiscUtilityFunctions.hpp>

class ChemTable
{

    MAKE_NONCOPYABLE(ChemTable)
        MAKE_NONMOVABLE(ChemTable)

public:

    ChemTable() = default;

    enum class Type { Basis, Complex, Mineral };

    // Used to copy ChemTables
    static constexpr int skip_col_bas() { return 0; };
    static constexpr int skip_col_complex() { return 14; };
    static constexpr int skip_col_min() { return 8; };

    /**
     * Shaves a table of elements that are not needed in the calculations, based on an input row vector.
     * If one of the existing row vectors has a 1 at a column index where the input vector has a value 0, it is to be taken out of the calculation.
     *
     * The index variables start_row and start_col are also required as input, to allow for the possibility to store more information in the M_ matrix.
     *
     * EXAMPLE: If start_row = 1,  start_col=1 and all the M_ values are to be considered:
     *
     *                  <col_name>       <col_name>           <col_name>          <col_name>           <col_name>             <col_name>
     *      <row_name>              M00                 M01                         M02                      M03                        M04                         M05
     *      <row_name>              M10                 M11                         M12                      M13                        M14                         M15
     *      etc ...
     *
     * */
    std::unique_ptr<ChemTable> create_reduced_table(const std::vector<int>& columns_to_keep, int start_row = 1, int start_col = 1) const;

    std::unique_ptr<ChemTable> create_reduced_table_rc(const std::vector<int>& rows_to_keep, const std::vector<int>& columns_to_keep, int start_row = 1, int start_col = 1) const;

    void read_csv(const std::string file_name, const std::string sep = "\t", bool nick_name = true);
    void write_csv(const std::string file_name, const std::string sep = "\t");
    void add_col_by_name(const std::string col_name, const std::string file_name, const std::string sep = "\t");
    void split_and_add_to_vector(std::vector<std::string>& elements, const std::string& buffer, int skip = 0, const std::string sep = "\t");
    void write(const std::string& file_name) const;

    void drop_col(std::vector<int>& columns_to_keep);

    void append(const ChemTable& CM);
    void sort_alphabetically(bool reverse = false);
    void rename_col(std::string old_name, std::string new_name);
    void rename_col(const std::vector<std::string>& old_names, const std::vector<std::string>& new_names);
    void swap_col(std::string col_name, int new_pos);
    void move_col(std::string col_name, int new_pos);

    int get_row_index(const std::string& name) const;
    int get_col_index(const std::string& name) const;
    bool name_is_in_row(const std::string& name, int row_index) const;

    std::vector<int> find_row_elements(const std::vector<std::string>& names) const;

    int check_for_duplicate_species(bool is_case_sensitive = true, bool check_nick_names = false);
    bool table_has_duplicate_row_names(bool case_sensitive = true) const;
    bool table_has_duplicate_nick_names(bool case_sensitive = true) const;

    void initialize_full_database_hkf(ChemTable::Type type, const std::vector < std::pair<std::string, std::string>>& replace_names = {});
    void init_tables();
    void add_cubic_EOS_parameters_for_gases();

    void remove_row_species(const std::vector<std::string>& species_to_remove);

    int add_col(const std::vector<std::string>& names, std::vector<double> default_values = {});
    int add_species_row(const std::string& name, const std::vector<double>& input_values);

    int add_reaction(const std::string& reac,
        int type,
        LogKModel model,
        double charge,
        double a0,
        double mol_volume,
        const std::vector<double>& hkf,
        const std::vector<std::string>& basis_specie_names,
        const std::vector<double>& stoic);

    int find_no_basis_species(int pos_buf, int* pos_bas, const std::vector<int>& ignore_list);

    void in_fmatrix(const int* row, int dim_row, const int* col, int dim_col, double** f_red);

    void resize(int no_rows, int no_columns, double default_col_val = 0.);

    void transpose();
    void update_number_of_analytical_species();
    void update_logK_analytical(double T,double P=1e5);
    void fix_hkf_units(ChemTable::Type type);
    void set_up_sparse_matrix();

    std::string get_reaction(const int i, int prec = 0) const;

    bool set_value(const std::string name_of_vector, const int index, const double new_value);
    bool get_value(const std::string name_of_vector, const int index, double& return_val);

    bool set_value(const std::string name_of_vector, const int index, const std::string new_value);
    bool get_value(const std::string name_of_vector, const int index, std::string & return_val);

public:

    double T_curr_{0.0};  // Current temperature
    double P_curr_{ 0.0 };  // Current pressure

    // Dimensions of matrix
    int noRows_{ 0 };
    int noColumns_{ 0 };

    std::vector<std::string> database_;

    std::vector<std::string> col_name_;
    std::vector<std::string> row_name_;   // Complexes or rock buffers
    std::vector<std::string> nick_name_;  // Complexes or rock buffers

    std::vector<double> T_;             // Temperature row vector

    std::vector<LogKModel> model_;      // Note: HKF first, then analytical
    int size_analytical_{ 0 };

    std::vector<double> log_m_;         // log concentrations
    std::vector<double> log_a_;         // log activities
    std::vector<double> log_g_;         // log activity coefficients

    std::vector<double> delta_;         // Amounts of precipitated minerals
    std::vector<double> fugacity_;      // Fugacities
    std::vector<double> log_QK_;        // log_Q/K value for mineral, usually 0

    std::vector<int> type_;             // Type of geochemical component
    std::vector<double> logK_;          // logK values for complexes / minerals
    std::vector<double> log_af_;        // log activity const-buffers
    std::vector<double> a0_;            // Extended Debye-Hückel constant
    std::vector<double> charge_;        // Charge of aqueous complex
    std::vector<double> scharge_;       // Charge of surface complex
    std::vector<double> mol_volume_;    // == molecular weight/density
    std::vector<double> mol_weight_;    // Molecular weight

    // HKF EOS parameters that are read directly from the input database
    std::vector<double> deltaG_;
    std::vector<double> deltaH_;
    std::vector<double> S_;
    std::vector<double> a1_;
    std::vector<double> a2_;
    std::vector<double> a3_;
    std::vector<double> a4_;
    std::vector<double> a7_;
    std::vector<double> c1_;
    std::vector<double> c2_;
    std::vector<double> omega_;

    // HKF EOS parameters that are calculated later
    std::vector<double> dG_TP_;
    std::vector<double> re_;  // Effective ion radius @ P, T
    std::vector<double> molar_volume_PT_;
    std::vector<double> Cp_PT_;

    // Peng-Robinson EOS State Parameter
    // Currently only used for Gases 
    std::vector<double> Tc_;
    std::vector<double> Pc_;
    std::vector<double> omega_acc_;


    // Stoichiometric matrix, both dense and sparse versions
    std::vector<std::vector<double> > M_;
    std::vector<std::vector<std::pair<int, double>>> sparseM_;

    // Transposed stoichiometric matrix
    std::vector<std::vector<double> > M_trans_;

    std::vector<int> abs_pos_;      // Absolute position in full database
    std::vector<int> abs_pos_bas_;  // Absolute position in full basis vector

    std::map<std::string, std::vector<double>&> double_vector_map_=
    {
        {"DeltaG", deltaG_},
        {"DeltaH", deltaH_},
        {"S", S_ },
        {"a1", a1_ },
        {"a2", a2_ },
        {"a3", a3_ },
        {"a4", a4_ },
        {"a7", a7_ },
        {"c1", c1_ },
        {"c2", c2_ },
        {"logK", logK_},
        {"log_af", log_af_ },
        {"a0", a0_ },
        {"charge", charge_ },
        {"scharge", scharge_},
        {"mol_volume", mol_volume_ },
        {"mol_weight", mol_weight_ },
        {"omega", omega_ },
        {"log_m", log_m_ },
        {"log_a", log_a_ },
        {"log_g", log_g_ },
        {"delta", delta_ },
        {"fugacity", fugacity_ },
        {"log_QK", log_QK_ },
        {"T", T_ },
        {"Tc", Tc_ },
        {"Pc", Pc_ },
        {"omega_acc", omega_acc_ }
    };
    std::map<std::string, std::vector<std::string>&> string_vector_map_ =
    { {"row_name", row_name_},
      {"nick_name", nick_name_},
    };

private:
    static constexpr int output_precision_{8};
    enum class CopyMode { Column = 0, Row = 1 };

    template<typename T>
    static void copy_ChemTable_vector(std::vector<T>& vec_out,
                                      const ChemTable& CM,
                                      const std::vector<T>& vec,
                                      CopyMode mode)
    {
        if(vec.empty()) return;

        const auto& pos = (mode == CopyMode::Column) ? CM.abs_pos_ : CM.abs_pos_bas_;
        const int vec_size = (mode == CopyMode::Column) ? CM.noRows_ : CM.noColumns_;

        vec_out.resize(vec_size);
        for(std::size_t i=0; i < vec_out.size(); ++i)
        {
            vec_out[i] = vec[pos[i]];
        }
    }

    template<typename T>
    static void copy_ChemTable_column(std::vector<T>& vec_out,
                                      const ChemTable& CM,
                                      const std::vector<T>& vec)
    {
        copy_ChemTable_vector(vec_out, CM, vec, CopyMode::Column);
    }

    template<typename T>
    static void copy_ChemTable_row(std::vector<T>& vec_out,
                                   const ChemTable& CM,
                                   const std::vector<T>& vec)
    {
        copy_ChemTable_vector(vec_out, CM, vec, CopyMode::Row);
    }

    /* Returns column vector if exists, and a zero vector if not found.
     * Can be used for both columns of ints and columns of doubles. */
    template<typename T, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
    std::vector<T> pick_column(const std::string& name)
    {
        std::vector<T> col;
        col.resize(noRows_, static_cast<T>(0));

        for(int j=0; j < noColumns_; ++j)
        {
            if(col_name_[j] == name)
            {
                for(int i=0; i < noRows_; ++i)
                {
                    col[i] = static_cast<T>(M_[i][j]);
                }
                return col;
            }
        }
        return col;
    }
    
private:
    
    void erase_row(std::size_t pos);
    void erase_row_elements(const std::vector<std::string>& names);
    
    void resize_col(int size, double default_val=0.);
    void resize_row(int size);

    void copy_data_to_new_table(ChemTable& CM) const;

    /* Helper function to get a string-representation of the Type enum. */
    std::string typeToString(ChemTable::Type type) const
    {
        switch (type)
        {
            case Type::Basis: return "Basis";
            case Type::Complex: return "Complex";
            case Type::Mineral: return "Mineral/gas buffer";
            default: return "Undefined type";
        }
    }
    
};

#endif
