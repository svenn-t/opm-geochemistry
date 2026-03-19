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
#include <opm/simulators/geochemistry/Thermo/thermodata.h>

#include <filesystem>
#include <iomanip>

#include <opm/simulators/geochemistry/IO/ErrorMsg.hpp>
#include <opm/simulators/geochemistry/Utility/HelpfulStringMethods.hpp>

int thnode::read_data(std::ifstream& dbf, int N,  std::set<std::string>& basis_set)
{
    std::istringstream iss;

    std::string line;
    std::getline(dbf, line);
    iss.str(line);

    std::string dname;
    if( !(iss >> dname >> stoi_) )
    {
        return 1;
    }

    // Bugfix (31/12-2022):
    //  Before, we searched for "AQ" by finding the index of the first
    //  comma. However, this method will NOT work in general, because the
    //  species name itself may contain comma(s).
    //
    //  Example: "1,1-DCE,AQ"
    auto remove_AQ_identifier = [](const std::string& str) -> std::string
    {
        static constexpr std::array<const char*, 2> search_strings = {
            ",AQ",
            ", AQ",
        };

        for(auto& str_search: search_strings)
        {
            std::size_t found = str.find(str_search);
            if(found != std::string::npos)
            {   // We assume that at most one of the search strings are present.
                return str.substr(0, found);
            }
        }
        return str;
    };

    name_ = remove_AQ_identifier(dname);  // Nick name can still have (AQ)

    iss.clear();
    std::getline(dbf, line);
    iss.str(line);

    // Note: We currently use the original nick name for secondary species.
    iss >> original_nick_name_ >> stoi_b_;

    // Hack: We keep the charge in the nick name for the electron, to be able
    //       to distinguish it from the upper case version of the surface
    //       complex specie E.
    if(name_ == "e-")   nick_name_ = original_nick_name_;
    else                nick_name_ = remove_charge(original_nick_name_);

    // At least for now, change (g) to (G), since that is what we are currently
    // doing in ChemDataBaseHKF_buffer.cpp.
    std::size_t found_lower_case_gas = nick_name_.find("(g)");
    if(found_lower_case_gas != std::string::npos)
    {
        nick_name_ = nick_name_.substr(0, found_lower_case_gas) + "(G)";
    }

    std::getline(dbf, line);
    iss.clear();
    iss.str(line);
    iss >> ref_ >> date_;

    for (int i = 0; i < N; ++i)
    {
        iss.clear();
        std::getline(dbf, line);
        iss.str(line);
        double a;
        while (iss >> a)
        {
            a_.push_back(a);
        }
    }
    iss.clear();

    std::string bn;
    for(const char stoi_b: stoi_b_)
    {
        if (stoi_b == '(' || stoi_b == ')')         bn += ' ';
        else                                        bn += stoi_b;
    }

    iss.clear();
    iss.str(bn);
    std::string basis;

    charge_ = 0.0;

    double d;
    while (iss >> basis >> d)
    {
        if (basis == "-")
        {
            charge_ = -d;
        }
        else if(basis == "+")
        {
            charge_ = d;
        }
        else
        {
            basis_.push_back(basis);
            dstoi_.push_back(d);
            basis_set.insert(basis);
        }

    }
    return 0;
}

std::string thnode::remove_charge(std::string full_name)
{
    // Remove charges from the full_name
    // Examples: Ca+2 -> Ca, Cl- ->Cl.
    std::size_t found = full_name.find_first_of("+-");
    if (found != std::string::npos)
    {
        if ((found + 2) == full_name.size() || (found+1)==full_name.size())
        {
            return full_name.substr(0, found);
        }
    }
    return full_name;
}

/******************************************************************/
/* Read supcrt database and write standard tables                 */
/* File from http://geopig.asu.edu/sites/default/files/slop07.dat */
/* File has been modified:                                        */
/* Line: 1694 added H2O as a dummy species                        */
/* Line: 1650 o-->O                                               */
/* Line: 765, 2235, 2253, 2259, 4329, 4539, 4545, 4653 fixed basis*/
/* Line: 3219 h-->H                                               */
/* BeOH+ wrong charge -2.-->1                                     */
/* Line:  U(Pent)+2           removed space                       */
/******************************************************************/
thermodata::thermodata() = default;

thermodata::~thermodata()
{
    const auto sizeB = static_cast<int>(new_basis_set_.size());
    if(Basis_Trans_)        freeMatrixMemory(Basis_Trans_, sizeB);
    if (Basis_Trans_inv_)   freeMatrixMemory(Basis_Trans_inv_, sizeB);
}

void thermodata::create_files_from_slop07_ah()
{
    supcrt2tables();
    read_new_basis();
    sub_basis_set();

    // TODO: Update these three functions so that the generated text files
    //       are consistent with the source code in GeoChemicalModel/Database.
    db2cpp_basis();
    db2cpp_aq();
    db2cpp_minerals(/*data_file_mineral =*/ "tmp/supcrt_min123.dat");
}

void thermodata::supcrt2tables()
{
    int imin = 0;
    int imin1 = 0;
    int imin2 = 0;
    int imin3 = 0;
    int igas = 0;
    int iaq = 0;

    static const std::string db_file = "slop07_ah.dat";
    if (!std::filesystem::is_regular_file(db_file))
    {
        error_and_exit("thermodata::supcrt2tables() requires the file\"{:s}\" to exist in the current working directory...", db_file);
    }

    std::ifstream dbf(db_file, std::ios::binary);

    if (dbf.is_open())
    {
        std::string buf;
        while (getline(dbf, buf))
        {
            if(buf.find('*') != std::string::npos); // Do nothing
            else if(string_contains(buf, "minerals that do not undergo phase transition") && !imin)
            {
                updateNode(min_, dbf, 3);
                ++imin;
            }
            else if(string_contains(buf, "minerals that undergo one phase transition") && !imin1)
            {
                updateNode(min1_, dbf, 4);
                ++imin1;
            }
            else if(string_contains(buf, "minerals that undergo two phase transitions") && !imin2)
            {
                updateNode(min2_, dbf, 5);
                ++imin2;
            }
            else if(string_contains(buf, "minerals that undergo three phase transitions") && !imin3)
            {
                updateNode(min3_, dbf, 6);
                ++imin3;
            }
            else if(string_contains(buf, "gases") && !igas)
            {
                updateNode(gas_, dbf, 3);
                ++igas;
            }
            else if(string_contains(buf, "aqueous species") && !iaq)
            {
                updateNode(aq_, dbf, 3);
                ++iaq;
            }
        }
    }
    else
    {
        std::cerr << "Could not open database file \"" << db_file << "\"\n";
    }
}

void thermodata::updateNode(std::vector<thnode*>& nodes, std::ifstream& dbf, int N)
{
    std::string buf;
    getline(dbf, buf);

    bool finished = false;
    do
    {
        thnode* node_dum = new thnode;
        const int err_code = node_dum->read_data(dbf, N, full_basis_set_);
        if(err_code > 0)
        {
            finished = true;
            delete node_dum;
        }
        else
        {
            nodes.push_back(node_dum);
        }
    } while (!finished);
}

/**
 *  We read a new basis set and express the old basis set in terms of the new, e.g.:
 *
 *      HCO3- = H + C + 3*O  (etc.)
 *
 * The charge is balanced by introducing aqueous electrons: pe = -log a_e.
 *
 * Returns True if successful. */
bool thermodata::read_new_basis()
{
    static const std::string basis_species_file = "basis_species.dat";
    if (!std::filesystem::is_regular_file(basis_species_file))
    {
        error_and_exit("thermodata::read_new_basis() requires the file\"{:s}\" to exist in the current working directory...", basis_species_file);
    }

    std::ifstream ibf(basis_species_file, std::ios::binary);
    if(!ibf.is_open())
    {
        std::cout << "Error: Could not read database file containing basis species!\n";
        return false;
    }

    // TODO: Why do we need these? B/c species can be repeated in the input?
    std::vector<std::string> basis_specie_names;
    std::vector<double> mw;

    std::string sbn;
    double mwi;
    while (ibf >> sbn >> mwi)
    {
        auto ret_ = new_basis_set_.insert(sbn);
        if (ret_.second)
        {   // New basis specie was inserted
            mw.push_back(mwi);
            basis_specie_names.push_back(sbn);
        }
    }

    new_basis_set_Mw_.resize(new_basis_set_.size());

    std::size_t i = 0;
    for(const auto& basis_specie: basis_specie_names)
    {
        // TO DO: Add explanatory note here?
        auto sit = new_basis_set_.find((basis_specie));
        auto pos = static_cast<std::size_t>(std::distance(new_basis_set_.begin(), sit));
        new_basis_set_Mw_[pos] = mw[i];
        ++i;
    }
  //  sub_basis_set();

    return true;
}

/* Fills in stoichiometric matrix in full data set, both for old basis and new basis. */
void thermodata::abs_basis_set(std::vector<thnode*>& nodes)
{
    std::size_t pos_e = pos_electron();

    for(auto* node: nodes)
    {
        bool found = true;
        for(const auto& basis_specie: node->basis_)
        {
            auto sit = sub_basis_set_.find(basis_specie);
            if (sit == sub_basis_set_.end())
            {
                found = false;
                node->in_list_ = false;
                node->in_basis_ = false;
            }
        }

        if (found)
        {
            node->in_list_ = true;
            auto sit = new_basis_set_.find(node->name_);
            if (sit == new_basis_set_.end())
            {
                node->in_basis_ = false;
            }
            else
            {
                //std::cout << node->name_ << "\n";
                node->in_basis_ = true;
            }

            node->dstoi_full_.resize(new_basis_set_.size(), 0.0);

            std::size_t basis_count = 0;
            for(const auto& basis_specie: node->basis_)
            {
                sit = sub_basis_set_.find(basis_specie);
                auto pos_i = static_cast<std::size_t>(std::distance(sub_basis_set_.begin(), sit));
                node->dstoi_full_[pos_i] = node->dstoi_[basis_count];
                ++basis_count;
            }

            node->dstoi_full_new_.resize(node->dstoi_full_.size());

            for (std::size_t i = 0; i < node->dstoi_full_.size(); ++i)
            {
                node->dstoi_full_new_[i] = 0.0;
                for (std::size_t j = 0; j < node->dstoi_full_.size(); ++j)
                {
                    node->dstoi_full_new_[i] += Basis_Trans_inv_[j][i] * node->dstoi_full_[j];
                }
            }

            // Value of "e-" is determined by requiring charge balance
            node->charge_res_ = node->charge_;
            for (std::size_t i = 0; i < node->dstoi_full_new_.size(); ++i)
            {
                node->charge_res_ -= node->dstoi_full_new_[i] * new_basis_set_charge_[i];
            }
            node->dstoi_full_new_[pos_e] = -node->charge_res_; //@ah wrong sign of e-
        }
    }
}

void thermodata::sub_basis_set()
{
    for (const auto& basis_specie : new_basis_set_)
    {
        thnode* th = get_thnode(basis_specie, aq_);
        if (!th)
        {
            std::cout << "No basis named " << basis_specie << " in database\n";
        }
        else
        {
            for (const auto& specie : th->basis_)
            {
                sub_basis_set_.insert(specie);
            }
        }
    }
    if (sub_basis_set_.size() != new_basis_set_.size())
    {
        error_and_exit("Wrong dimension for basis set!");
    }

    const auto sizeB = static_cast<int>(new_basis_set_.size());

    // Basis transformation matrix and inverse
    Basis_Trans_ = allocateMemoryForMatrix(sizeB);
    Basis_Trans_inv_ = allocateMemoryForMatrix(sizeB);
    fillMatrix(Basis_Trans_, sizeB, 0.0);
    fillMatrix(Basis_Trans_inv_, sizeB, 0.0);

    std::size_t pos_i_exist = 0;
    std::size_t pos_i = 0;
    for (const auto& basis_specie : new_basis_set_)
    {
        thnode* th = get_thnode(basis_specie, aq_);
        if (th == nullptr)
        {
            std::cout << "No basis named " << basis_specie << "in database.\n";
        }
        else
        {
            new_basis_set_charge_.push_back(th->charge_);
            th->Mw_ = new_basis_set_Mw_[pos_i_exist];
            ++pos_i_exist;

            std::size_t pos_th = 0;
            for (const auto& specie : th->basis_)
            {
                auto it = sub_basis_set_.find(specie);
                auto pos_j = static_cast<std::size_t>(std::distance(sub_basis_set_.begin(), it));
                Basis_Trans_[pos_i][pos_j] = th->dstoi_[pos_th];
                ++pos_th;
            }
        }
        ++pos_i;
    }

    double determinant;
    invert_matrix(Basis_Trans_,
        sizeB,
        Basis_Trans_inv_,
        determinant);
    remove_specie("BaCO3", aq_); // HKF contains both BaCO3 and Ba(CO3)
    remove_specie("CaCO3", aq_); // HKF contains both CaCO3 and Ca(CO3)
    remove_specie("MgCO3", aq_); // HKF contains both MgCO3 and Mg(CO3)
    remove_specie("SrCO3", aq_); // HKF contains both SrCO3 and Sr(CO3)
    abs_basis_set(aq_);
    abs_basis_set(min_);
    abs_basis_set(min1_);
    abs_basis_set(min2_);
    abs_basis_set(min3_);
    abs_basis_set(gas_);

    if (!std::filesystem::is_directory("tmp"))
    {
        error_and_exit("No directory \"tmp\" exists within the current working directory...");
    }

    write_basis_new("tmp/supcrt_basis", aq_);
    write_basis("tmp/aq_supcrt_old", aq_);
    write_mineral("tmp/supcrt_min", min_);
    write_mineral("tmp/supcrt_min1", min1_);
    write_mineral("tmp/supcrt_min2", min2_);
    write_mineral("tmp/supcrt_min3", min3_);
    write_mineral("tmp/supcrt_gas", gas_);

    // add missing col
    ChemTable* CM, * CM1;
    CM = new ChemTable;
    CM->read_csv("tmp/supcrt_basis_basis.dat");
    std::vector<std::string> old_col_b = { "Mw","a0","a1","a2","a3","a4","a5","a6","a7","a8","a9","a10" };
    std::vector<std::string> new_col_b = { "mol_weight","DeltaG","DeltaH","S","a1","a2","a3","a4","c1","c2","omega","charge" };
    CM->rename_col(old_col_b, new_col_b);
    CM->add_col_by_name("a0", "a0_vals.dat");
    CM->add_col(std::vector<std::string> {"scharge", "type"}, std::vector<double>{0., 0.});
    std::vector<std::string> col_nb = { "a0","mol_weight","scharge", "charge", "type" };

    for (auto& nn : col_nb)
    {
        CM->move_col(nn, 0);
    }

    CM->write_csv("tmp/supcrt_basis_basis.dat");
    delete CM;

    CM = new ChemTable;
    CM->read_csv("tmp/supcrt_basis_aq_complexes.dat");
    std::vector<std::string> old_col_aq = { "a0","a1","a2","a3","a4","a5","a6","a7","a8","a9","a10"};
    std::vector<std::string> new_col_aq = { "DeltaG","DeltaH","S","a1","a2","a3","a4","c1","c2","omega","charge"};
    CM->rename_col(old_col_aq, new_col_aq);
    std::vector<std::string> col_names = { "a0", "scharge","type"};
    std::vector<double> col_vals = {DEFAULT_ION_SIZE,0.,0.};
    CM->add_col(col_names, col_vals);
    col_names.insert(col_names.begin() + 2, "charge");
    for (auto& nn : col_names)
    {
        CM->move_col(nn, 0);
    }
    CM->write_csv("tmp/supcrt_basis_aq_complexes.dat");

    delete CM;

    // Need to combine min, min1, min2, min3 files - min1-min3 have more parameters that
    // are (currently) not used in the HKF EOS

    std::vector<std::string>fnames = { "tmp/supcrt_min1.dat","tmp/supcrt_min2.dat","tmp/supcrt_min3.dat" };

    CM = new ChemTable;
    CM->read_csv("tmp/supcrt_min.dat");
    for (auto it = fnames.begin(); it != fnames.end(); ++it)
    {
        CM1 = new ChemTable;
        CM1->read_csv(*it);
        CM->append(*CM1);
        delete CM1;
    }
    CM->move_col("a3", 0);
    std::vector<std::string> old_col = { "a3","a0","a1","a2","a4","a5","a6" };
    std::vector<std::string> new_col = { "mol_volume","DeltaG","DeltaH","S","a1","a2","a3" };
    CM->rename_col(old_col, new_col);
    CM->sort_alphabetically();
    CM->write_csv("tmp/supcrt_min123.dat");
    delete CM;

    // Do the same for gas:
    CM = new ChemTable;
    CM->read_csv("tmp/supcrt_gas.dat");
    CM->move_col("a3", 0);
    CM->rename_col(old_col, new_col);
    CM->write_csv("tmp/supcrt_gas.dat");
    delete CM;
    //----------------------------------------
    //
    // Debug files
    std::ofstream of_dbg("tmp/debug.dat", std::ios::out);
    for (const auto& bas : sub_basis_set_)
    {
        of_dbg << bas << "\n";
    }

    std::ofstream of_dbg2("tmp/debug2.dat", std::ios::out);
    for (const auto& bas : new_basis_set_)
    {
        of_dbg2 << bas << "\n";
    }

    std::ofstream of_dbg3("tmp/debug3.dat", std::ios::out);
    of_dbg3 << "NewBasis\t";
    for (const auto& bas : sub_basis_set_)
    {
        of_dbg3 << bas << "\t";
    }
    of_dbg3 << "\n";

    pos_i = 0;
    for(const auto& bas: new_basis_set_)
    {
        of_dbg3 << bas <<"\t";

        std::size_t pos_j = 0;
        for([[maybe_unused]] const auto& bas_sub: sub_basis_set_)
        {
            of_dbg3 << Basis_Trans_[pos_i][pos_j] << "\t";
            ++pos_j;
        }
        of_dbg3 << "\n";
        ++pos_i;
    }

    std::ofstream of_dbg4("tmp/debug4.dat", std::ios::out);
    of_dbg4 << "NewBasisInv\t";
    for (const auto& bas : new_basis_set_)
    {
        of_dbg4 << bas << "\t";
    }
    of_dbg4 << "\n";

    pos_i = 0;
    for(const auto& bas_sub: sub_basis_set_)
    {
        of_dbg4 << bas_sub << "\t";

        std::size_t pos_j = 0;
        for([[maybe_unused]] const auto& new_bas: new_basis_set_)
        {
            of_dbg4 << Basis_Trans_inv_[pos_i][pos_j] << "\t";
            ++pos_j;
        }
        of_dbg4 << "\n";
        ++pos_i;
    }
}

/** Searches for element using both name and nick_name. */
thnode* thermodata::get_thnode(const std::string& name, std::vector<thnode*>& nodes)
{
    for(auto* node: nodes)
    {
        if(node->name_ == name || node->nick_name_ == name) return node;
    }
    return nullptr;
}

void thermodata::write_basis_new(const std::string& file_name,
                                 std::vector<thnode*>& nodes)
{
    if(nodes.empty())
    {
        std::cout << "thermodata::write_basis_new(): There are no thnodes.\n";
        return;
    }

    std::ofstream of_basis(file_name + "_basis" + ".dat", std::ios::out);
    if(!of_basis.is_open())
    {
        std::cout << "thermodata::write_basis_new(): Cannot create file \"";
        std::cout << file_name << "\", have you checked that the directory exists?\n";
        return;
    }

    std::ofstream of_aq(file_name + "_aq_complexes" + ".dat", std::ios::out);

    of_basis << "Name\tNickName\tMw" << std::setprecision(output_precision_);
    of_aq << "Name\tNickName" << std::setprecision(output_precision_);
    for (std::size_t i = 0; i < nodes[0]->a_.size(); ++i)
    {
        of_basis << "\t" << "a" << i;
        of_aq << "\t" << "a" << i;
    }
    for (const auto& bas : new_basis_set_)
    {
        of_aq << "\t" << bas;
    }
    of_basis << "\n";
    of_aq << "\n";

    for(const auto& bas_specie: new_basis_set_)
    {
        for(const auto& node: nodes)
        {
            if (bas_specie == node->name_ || bas_specie == node->nick_name_)
            {
                of_basis << node->name_ << "\t" << node->nick_name_ << "\t" << node->Mw_;

                for (double a_i : node->a_)
                {
                    of_basis << "\t" << a_i;
                }
                of_basis << "\n";
                break;
            }
        }
    }

    const auto pos_e = pos_electron();

    for(const auto& node: nodes)
    {
        if(node->in_list_
           && !node->in_basis_
           // Skip oxidation states?
           && (node->dstoi_full_new_[pos_e] == 0 || OXIDATION_STATES_))
        {
            of_aq << node->name_ << "\t" << node->nick_name_;

            for (double a_i : node->a_)
            {
                of_aq << "\t" << a_i;
            }
            for (double dstoi : node->dstoi_full_new_)
            {
                of_aq << "\t" << dstoi;
            }
            of_aq << "\n";
        }
    }
}

void thermodata::write_mineral(const std::string& file_name, std::vector<thnode*>& nodes)
{
    std::ofstream of_min(file_name + ".dat", std::ios::out);

        of_min << "Name\tNickName" << std::setprecision(output_precision_);
        for (std::size_t i = 0; i < nodes[0]->a_.size(); ++i)
        {
            of_min << "\t" << "a" << i;
        }

    for (const auto& bas : new_basis_set_)
    {
        of_min << "\t" << bas;
    }
    of_min << "\n";

    auto pos_e = pos_electron();

    for(const auto* node: nodes)
    {
        if (node->in_list_ && (node->dstoi_full_new_[pos_e] == 0 || OXIDATION_STATES_))
        {
            of_min << node->name_ << "\t" << node->nick_name_;
            for (double a_i : node->a_)
            {
                of_min << "\t" << a_i;
            }
            for (double dstoi_i : node->dstoi_full_new_)
            {
                of_min << "\t" << dstoi_i;
            }
            of_min << "\n";
        }
    }
}

void thermodata::write_basis(const std::string& file_name, std::vector<thnode*>& nodes)
{
    if(nodes.empty())
    {
        std::cerr << "thermodata::write_basis: There are no thnodes.\n";
        return;
    }

    std::ofstream of_bas(file_name + ".dat", std::ios::out);

    of_bas << "Name\tNickName" << std::setprecision(output_precision_);
    for (std::size_t i = 0; i < nodes[0]->a_.size(); ++i)
    {
        of_bas << "\t" << "a" << i;
    }
    for (const auto& bas : sub_basis_set_)
    {
        of_bas << "\t" << bas;
    }
    of_bas << "\n";

    for(const auto* node: nodes)
    {
        of_bas << node->name_ << "\t" << node->nick_name_;
        for (double a_i : node->a_)
        {
            of_bas << "\t" << a_i;
        }
        for (double dstoi_i : node->dstoi_full_)
        {
            of_bas << "\t" << dstoi_i;
        }
        of_bas << "\n";
    }
    of_bas << "\n";
}

/** Basis species database. */
void thermodata::db2cpp_basis(const std::string& data_file)
{
    std::ifstream dbf(data_file, std::ios::in);
    if (!dbf.is_open())
    {
        error_and_exit("thermodata::db2cpp_basis(): Cannot open file \"{}\"...", data_file);
    }

    std::ofstream cpp_file("tmp_HKF_basis.cpp", std::ios::out);
    cpp_file << "#include \"ChemDatabaseHKF.hpp\"\n\n";
    cpp_file << "auto create_basis_species_table_hkf() -> std::vector<std::string>\n{\n";
    cpp_file << separator_ << "std::vector<std::string> basis_db;\n\n";

    std::string buf;
    while (std::getline(dbf, buf))
    {
        cpp_file << separator_ << "basis_db.push_back(\"";
        extract_species_tokens_and_write_to_stream(cpp_file, buf);
        cpp_file << "\");\n";
    }

    // TODO: At least for now, surface species are hard-coded.
    //
    // Inconsistency: Here, we use 4 whitespace characters to separate entries...(instead of tabs)
    //auto str_sc_io_basis = R"(    basis_db.push_back("X-  X   2   -1  0   0   0   0   0   0   0   0   0   0   0   0   0");
    //basis_db.push_back("GCO3-   GCO3    1   0   -1  0   0   -140282 -164898 23.53   7.5621  1.1505  1.2346  -2.8266 12.9395 -4.7579 1.2733");
    //basis_db.push_back("GCa+    GCa 1   0   1   0   0   -132120 -129800 -13.5   -0.1947 -7.252  5.2966  -2.4792 9   -2.522  1.2366");
    //basis_db.push_back("E   E   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0");
    //basis_db.push_back("GSiOH   GSiOH   1   0   0   0   0   -199190 -209775 18  1.9 1.7 20  -2.7    29.1    -51.2   0.1291");)";

    //cpp_file << "\n";
    //cpp_file << separator_ << "// Next, add surface species.\n";
    //cpp_file << separator_ << "// Previously used a separate function add_sc_io_basis_species().\n";
    //cpp_file << str_sc_io_basis;
    //cpp_file << "\n";

    cpp_file << "\n";
    cpp_file << separator_ << "return basis_db;\n";
    cpp_file << "}\n";
}

void thermodata::db2cpp_aq(const std::string& data_file)
{
    std::ifstream dbf_aq(data_file, std::ios::in);
    if (!dbf_aq.is_open())
    {
        error_and_exit("thermodata::db2cpp_aq(): Cannot open file \"{}\"...", data_file);
    }

    std::ofstream cpp_file("tmp_HKF_aq.cpp", std::ios::out);
    cpp_file << "#include \"ChemDatabaseHKF.hpp\"\n\n";
    cpp_file << "auto create_complex_species_table_hkf() -> std::vector<std::string>\n{\n";
    cpp_file << separator_ << "std::vector<std::string> aq_db;\n\n";

    std::string buf;
    while (std::getline(dbf_aq, buf))
    {
        cpp_file << separator_ << "aq_db.push_back(\"";  // << buf;
        extract_species_tokens_and_write_to_stream(cpp_file, buf);
        cpp_file << "\");\n";
    }

    cpp_file << "\n";
    cpp_file << separator_ << "return aq_db;\n";
    cpp_file << "};\n";
}

void thermodata::db2cpp_minerals(const std::string& data_file_mineral, const std::string& data_file_gas, int skip_gas_header)
{
    std::ifstream dbf_min(data_file_mineral, std::ios::in);
    if (!dbf_min.is_open())
    {
        error_and_exit("thermodata::db2cpp_minerals(): Cannot open file \"{}\"...", data_file_mineral);
    }
    std::ifstream dbf_gas(data_file_gas, std::ios::in);
    if (!dbf_gas.is_open())
    {
        error_and_exit("thermodata::db2cpp_minerals(): Cannot open file \"{}\"...", data_file_gas);
    }

    // Now we know the input files could be opened...
    std::ofstream cpp_file("tmp_HKF_buffer.cpp", std::ios::out);
    cpp_file << "#include \"ChemDatabaseHKF.hpp\"\n\n";
    cpp_file << "auto create_minerals_table_hkf() -> std::vector<std::string>\n{\n";
    cpp_file << separator_ << "std::vector<std::string> minerals_db;\n\n";

    std::string buf;
    while (std::getline(dbf_min, buf))
    {
        cpp_file << separator_ << "minerals_db.push_back(\"";  // << buf;
        extract_species_tokens_and_write_to_stream(cpp_file, buf);
        cpp_file << "\");\n";
    }
    int lines = 0;
    while (std::getline(dbf_gas, buf))
    {
        lines++;
        if (lines > skip_gas_header)
        {
            cpp_file << separator_ << "minerals_db.push_back(\"";  // << buf;
            extract_species_tokens_and_write_to_stream(cpp_file, buf);
            cpp_file << "\");\n";
        }
    }

    cpp_file << "\n";
    cpp_file << separator_ << "return minerals_db;\n";
    cpp_file << "};\n";
}

std::size_t thermodata::pos_electron() const
{
    auto it = new_basis_set_.find("e-");
    return static_cast<std::size_t>(std::distance(new_basis_set_.begin(), it));
}

void thermodata::extract_species_tokens_and_write_to_stream(std::ostream& os, const std::string& buffer, char delim)
{
    const auto tokens = split_by_delimiter(buffer, delim);
    os << tokens[0];
    for(std::size_t i=1; i < tokens.size(); ++i)
    {
        os << separator_ << tokens[i];
    }
}

int thermodata::find_specie(const std::string name, std::vector<thnode*>& nodes)
{
    int pos = -1;
    for (auto& node : nodes)
    {
        pos++;
        if (name == node->name_ || name == node->nick_name_ || name == node->original_nick_name_)
            return pos;
    }
    return -1;
}

void thermodata::remove_specie(const std::string name, std::vector<thnode*>& nodes)
{
    int pos = find_specie(name,nodes);
    if (pos > -1)
    {
        nodes.erase(nodes.begin() + pos);
    }
}
