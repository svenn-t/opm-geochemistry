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
#include <opm/simulators/geochemistry/Database/ChemDatabaseHKF.hpp>

void add_surface_species_to_basis_table(std::vector<std::string>& basis_db)
{
    basis_db.push_back("X-  X   2   -1  0   0   0   0   0   0   0   0   0   0   0   0   0");
    basis_db.push_back("GCO3-   GCO3    1   0   -1  0   0   -140282 -164898 23.53   7.5621  1.1505  1.2346  -2.8266 12.9395 -4.7579 1.2733");
    basis_db.push_back("GCa+    GCa 1   0   1   0   0   -132120 -129800 -13.5   -0.1947 -7.252  5.2966  -2.4792 9   -2.522  1.2366");
    basis_db.push_back("E   E   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0");
    basis_db.push_back("GSiOH   GSiOH   1   0   0   0   0   -199190 -209775 18  1.9 1.7 20  -2.7    29.1    -51.2   0.1291");
}


void add_surface_species_to_complex_table(std::vector<std::string>& aq_db)
{
    aq_db[0].append("\tX-\tGCO3-\tGCa+\tE\tGSiOH");
    for(auto it = aq_db.begin()+1; it != aq_db.end(); ++it)
    {
        (*it).append("\t0\t0\t0\t0\t0");
    }

    // Name    NickName type   charge  scharge a0  DeltaG  DeltaH  S   a1  a2  a3  a4  c1  c2  omega   Ag+ Al+3    Am+3    Ar  Au+ BO2-    Ba+2    Be+2    Bi+3    Br- Ca+2    Cd+2    Ce+3    Cl- Co+3    Cr+3    Cs+ Cu+ Dy+3    Er+3    Eu+3    F-  Fe+3    Fr+ Ga+3    GdO+    H+  H2O HAsO4-2 HCO3-   HPO4-2  He  Hf+4    Hg+2    HoO+    I-  In+3    K+  Kr  La+3    Li+ Lu+3    Mg+2    Mn+2    MoO4-2  NO3-    Na+ NbO3-   Nd+3    Ne  Ni+2    Pb+2    Pd+2    Pm+3    Pr+3    Pt+2    Ra+2    Rb+ Rh+3    Rn  RuO4-2  SO4-2   Sc+3    SeO4-2  SiO2    Sm+3    Sn+2    Sr+2    Tb+3    TcO4-   Tl+3    Tm+3    U+3 VO4-3   WO4-2   Xe  Y+3 Yb+3    Zn+2    Zr+4    e- X- GCO3- GCa+ E GSiOH");
    aq_db.push_back("GCaSO4-  GCaSO4  1   0   -1  0   -312930 -345900 5   2.4079  -1.8992 6.4895  -2.7004 -8.4942 -8.1271 -0.001  0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   2   0");
    aq_db.push_back("GCaCO3-  GCaCO3  1   0   -1  0   -262850 -287390 2.5 -0.3907 -8.7325 9.1753  -2.4179 -11.5309    -9.0641 -0.038  0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   -1  0   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   2   0");
    aq_db.push_back("GCaOH    GCaOH   1   0   0   0   -171300 -179600 6.7 2.7243  -1.1303 6.1958  -2.7322 11.1286 -2.7493 0.4496  0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   -1  1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   1   0");
    aq_db.push_back("GCaHCO3  GCaHCO3 1   0   0   0   -273830 -294350 16  3.706   1.267   5.252   -2.831  41.722  8.336   0.308   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   1   0");
    //
    // Version 1: GCO3H = GCO3 + H+.
    aq_db.push_back("GCO3H  GCO3H   1   0   0   0   -146940.8956    -98900  28.1    6.2466  7.4711  2.8136  -3.0879 38.4529 5.3534  -0.1934 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   -1  0");
    // Version 2: correct for different ratio of H2CO3 and CO2(aq) and added H2O since CO2 = HCO3- -H2O + H+ ---> GCO3H = GCO3- + H+ - H2O.
//    aq_db.push_back("GCO3H GCO3H   1   0   0   0   -92250  -98900  28.1    15.2964 -10.0382    -55.4493    7.4092  36.8069 3.5851  -0.3107 0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   -1  0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   -1  0");
    //
    aq_db.push_back("GCO3Ca+  GCO3Ca  1   0   1   0   -274776.5694    -294350 16  3.706   1.267   5.252   -2.831  41.722  8.336   0.308   0   0   0   0   0   0   0   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   -2  0");
    aq_db.push_back("GCO3Mg+  GCO3Mg  1   0   1   0   -251147.8572    -275750 -3  3.271   0.206   5.669   -2.788  47.284  9.34    0.599   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   -2  0");
    aq_db.push_back("GCO3Ba+  GCO3Ba  1   0   1   0   -275700 -294350 46.83   7.882   11.464  1.244   -3.253  26.273  4.46    -0.158  0   0   0   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   -2  0");
    aq_db.push_back("GCO3Sr+  GCO3Sr  1   0   1   0   -276720 -291688 38  3.694   1.237   5.264   -2.83   35.345  7.181   -0.023  0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   -2  0");
    aq_db.push_back("NaX  NaX 2   1   0   4   -89844.71389    0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   0   0   0");
    aq_db.push_back("KX   KX  2   1   0   3   -95714.51176    0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   0   0   0");
    aq_db.push_back("CaX2 CaX2    2   2   6   0   -187718.0261    0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   2   0   0   0   0");
    aq_db.push_back("MgX2 MgX2    2   2   8   0   -163842.6066    0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   2   0   0   0   0");
    aq_db.push_back("BaX2 BaX2    2   2   5   0   -189776.2092    0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   2   0   0   0   0");
    aq_db.push_back("SrX2 SrX2    2   2   5   0   -190506.6554    0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0   2   0   0   0   0");

    //Original HKF values
    aq_db.push_back("GSiOMg+  GSiOMg  1   0   1   0   -353025 -385700 -23.78  0.6289  -6.2428 8.1967  -2.5209 36.7882 4.6702  0.9177  0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   -1  1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   -1  1");
    aq_db.push_back("GSiOCa+  GSiOCa  1   0   1   0   -376299 -403109 -1.99   1.0647  -5.1787 7.7785  -2.5649 30.8048 3.6619  0.5831  0   0   0   0   0   0   0   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   -1  1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   -1  1");
    aq_db.push_back("GSiO-    GSiO    1   0   -1  0   -246000 -273872 5   2.9735  -0.5158 5.9467  -2.7575 8.1489  -7.3123 1.5511  0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   -1  1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   1");
    aq_db.push_back("GSiONa   GSiONa  1   0   0   0   -314000 -333894 10  3.4928  0.75    5.4483  -2.81   20.2395 1.9785  -0.038  0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   -1  1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1");

    // "Matched" Revil et al.
    // aq_db.push_back("GSiOMg+    GSiOMg  1   0   1   0   -353025 0.  0.  0.  0. 0.   0. 0.   0.  0.  0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   -1  1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   -1  1");
    // aq_db.push_back("GSiOCa+    GSiOCa  1   0   1   0   -376299 0.  0.  0.  0. 0.   0. 0.   0.  0.  0   0   0   0   0   0   0   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   -1  1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   -1  1");
    // aq_db.push_back("GSiO-  GSiO    1   0   -1  0   -246000 0.  0.  0.  0. 0.   0. 0.   0.  0.  0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   -1  1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   1");
    // aq_db.push_back("GSiONa GSiONa  1   0   0   0   -314000 0.  0.  0.  0. 0.   0. 0.   0.  0.  0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   -1  1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1");
}


void add_surface_species_to_minerals_table(std::vector<std::string>& mineral_db)
{
    mineral_db[0].append("\tX-\tGCO3-\tGCa+\tE\tGSiOH");
    for (auto it = mineral_db.begin() + 1; it != mineral_db.end(); ++it)
    {
        (*it).append("\t0\t0\t0\t0\t0");
    }
}
