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
#ifndef HKF_H
#define HKF_H

#include <opm/simulators/geochemistry/Thermo/eps_JN.h>
#include <opm/simulators/geochemistry/Thermo/ions.h>
#include <opm/simulators/geochemistry/Thermo/thermodata.h>
#include <opm/simulators/geochemistry/Thermo/water.h>

class hkf
{
public:
    
    hkf();
    
    void WaterProp(double T, double P);
    
    void dGMineral(double T, double P, double* G, double* S, double* a, double* c, double* b, double* V, int size, double* logK);
    
    void dGIons(double T, double P, double* G, double* S,
                double* a1, double* a2, double* a3, double* a4,
                double* c1, double* c2,
                double* omega, double* Z, double* re_ref, int size,
                int skip, double* dG, double* MV);
    
    double rhow_;
    double epsw_;
    double G_;  // Gibbs free energy ref Tripple point of H2O
    double H_;  // Entalpy  ref Triple point of H2O
    double S_;  // Entropy  ref Triple point of H2O
    
private:
    
    std::unique_ptr<water> W_;
    std::unique_ptr<eps_JN> EPS_;
    std::unique_ptr<ions> BORN_;
    
    static constexpr double Rg_inv_ = 0.12027239580856473682150703284415;
    static constexpr double loge_ = 0.43429448190325182765112891891661;
    static constexpr double ln10_ = 2.3025850929940456840179914546844;
    
    double Tref_;
    double Tref2_inv_;
    double Pref_;
    double Pref_inv_;
    double Theta_;
    double Psi_;
    
    double Y_ref_;
    double Z_ref_;
    
    double Gtr_;
    double Htr_;  // not used?
    double Str_;
    double Ttr_;
    
};

#endif
