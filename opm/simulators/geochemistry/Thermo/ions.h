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
#ifndef IONS_H
#define IONS_H

#include <opm/simulators/geochemistry/Thermo/eps_JN.h>
#include <opm/simulators/geochemistry/Thermo/water.h>

class ions
{
    
public:
    
    ions(water* W);
    
    void compare_analytical_and_numerical_derivatives(double T = 200 + 273.15, double P = 500e5, double dT = .01, double dP = 1.0);
    
    void born(double Z, double& re_ref, double& w, double& w_T, double& w_TT, double& w_P);
    void born(double Z, double& re_ref, double& w);
    void born_f(double T);
    void born_df(double T, double P);
    
private:
    
    water* W_;
    
    double eta_;
    double born_f_;
    double born_f_P_;
    double born_f_T_;
    double born_f_TT_;
    double born_g_;
    double born_g_P_;
    double born_g_T_;
    double born_g_TT_;
    double born_gp2_;
    double u_;
    double u_T_;
    double u_TT_;
    
};


#endif
