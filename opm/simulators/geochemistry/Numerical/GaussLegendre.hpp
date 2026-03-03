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
#ifndef GAUSS_LEGENDRE_IS_DEF_H
#define GAUSS_LEGENDRE_IS_DEF_H

#include <opm/simulators/geochemistry/Numerical/GaussLegendrePrecomputedTables.hpp>

/* Class for doing numerical integration with Gauss-Legendre quadrature.
*  Strongly inspired by the implementation in the BOOST library.
*/
template<typename Real, unsigned N>
struct GaussLegendreIntegrator: public GaussLegendreTables<Real, N>
{
    using base = GaussLegendreTables<Real, N>;
    using value_type = Real;
    
public:
    
    template<class Function>
    static Real integrate(Function& f, Real lower_integration_limit, Real upper_integration_limit){
        
        const value_type A = 0.5*(upper_integration_limit - lower_integration_limit);
        const value_type B = 0.5*(upper_integration_limit + lower_integration_limit);
        
        value_type result = 0.0;
        
        // Must handle the case of odd order separately
        const bool isOddNumberOfPoints = (N&1);
        unsigned startIndexForLoop = 0;
        if(isOddNumberOfPoints){
            result += f(B)*base::weights()[0];
            ++startIndexForLoop;
        }
        
        for(unsigned i=startIndexForLoop; i < base::abscissas().size(); ++i){
            const value_type Ax = A * base::abscissas()[i];
            const value_type f_pos = f(B+Ax);
            const value_type f_min = f(B-Ax);
            result += (f_pos + f_min) * base::weights()[i];
        }
        result *= A;  // scale the answer to the range of integration
        return result;
    }
    
};

#endif
