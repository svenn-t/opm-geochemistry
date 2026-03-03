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
//
//  IMPLEMENTATION OF SECANT METHOD FOR FINDING THE ROOT OF A FUNCTION
//
#ifndef SECANT_METHOD_IS_DEF_H
#define SECANT_METHOD_IS_DEF_H

#include <cassert>
#include <cmath>
#include <exception>
#include <string>

#include <opm/simulators/geochemistry/Numerical/NumericalErrorPolicy.hpp>
#include <opm/simulators/geochemistry/Numerical/NumericalHelperFunctions.hpp>

// See Chapter 1.5.1 of Andrei Alexandrescu: Modern C++ design
// ("template-template parameters" technique)
template<class Real, template<class dummy> class ErrorPolicy=WarningOnError>
class SecantMethodReal : public ErrorPolicy<Real>
{

public:

    template<class Function>
    static constexpr Real solve(const Function& func,
                                const Real a,
                                const Real b,
                                int maxNoIter,
                                double tolerance,
                                int& noIter)
    {
        Real x0 = a;
        Real f0 = func(x0);
        Real x1 = b;
        Real f1 = func(x1);

        if (absoluteValue(f1 - f0) < tolerance) return x1;

        noIter = 0;
        Real x;
        Real fx;
        while(true)
        {
            ++noIter;

            // Compute step length & update solution
            Real dx = f1 * (x1 - x0) / (f1 - f0);
            x = x1 - dx;
            fx = func(x);

            // Check if we are done...
            const auto df = absoluteValue(fx - f1);
            const auto absolute_error = absoluteValue(fx);
            if (absolute_error < tolerance || df < tolerance || noIter == maxNoIter) {
                break;
            }

            // ...otherwise, update iteration variables for the next iteration
            x0 = x1;
            f0 = f1;
            x1 = x;
            f1 = fx;
        }

        const auto absolute_error = absoluteValue(fx);
        if (absolute_error > tolerance) ErrorPolicy<Real>::handleTooManyIterations();
        return x;
    }

    // In case we need more complex convergence criteria, we can call the next two functions with different inputs...

    /**
     * Compute next step size for secant method  for the case when we need to calculate both f0=f(x0) and f1=f(x1); we store those values into the non-const reference arguments.
     */
    template<class Function>
    static constexpr Real computeStepSize(const Function& func,
                                          const Real& x0,
                                          const Real& x1,
                                          Real& f0,
                                          Real& f1)
    {
        f0 = func(x0);
        f1 = func(x1);
        Real dx = f1 * (x1 - x0) / (f1 - f0);
        return dx;
    }


    /** Compute next step size for secant method in the situation when we already know f0=f(x0) and f1=f(x1).  */
    static constexpr Real computeStepSize(const Real& x0,
                                          const Real& x1,
                                          const Real& f0,
                                          const Real& f1)
    {
        return f1 * (x1 - x0) / (f1 - f0);
    }

};

#endif
