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
// IMPLEMENTATION OF BRENT'S METHOD FOR FINDING THE ROOT OF A FUNCTION.
// ESSENTIALLY, THE SAME IMPLEMENTATION AS IN "NUMERICAL RECIPES IN C/C++".
//
//
// IMPORTANT NOTE (!!!): THE CODE HAS SO FAR BARELY BEEN TESTED...
//
#ifndef BRENT_METHOD_IS_DEF_H
#define BRENT_METHOD_IS_DEF_H

#include <limits>

#include <opm/simulators/geochemistry/Numerical/NumericalErrorPolicy.hpp>
#include <opm/simulators/geochemistry/Numerical/NumericalHelperFunctions.hpp>

namespace Brent
{
    template<class Float, class Real>
    auto static constexpr SIGN(const Float& a, const Real& b) -> Float
    // Use case: a=tol, b=xm
    {
        if(b >= 0)
        {
            return (a >= 0 ? a : -a);
        }
        else
        {
            return (a >= 0 ? -a : a);
        }
    }
}

// See Chapter 1.5.1 of Alexandrescu: Modern C++ ("template-template parameters" technique)
template<class Real, template<class dummy> class ErrorPolicy=WarningOnError>
class BrentMethodReal : public ErrorPolicy<Real>
{

public:

    template<class Function>
    static constexpr Real solve(const Function& func,
                                const Real x1,
                                const Real x2,
                                int maxNoIter,
                                double tolerance,
                                int& noIter)
    {
        using namespace Brent;

        // NB: Precludes using integer type for the template argument "Real"
        constexpr auto machine_eps = std::numeric_limits<Real>::epsilon();

        Real a = x1;
        Real b = x2;
        Real c = x2;
        Real fa = func(a);
        Real fb = func(b);

        if(fa*fb > 0.0)  // f(a) and f(b) have the same sign
        {
            ErrorPolicy<Real>::handleBracketingFailure(a, b, fa, fb);
            return b;  // Do nothing, but choose x2 as "best guess"...
        }

        Real fc = fb;
        Real d = 0.0;
        Real e = 0.0;
        while(noIter < maxNoIter)
        {
            if(fb*fc > 0.0)  // f(b) and f(c) have the same sign
            {
                // Note: This will always happen during the very first iteration,
                //       ensuring that b != c before proceeding.
                c = a;
                fc = fa;
                d = b-a;
                e = d;
            }

            if (absoluteValue(fc) < absoluteValue(fb))  // |f(c)| < |f(b)|
            {
                a = b;
                b = c;
                c = a;
                fa = fb;
                fb = fc;
                fc = fa;
            }

            const auto tol = 2.0*machine_eps*absoluteValue(b) + 0.5*tolerance;
            const auto xm = 0.5*(c-b);

            if (absoluteValue(xm) <= tol || fb == 0.0) return b;

            if (absoluteValue(e) >= tol && absoluteValue(fa) > absoluteValue(fb))
            {
                const Real s = fb/fa;

                Real p;
                Real q;
                Real r;

                if(a == c)
                {
                    p = 2.0*xm*s;
                    q = 1.0 - s;
                }
                else
                {
                    q = fa/fc;
                    r = fb/fc;
                    p = s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
                    q = (q-1.0)*(r-1.0)*(s-1.0);
                }

                if (p > 0.0) q = -q;
                p = absoluteValue(p);

                const Real min1 = 3.0*xm*q - absoluteValue(tol*q);
                const Real min2 = absoluteValue(e*q);
                if (2.0*p < (min1 < min2 ? min1 : min2))
                {
                    e = d;
                    d = p/q;
                } else
                {
                    d = xm;
                    e = d;
                }
            }
            else
            {
                d = xm;
                e = d;
            }

            a = b;
            fa = fb;

            if (absoluteValue(d) > tol) b += d;
            else
            {
                b += SIGN(tol, xm);
            }
            // NB: In the original code, fb=func(b) had the same level of
            //     indentation as the statement "b += SIGN(tol, xm)"...
            //     However, without braces, this is highly misleading, because
            //     it falsely suggests that the function evaluation takes place
            //     inside the else block.
            fb = func(b);

            ++noIter;
            //std::cout << "Finished iteration #" << noIter << "...\n";
        }  // end while

        ErrorPolicy<Real>::handleTooManyIterations();
        return b;  // current best iterate

    }  // end solve
};

#endif
