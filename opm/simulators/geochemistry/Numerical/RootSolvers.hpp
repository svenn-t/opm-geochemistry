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
#ifndef ROOT_SOLVERS_IS_DEF_H
#define ROOT_SOLVERS_IS_DEF_H

#include "NumericalErrorPolicy.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <stdexcept>
#include <vector>


template<class Real, template<class dummy> class ErrorPolicy=WarningOnError>
class LinearRootSolver : public ErrorPolicy<Real>
{
public:
    static Real solve(Real a, Real b)
    {
        // If a = 0, return error
        if ( std::abs(a) < 1e-16 ) {
            ErrorPolicy<Real>::handleNoRoots(0.0, 0.0, a, b, 0.0, 0.0);
        }

        return -b / a;
    }
}; // class LinearRootSolver

template<class Real, template<class dummy> class ErrorPolicy=WarningOnError>
class QuadraticRootSolver : public ErrorPolicy<Real>
{
public:
    static std::vector<Real> solve(Real a, Real b, Real c)
    {
        // Init. output vector
        std::vector<Real> roots;

        // Check if it's really a linear equation
        if ( std::abs(a) < 1e-16 ) {
            Real z1 = LinearRootSolver<Real, ErrorPolicy>::solve(b, c);
            roots.push_back(z1);    
            return roots;
        }

        // Check discriminant for two, one (two repeated), or zero solutions
        Real disc = b * b - 4 * a * c;
        if ( disc > 0 ) {
            // Quadratic formula
            Real disc_sqrt = std::sqrt(disc);
            Real z1 = (-b + disc_sqrt) / (2 * a); 
            Real z2 = (-b - disc_sqrt) / (2 * a);

            // Insert in roots
            roots.push_back(z1);
            roots.push_back(z2);
        }
        else if ( disc < 0 ) {
            ErrorPolicy<Real>::handleNoRoots(0.0, a, b, c, 0.0, 0.0);
        }
        else {
            // disc = 0, one simple solution
            Real z1 = -b / (2 * a);
            roots.push_back(z1);
        }

        return roots;
    }
};  // class QuadraticRootSolver

/*
* Cubic root solver from paper:
* Nickalls R. W. D. (1993) A new approach to solving the cubic: Cardan’s solution revealed. Math. Gazette 77, 354–359
*/
template<class Real, template<class dummy> class ErrorPolicy=WarningOnError>
class CubicRootSolver : public ErrorPolicy<Real>
{
public:
    static constexpr Real PI = 3.14159265358979323846;
    static std::vector<Real> solve(Real a, Real b, Real c, Real d)
    {
        // Check if it's really a quadratic equation
        if ( std::abs(a) < 1e-16 ) {
            return QuadraticRootSolver<Real, ErrorPolicy>::solve(b, c, d);
        }

        // Helper function for cubic equation eval
        auto f = [&] (Real x) { return a * x * x * x + b * x * x + c * x + d; };

        // Temp. calculations for later checks
        Real xn = -b / (3 * a);
        Real yn = f(xn);
        Real yn_sq = yn * yn;
        Real delta_sq = (b * b - 3 * a * c) / (9 * a * a);
        Real h_sq = 4 * a * a * delta_sq * delta_sq * delta_sq;

        // Init. output vector
        std::vector<Real> roots;

        // Three checks with different outcomes:
        // (i)   yn_sq > h_sq: 1 real root,
        // (ii)  yn_sq = h_sq: 3 real roots where two or all are equal
        // (iii) yn_sq < h_sq: 3 distinct real roots
        if ( yn_sq > h_sq ) {
            // Intermediate calculations
            Real yn_h_sq = std::sqrt(yn_sq - h_sq);
            Real inv_2a = 1 / (2 * a);
            Real cub_term_1 = inv_2a * (-yn + yn_h_sq); 
            Real cub_term_2 = inv_2a * (-yn - yn_h_sq);
            
            // Root
            Real z1 = xn + std::cbrt(cub_term_1) + std::cbrt(cub_term_2);
            roots.push_back(z1);
        }
        else if ( yn_sq < h_sq ) {
            // To get real roots yn / h must be in [-1, 1] since cos(3*theta) = -yn / h
            Real h = std::sqrt(h_sq);
            Real delta = std::sqrt(delta_sq);
            Real yn_div_h = yn / h;
            if ( yn_div_h >= -1 && yn_div_h <= 1 ) {
                // Calculate angle theta
                Real theta = std::acos(-yn / h) / 3;

                // Three distinct roots
                Real z1 = xn + 2 * delta * std::cos(theta);
                Real z2 = xn + 2 * delta * std::cos(theta + 2 * PI / 3);
                Real z3 = xn + 2 * delta * std::cos(theta + 4 * PI / 3);

                // Insert roots in vector
                roots.push_back(z1);
                roots.push_back(z2);
                roots.push_back(z3);
            }
            else {
                ErrorPolicy<Real>::handleNoRoots(a, b, c, d, yn, delta);
            }
        }
        else {
            // Roots are xn + delta twice and xn - 2 * delta. If h = 0 then delta = 0, thus all roots are equal to xn
            if ( h_sq < 1e-12 ) {
                roots.push_back(xn);
            }
            else {
                // delta must be using yn (see paper)
                Real delta = std::cbrt(yn / (2 * a));

                // Insert roots in vector
                roots.push_back(xn + delta);
                roots.push_back(xn - 2 * delta);
            }
        }
        
        // Return roots in descending order if more than one root
        if ( roots.size() > 1 ) {
            std::sort(roots.begin(), roots.end(), std::greater<Real>());
        }

        return roots;
    }
};  // class CubicRootSolver

#endif