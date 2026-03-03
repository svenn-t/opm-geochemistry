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
// ****************************************************************************
//
//                              POLICY STRUCTS
//
//          See Chapter 1.5 of Andrei Alexandrescu: Modern C++ design.
// ****************************************************************************
#ifndef NUMERICAL_ERROR_POLICY_IS_DEF_H
#define NUMERICAL_ERROR_POLICY_IS_DEF_H

#include <iostream>

template<class Real>
struct DoNothingOnError {

    static constexpr void handleTooManyIterations() {};

    static constexpr void handleBracketingFailure(const Real& a,
                                                  const Real& b,
                                                  const Real& fa,
                                                  const Real& fb)
    {
        // Just to avoid complaints from the compiler about unused function arguments...
        static_assert(a == a);
        static_assert(b == b);
        static_assert(fa == fa);
        static_assert(fb == fb);
    };

    static constexpr void handleNoRoots(const Real& a,
                                        const Real& b,
                                        const Real& c,
                                        const Real& d,
                                        const Real& yn, 
                                        const Real& delta)
    {
        static_assert(a == a);
        static_assert(b == b);
        static_assert(c == c);
        static_assert(d == d);
        static_assert(yn == yn);
        static_assert(delta == delta);
    };

};

template<class Real>
struct ThrowOnError {

    static void handleTooManyIterations() {
        throw std::runtime_error("Too many iterations without convergence!\n");
    };

    static void handleBracketingFailure(const Real& a, const Real& b, const Real& fa, const Real& fb) {
        if (static_cast<double>(fa * fb) > 0) {
            double x0 = static_cast<double>(a);
            double x1 = static_cast<double>(b);
            const std::string err_msg = "Interval [" + std::to_string(x0) + ", " + std::to_string(x1) + "] does not enclose a root of the function\n";
            //const std::string err_msg = "Interval does not enclose a root of the function\n";
            throw std::runtime_error(err_msg);
        }
    };

    static constexpr void handleNoRoots(const Real& a,
                                        const Real& b,
                                        const Real& c,
                                        const Real& d,
                                        const Real& yn, 
                                        const Real& delta)
    {
        const std::string err_msg = "No real roots for a = " + std::to_string(a) + ", " +
                                    "b = " + std::to_string(b) + ", " +
                                    "c = " + std::to_string(c) + ", " +
                                    "d = " + std::to_string(d) + ", " +
                                    "yn = " + std::to_string(yn) + " and " +
                                    "delta = " + std::to_string(delta);
        throw std::runtime_error(err_msg);
    };

};

template<class Real>
struct WarningOnError {

    static void handleTooManyIterations() {
        std::cout << "Warning: Too many iterations without convergence!\n";
    };

    static void handleBracketingFailure(const Real& a,
                                        const Real& b,
                                        const Real& fa,
                                        const Real& fb)
    {
        if (static_cast<double>(fa * fb) > 0)
        {
            double x0 = static_cast<double>(a);
            double x1 = static_cast<double>(b);
            const std::string err_msg = "Interval [" + std::to_string(x0) + ", " + std::to_string(x1) + "] does not enclose a root of the function\n";
            std::cout << "Warning: " + err_msg;
        }
    };

    static constexpr void handleNoRoots(const Real& a,
                                        const Real& b,
                                        const Real& c,
                                        const Real& d,
                                        const Real& yn, 
                                        const Real& delta)
    {
        const std::string err_msg = "No real roots for a = " + std::to_string(a) + ", " +
                                    "b = " + std::to_string(b) + ", " +
                                    "c = " + std::to_string(c) + ", " +
                                    "d = " + std::to_string(d) + ", " +
                                    "yn = " + std::to_string(yn) + " and " +
                                    "delta = " + std::to_string(delta);
            std::cout << "Warning: " + err_msg;
    };

};

#endif
