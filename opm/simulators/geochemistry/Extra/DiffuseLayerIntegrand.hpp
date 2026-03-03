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
#ifndef DIFFUSE_LAYER_INTEGRAND_IS_DEF_H
#define DIFFUSE_LAYER_INTEGRAND_IS_DEF_H

#include <iostream>

/*
* Class for evaluating the integrals used in the explicit diffuse layer model,
* implemented using the Curiously Recurring Template Pattern (CRTP).
*
* Note that this class is "standalone": it is not used anywhere inside the
* geochemical solver, only for verifying the integration formulas by comparing
* them to numerical integration.
*/
template<typename Derived>
class DiffuseLayerIntegrandBase{

public:

    DiffuseLayerIntegrandBase(double n1, double n2, double n_2)
    : n1_(n1)
    , n2_(n2)
    , n_2_(n_2)
    {
        n_1_ = n1_ + 2.0*(n2_ - n_2_);
        A_ = n_2_;
        B_ = n1_ + 2.0*n2_;
        C_ = n2_;
        discriminant_ = B_*B_ - 4.0*A_*C_;
        only_monovalent_ = (n2_ == 0.0) && (n_2_ == 0.0);
    }

    double operator()(double E){
        return derived(E);
    }

    double integrate(double E1, double E2){
        return derived(E1, E2);
    }

    double discriminant(){
        return discriminant_;
    }

private:

    friend class DiffuseLayerIntegrandOne;
    friend class DiffuseLayerIntegrandTwo;
    friend class DiffuseLayerIntegrandMinusOne;
    friend class DiffuseLayerIntegrandMinusTwo;
    DiffuseLayerIntegrandBase(){}; // prevent default constructor

    // Concentration variables
    double n1_;
    double n2_;
    double n_2_;
    double n_1_;  // NB: calculated by assuming charge-balance

    // Helper variables
    double A_;
    double B_;
    double C_;
    double discriminant_;
    bool only_monovalent_;

    Derived& derived(){
        return *static_cast<Derived*>(this);
    }

};

/* Case when Z=+1. */
class DiffuseLayerIntegrandOne: public DiffuseLayerIntegrandBase<DiffuseLayerIntegrandOne>
{
    using base = DiffuseLayerIntegrandBase<DiffuseLayerIntegrandOne>;

public:

    DiffuseLayerIntegrandOne(double n1, double n2, double n_2)
    : base(n1, n2, n_2)
    {
    }

    double operator()(double E){
        if (E == 1 || n1_ == 0.0){
            return 0.0;
        }
        else{
            return -1. / (E*std::sqrt(n1_*E + n2_*(2.0*E + 1.0) + n_2_*E*E));
        }
    }

    double integrate(double E1, double E2){

        const auto& discriminant = discriminant_;
        const auto& A = A_;
        const auto& B = B_;
        const auto& C = C_;

        if(n1_ == 0.0){
            return 0.0;
        }
        else if(only_monovalent_){  // B == n1_
            return (2.0/std::sqrt(B))*(std::sqrt(1/E2)-std::sqrt(1/E1));
        }
        else if(A == 0){  // no ions with charge -2 (striclyy speaking, not needed...)
            const auto lnarg_2 = B + 2.0*C/E2  + 2.0*std::sqrt(B*E2 + C)*std::sqrt(C)/E2;
            const auto lnarg_1 = B + 2.0*C/E1  + 2.0*std::sqrt(B*E1 + C)*std::sqrt(C)/E1;
            return std::log(lnarg_2/lnarg_1)/sqrt(C);
        }
        else if (discriminant > 0 && C > 0){
            // ln(B + 2*sqrt(A*E^2 + B*E + C)*sqrt(C)/E + 2*C/E)/sqrt(C)
            const auto sqrt_quad_2 = std::sqrt(A*E2*E2 + B*E2 + C);
            const auto sqrt_quad_1 = std::sqrt(A*E1*E1 + B*E1 + C);
            const auto lnarg_2 = B + 2.0*C/E2 + 2.0*sqrt_quad_2*std::sqrt(C)/E2;
            const auto lnarg_1 = B + 2.0*C/E1 + 2.0*sqrt_quad_1*std::sqrt(C)/E1;
            return std::log(lnarg_2/lnarg_1)/sqrt(C);
        }
        else if(discriminant > 0 && C == 0){
            return (2.0/B)*(std::sqrt(A*E2*E2 + B*E2)/E2 - std::sqrt(A*E1*E1+B*E1)/E1);
        }
        else{
            // arcsinh(B/sqrt(-B^2 + 4*A*C) + 2*C/(sqrt(-B^2 + 4*A*C)*E))/sqrt(C)
            const auto asinh_2 = std::asinh(B/std::sqrt(-discriminant) + 2.0*C/E2/std::sqrt(-discriminant));
            const auto asinh_1 = std::asinh(B/std::sqrt(-discriminant) + 2.0*C/E1/std::sqrt(-discriminant));
            return (asinh_2-asinh_1)/sqrt(C);
        }
    }
};

/* The case when Z=+2. */
class DiffuseLayerIntegrandTwo: public DiffuseLayerIntegrandBase<DiffuseLayerIntegrandTwo>
{
    using base = DiffuseLayerIntegrandBase<DiffuseLayerIntegrandTwo>;

public:

    DiffuseLayerIntegrandTwo(double n1, double n2, double n_2)
    : base(n1, n2, n_2)
    {
    }

    double operator()(double E){
        if (E == 1 || n2_ == 0.0){
            return 0.0;
        }
        else{
            return -(E + 1.0) / (E*E*std::sqrt(n1_*E + n2_*(2.0*E + 1.0) + n_2_*E*E));
        }
    }

    double integrate(double E1, double E2){

        const auto& discriminant = discriminant_;
        const auto& A = A_;
        const auto& B = B_;
        const auto& C = C_;

        const auto sqrt_quad_2 = std::sqrt(A*E2*E2 + B*E2 + C);
        const auto sqrt_quad_1 = std::sqrt(A*E1*E1 + B*E1 + C);

        if(n2_ == 0.0){  // C == 0
            return 0.0;
        }
        else if(only_monovalent_){
            return 0.0;
        }
        else if(A == 0){ // no ions of charge -2 (strictly speaking, not needed...)
            // First compute f1...
            const auto lnarg_2 = B + 2.0*C/E2  + 2.0*std::sqrt(B*E2 + C)*std::sqrt(C)/E2;
            const auto lnarg_1 = B + 2.0*C/E1  + 2.0*std::sqrt(B*E1 + C)*std::sqrt(C)/E1;
            const auto f1 = std::log(lnarg_2/lnarg_1)/sqrt(C);
            //... then use that to compute the final result:
            return (1.0/C)*(std::sqrt(B*E2+C)/E2 - std::sqrt(B*E1+C)/E1) + (1.0-0.5*B/C)*f1;
        }
        else if(discriminant > 0 && C > 0.0){
            // (1-0.5*B/C)/sqrt(C) * ln(B + 2*sqrt(A*E^2 + B*E + C)*sqrt(C)/E + 2*C/E) + sqrt(A*E^2 + B*E + C)/(C*E)
            const auto lnarg_2 = B + 2.0*sqrt_quad_2*std::sqrt(C)/E2 + 2.0*C/E2;
            const auto lnarg_1 = B + 2.0*sqrt_quad_1*std::sqrt(C)/E1 + 2.0*C/E1;
            return std::log(lnarg_2/lnarg_1)*(1.0-0.5*B/C)/std::sqrt(C) + (sqrt_quad_2/E2/C-sqrt_quad_1/E1/C);
        }
        else{
            // (1-0.5*B/C)/sqrt(C) * arcsinh(B/sqrt(-B^2 + 4*A*C) + 2*C/(sqrt(-B^2 + 4*A*C)*E)) + sqrt(A*E^2 + B*E + C)/(C*E)
            const auto asinh_2 = std::asinh(B/std::sqrt(-discriminant) + 2.0*C/E2/std::sqrt(-discriminant));
            const auto asinh_1 = std::asinh(B/std::sqrt(-discriminant) + 2.0*C/E1/std::sqrt(-discriminant));
            // Note that C cannot be zero when the discriminant is negative.
            return (asinh_2-asinh_1)*(1.0-0.5*B/C)/std::sqrt(C) + sqrt_quad_2/E2/C - sqrt_quad_1/E1/C;
        }
    }
};

/* The case when Z=-1. */
class DiffuseLayerIntegrandMinusOne: public DiffuseLayerIntegrandBase<DiffuseLayerIntegrandMinusOne>
{
    using base = DiffuseLayerIntegrandBase<DiffuseLayerIntegrandMinusOne>;

public:

    DiffuseLayerIntegrandMinusOne(double n1, double n2, double n_2)
    : base(n1, n2, n_2)
    {
    }

    double operator()(double E){
        if (E == 1 || n_1_ == 0.0){
            return 0.0;
        }
        else{
            return 1.0 / std::sqrt(n1_*E + n2_*(2.0*E + 1.0) + n_2_*E*E);
        }
    }

    double integrate(double E1, double E2){

        const auto& discriminant = discriminant_;
        const auto& A = A_;
        const auto& B = B_;
        const auto& C = C_;

        const auto sign = discriminant > 0.0 ? 1.0 : -1.0;
        const auto D = 0.5*sqrt(sign*discriminant)/A;
        const auto u2 = E2 + 0.5*B/A;
        const auto u1 = E1 + 0.5*B/A;

        if(n_1_ == 0.0){
            return 0.0;
        }
        else if(only_monovalent_){  // B == n1_
            return (2.0/std::sqrt(B))*(std::sqrt(E2)-std::sqrt(E1));
        }
        else if(A == 0){  // Ions of charge -2
            // Note: In this case, a separate if clause is needed, since we
            //       divide by A otherwise...
            return (2.0/B)*(std::sqrt(B*E2+C)-std::sqrt(B*E1+C));
        }
        else if(discriminant > 0 && C > 0.0){
            return (std::acosh(u2/D)-std::acosh(u1/D))/std::sqrt(A);
        }
        else if(discriminant > 0 && C == 0.0){
            const auto numer = 2.0*A*E2 + B + 2.0*std::sqrt(A*(A*E2*E2+B*E2));
            const auto denom = 2.0*A*E1 + B + 2.0*std::sqrt(A*(A*E1*E1+B*E1));
            return std::log(numer/denom)/std::sqrt(A);
        }
        else{
            return (std::asinh(u2/D)-std::asinh(u1/D))/std::sqrt(A);
        }
    }
};

/* The case when Z=-2. */
class DiffuseLayerIntegrandMinusTwo: public DiffuseLayerIntegrandBase<DiffuseLayerIntegrandMinusTwo>
{
    using base = DiffuseLayerIntegrandBase<DiffuseLayerIntegrandMinusTwo>;

public:

    DiffuseLayerIntegrandMinusTwo(double n1, double n2, double n_2)
    : base(n1, n2, n_2)
    {
    }

    double operator()(double E){
        if (E == 1 || n_2_ == 0.0){
            return 0.0;
        }
        else{
            return (E + 1.0) / std::sqrt(n1_*E + n2_*(2.0*E + 1.0) + n_2_*E*E);
        }
    }

    double integrate(double E1, double E2){

        const auto& discriminant = discriminant_;
        const auto& A = A_;
        const auto& B = B_;
        const auto& C = C_;

        const auto sign = discriminant > 0.0 ? 1.0 : -1.0;
        const auto D = 0.5*sqrt(sign*discriminant)/A;

        const auto u2 = E2 + 0.5*B/A;
        const auto u1 = E1 + 0.5*B/A;
        const auto sqrt_diff_2 = std::sqrt(u2*u2 - sign*D*D);
        const auto sqrt_diff_1 = std::sqrt(u1*u1 - sign*D*D);

        if(n_2_ == 0.0){
            return 0.0;  // no ions of charge -2 (A==0)
        }
        else if(only_monovalent_){
            return 0.0;
        }
        else if(discriminant > 0 && C > 0.0){
            const auto d_acosh = std::acosh(u2/D)-std::acosh(u1/D);
            return (sqrt_diff_2-sqrt_diff_1 + (1.0-0.5*B/A)*d_acosh)/std::sqrt(A);
        }
        else if(discriminant > 0 && C == 0.0){
            // First calculate f(-1)....
            const auto numer = 2.0*A*E2 + B + 2.0*std::sqrt(A*(A*E2*E2+B*E2));
            const auto denom = 2.0*A*E1 + B + 2.0*std::sqrt(A*(A*E1*E1+B*E1));
            const auto f_m1 = std::log(numer/denom)/std::sqrt(A);
            // ... and use it to compute the final result:
            return (1.0/A)*(std::sqrt(A*E2*E2+B*E2)-std::sqrt(A*E1*E1+B*E1)) + (1-0.5*B/A)*f_m1;
        }
        else{
            const auto d_asinh = std::asinh(u2/D)-std::asinh(u1/D);
            return (sqrt_diff_2-sqrt_diff_1 + (1.0-0.5*B/A)*d_asinh)/std::sqrt(A);
        }
    }
};

#endif
