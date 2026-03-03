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
#ifndef WATER_H
#define WATER_H

#include <array>
#include <cmath>
#include <iostream>

/*
 * The International Association for the Properties of Water and Steam Lucerne, Switzerland
 * August 2007, Revised Release on the IAPWS Industrial Formulation 1997 (IAPWS-97).
 * (The revision only relates to the extension of region 5 to 50 MPa)
 */
class water
{
    
public:
    
    water();

	void gibbsIAPWS(double T, double P);
	static double PsatIAPWS(double T);

    void printProperties() const;
    
    // ****************************************************************************************
    //                                  PUBLIC VARIABLES
    // ****************************************************************************************
    //                              (exposed to, e.g., eps_JN)
    //
    double v_;  // specific volume [m^3/kg]
    double u_;  // specific internal energy [J/kg]
    double s_;  // specific entropy [J/kg]
    double h_;  // specific enthalpy [J/kg/K]
    double cp_;  // specific isobaric heat capacity [J/kg]
    double cv_;  // specific isochoric heat capacity [J/kg]
    double w_;  // speed of sound [m/s]
    double g_;  // specific Gibbs free energy [J/kg/K]
    
    double G_;  // specific Gibbs free energy [J/kg/K]
    double H_;  // specific enthalpy [J/mol/K]
    double denst_;  // 1/v_
    
    double alpha_;  // Isobaric thermal expansion [1/K]
    double alpha_t_;  //  d(alpha)/dt
    double beta_;  // Isothermal compressibility [1/Pa]
    double Psat_;  //  saturation pressure for given T
    

private:
    
    double P_;  // Pa
    double T_;  // K
    double R_;  // ideal gas constant kJ/kg/K
    double Mw_;  // mol weight H2O
    double Tcrit_;  // K
    double Pcrit_;  // Pa*
    double rho_crit_;  // kg/m^3
    
    /* Calculates Gibbs free energy and thermodynamic constants:
    *
    *   - specific volume [m^3/kg]
    *   - specific internal energy J/kg, specific entropy [J/kg]
    *   - specific enthalpy [J/kg/K]
    *   - specific isobaric heat capacity [J/kg]
    *   - specific isochoric heat capacity [J/kg/K]
    *   - speed of sound [m/s]
    *
    * Units for input pressure and temperature: [P]=Pa, [T]=K.
    */
    void gibbsIAPWSlocal(double T, double P);
    
    /*
    * Calculates integer powers x^n, as well as the first and second derivatives.
    * Stores the results in variables passed in by reference.
    */
    static void mypow(double x, int n, double& xn, double& dxn, double& ddxn);
    
    template<typename real>
    void nth_power(real x, int n, real& xn, real& dxn, real& ddxn)
    {
        xn = std::pow(x, n);
        const auto dbl_n = static_cast<real>(n);
        dxn = dbl_n*xn / x;
        ddxn = (dbl_n - 1.0)*dxn / x;
    }
    
};

#endif
