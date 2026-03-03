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
#ifndef SURFACE_CHEMISTRY_IS_DEF_HPP
#define SURFACE_CHEMISTRY_IS_DEF_HPP

// Note: All three options assume that there can only be -1, +1, -2, +2 charged species
enum class GrahameEquation { GENERAL = 0, ASSUME_CHARGE_BALANCE = 1, SYMMETRICAL_ELECTROLYTE = 2 };

enum class DiffuseLayerIntegrationMethod { MEAN_POTENTIAL_APPROXIMATION = 0, BORKOVEC_WESTALL = 1, ANALYTICAL=2, AKSEL_HIORTH = 3 };

struct SurfaceChemistryOptions
{
    int smethod_{0};

    bool INCLUDE_SURFACE_EXCESSES_IN_SURFACE_EQUATION = true;  // !TESTMODE

    GrahameEquation GRAHAME_EQ_VERSION = GrahameEquation::ASSUME_CHARGE_BALANCE;
    double VALENCE_OF_SYMMETRICAL_ELECTROLYTE = 1.0;  // only used for the symmetrical electrolyte case

    bool CALCULATE_DIFFUSE_LAYER = false;
    bool ONLY_COUNTER_IONS = false;
    bool EXCHANGE_MODEL_DL = false;

    // Only relevant if CALCULATE_DIFFUSE_LAYER == true
    DiffuseLayerIntegrationMethod INTEGRATION_METHOD = DiffuseLayerIntegrationMethod::AKSEL_HIORTH;
};

/** For holding total concentrations of species of a given charge in the
 *  diffuse layer, as well as surface excesses and the volume fraction
 *  of the diffuse layer. */
struct DiffuseLayerProperties
{
    double n1_{ 0.0 };      // Total (bulk) concentration of +1 species.
    double n_1_{ 0.0 };     // Total (bulk) concentration of -1 species.
    double n2_{ 0.0 };      // Total (bulk) concentration of +2 species.
    double n_2_{ 0.0 };     // Total (bulk) concentration of -2 species.
    double fDL1_{ 0.0 };    // g(1) - Ratio of surface excess to bulk concentration for +1 species.
    double fDL_1_{ 0.0 };   // g(-1) - Ratio of surface excess to bulk concentration for -1 species.
    double fDL2_{ 0.0 };    // g(2) - Ratio of surface excess to bulk concentration for +2 species.
    double fDL_2_{ 0.0 };   // g(-2) - Ratio of surface excess to bulk concentration for -2 species.
    double frac_DL_{ 1.0 }; // Volume fraction of the diffuse layer.

    /** @return Ionic strength calculated from n(1), n(2), n(-1), and n(-2).*/
    double calcIonicStrength() const;

    /**
     *  n(-1) is eliminated by assuming charge-balance.
     *  @return Ionic strength calculated from n(1), n(2), and n(-2).
     *  */
    double calcIonicStrengthReplacementFactor(double E) const;

    double calcSquareOfDiffuseLayerIntegrandDenominator(double E) const;
    double calcTotalSurfaceExcessConcentration() const;

    void setConcentrationsToZero();
    void setSurfaceExcessesOfNegativeIonsToZero();
    void setSurfaceExcessesOfPositiveIonsToZero();
    void setSurfaceExcessesToZero();
};

#endif