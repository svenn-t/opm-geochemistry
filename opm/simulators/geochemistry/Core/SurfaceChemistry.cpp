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
#include <opm/simulators/geochemistry/Core/SurfaceChemistry.hpp>

double DiffuseLayerProperties::calcIonicStrength() const
{
    return 0.5 * (n1_ + n_1_) + 2.0 * (n2_ + n_2_);
}

double DiffuseLayerProperties::calcIonicStrengthReplacementFactor(double E) const
{
    return n1_ + n2_*(2.0 + 1.0 / E) + n_2_*E;
}

double DiffuseLayerProperties::calcSquareOfDiffuseLayerIntegrandDenominator(double E) const
{
    return E*E*(n1_ * (1.0 / E - 1.0) + n2_ * (1.0 / E / E - 1.0) + n_1_ * (E - 1.0) + n_2_ * (E * E - 1.0));
}

double DiffuseLayerProperties::calcTotalSurfaceExcessConcentration() const
{
    // Until 18/7-22: Also multiplied with frac_DL_...
    return fDL1_*n1_ + 2.0*fDL2_*n2_ - fDL_1_*n_1_ - 2.0*fDL_2_*n_2_;
}

void DiffuseLayerProperties::setConcentrationsToZero()
{
    n1_ = 0.0;
    n_1_ = 0.0;
    n2_ = 0.0;
    n_2_ = 0.0;
}

void DiffuseLayerProperties::setSurfaceExcessesOfNegativeIonsToZero()
{
    fDL_1_ = 0.0;
    fDL_2_ = 0.0;
}

void DiffuseLayerProperties::setSurfaceExcessesOfPositiveIonsToZero()
{
    fDL1_ = 0.0;
    fDL2_ = 0.0;
}

void DiffuseLayerProperties::setSurfaceExcessesToZero()
{
    setSurfaceExcessesOfNegativeIonsToZero();
    setSurfaceExcessesOfPositiveIonsToZero();
}