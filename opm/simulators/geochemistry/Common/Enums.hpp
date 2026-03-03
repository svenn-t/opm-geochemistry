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
/**
 * @file
 */
#ifndef ENUMS_DEF_HPP
#define ENUMS_DEF_HPP

#include <array>
#include <string>

/**
 * Enum struct to differentiate between geochemical components.
 *
*  Note that the (positive) integer values should never be changed, because they are
 * also used as indices for storing numerical data during geochemical simulations
 * (e.g., index 0 is reserved for aqueous complexes and index 1 corresponds to
 * surface complexes).
 */
struct GeochemicalComponentType
{
    static constexpr int NUMBER_OF_NON_MINERAL_TYPES_ = 3;

    static constexpr int MINERAL = -1;
    static constexpr int AQUEOUS_COMPLEX = 0;
    static constexpr int SURFACE_COMPLEX = 1;
    static constexpr int ION_EXCHANGE = 2;
};

struct gFluidPhase  // "g" to avoid nameclash
{
    // NB: Must not change the order.
    static constexpr int UNDEFINED = -1;
    static constexpr int WATER = 0;
    static constexpr int OIL = 1;
    static constexpr int GAS = 2;
};

/**
 * Keywords used to define new geochemical reactions, or modify existing ones.
 */
struct GeochemicalDatabaseKeyword
{
    static inline const std::string BASIS_SPECIES = "BASIS_SPECIES";
    static inline const std::string SECONDARY_SPECIES = "SECONDARY_SPECIES";
    static inline const std::string SURFACE_SPECIES = "SURFACE_SPECIES";
    static inline const std::string EXCHANGE_SPECIES = "EXCHANGE_SPECIES";
    static inline const std::string MINERAL_PHASES = "MINERAL_PHASES";
    static inline const std::string REMOVE_SPECIES = "REMOVE_SPECIES";

    static inline const std::array<std::string, 6> ALL_DATABASE_KEYWORDS = {
        BASIS_SPECIES,
        SECONDARY_SPECIES,
        SURFACE_SPECIES,
        EXCHANGE_SPECIES,
        MINERAL_PHASES,
        REMOVE_SPECIES
    };
};

struct DebugInfoLevel{

    static constexpr int NONE = 0;
    // Not the best names...
    static constexpr int SOME = 1;
    static constexpr int ALOT = 2;
    static constexpr int MAX = 3;

};

enum class DebugInfoTypeMBAL{ BEFORE = 0, FAILED_LINSOLV=1, NO_CONVERGENCE=2 };

#endif  // ENUMS_DEF_HPP
