/*
  Copyright 2026 NORCE Norwegian Research Centre AS

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef OPM_UPWINDTRANSPORTSOLVER_HEADER_INCLUDED
#define OPM_UPWINDTRANSPORTSOLVER_HEADER_INCLUDED

#include <opm/common/utility/SparseTable.hpp>

#include <opm/simulators/utils/Tarjan.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

/// \file
///
/// Solver for the linear system of an implicit upwind transport scheme, cell
/// by cell in the order of the flow (a "reordering" solver).
///
/// For every cell `c` and every transported quantity `q` (e.g. a species)
/// the system is
///
///     diagonal[c] * x_q[c] = rhs_q[c] + weightFactor * sum_k fluxes[c][k] * x_q[upwind[c][k]]
///
/// where `upwind[c]` are the cells that `c` receives flow from and
/// `fluxes[c]` the flux from each of them into `c`, in the layout of
/// Opm::UpwindGraphBuilder::build() and buildFluxes(). A transport scheme
/// that is upwind in the new time level has this form, with `diagonal` the
/// storage plus the outflow, and `weightFactor` the time step. All quantities
/// share the matrix, so they are solved together.
///
/// The cells are ordered with Opm::tarjan(). Solving the components in that
/// order means every component only needs values that are already final:
///
///  - A cell that is not part of a flow loop is solved with one division.
///  - Cells in a flow loop depend on each other. Small loops are solved
///    directly, large ones with Gauss-Seidel sweeps. The sweeps follow the
///    flow, since the cells of a component come out of Opm::tarjan() with the
///    upstream ones first.
///
/// With non-negative fluxes and a positive diagonal that exceeds the total
/// flux out of each cell, the matrix is an M-matrix. The solution is then
/// non-negative for non-negative `rhs`, and Gauss-Seidel converges.

namespace Opm
{

struct UpwindTransportSettings
{
    /// Cells in a flow loop are solved directly if there are at most this
    /// many of them, and by iteration if there are more.
    int maxDirectSize = 100;

    /// The iteration stops when the largest change of a quantity in a sweep
    /// is at most `tolerance` times its largest value in the component.
    double tolerance = 1e-10;

    /// Most sweeps over a flow loop before giving up.
    int maxIterations = 1000;
};

struct UpwindTransportStatistics
{
    /// Number of components, i.e. cells or flow loops.
    int numComponents = 0;

    /// Number of components with more than one cell, i.e. flow loops.
    int numLoops = 0;

    /// Number of cells in the largest component.
    int largestComponent = 0;

    /// Most sweeps used for any flow loop that was solved by iteration.
    int maxIterations = 0;
};

namespace detail
{

    /// Solve the dense system `a x = b` with `a` of size m x m and `b` of size
    /// m x numRhs, both row-major. The solution replaces `b`, and `a` is destroyed.
    template <class Scalar>
    void solveDenseSystem(std::vector<Scalar>& a,
                          std::vector<Scalar>& b,
                          std::size_t m,
                          std::size_t numRhs)
    {
        for (std::size_t col = 0; col < m; ++col) {
            // Gaussian elimination with partial pivoting.
            std::size_t pivot = col;
            for (std::size_t row = col + 1; row < m; ++row) {
                if (std::abs(a[row * m + col]) > std::abs(a[pivot * m + col])) {
                    pivot = row;
                }
            }
            if (a[pivot * m + col] == Scalar {0}) {
                throw std::runtime_error(
                    "Opm::solveUpwindTransport: singular matrix for a flow loop");
            }
            if (pivot != col) {
                for (std::size_t k = 0; k < m; ++k) {
                    std::swap(a[pivot * m + k], a[col * m + k]);
                }
                for (std::size_t k = 0; k < numRhs; ++k) {
                    std::swap(b[pivot * numRhs + k], b[col * numRhs + k]);
                }
            }
            for (std::size_t row = col + 1; row < m; ++row) {
                const Scalar factor = a[row * m + col] / a[col * m + col];
                if (factor == Scalar {0}) {
                    continue;
                }
                for (std::size_t k = col; k < m; ++k) {
                    a[row * m + k] -= factor * a[col * m + k];
                }
                for (std::size_t k = 0; k < numRhs; ++k) {
                    b[row * numRhs + k] -= factor * b[col * numRhs + k];
                }
            }
        }
        for (std::size_t col = m; col-- > 0;) {
            for (std::size_t k = 0; k < numRhs; ++k) {
                Scalar sum = b[col * numRhs + k];
                for (std::size_t j = col + 1; j < m; ++j) {
                    sum -= a[col * m + j] * b[j * numRhs + k];
                }
                b[col * numRhs + k] = sum / a[col * m + col];
            }
        }
    }


    /// The system and the scratch space for solving it, see solveUpwindTransport().
    template <class Scalar, class Vectors>
    class UpwindTransportSystem
    {
    public:
        UpwindTransportSystem(const SparseTable<int>& upwind,
                              const SparseTable<Scalar>& fluxes,
                              const Scalar weightFactor,
                              const std::vector<Scalar>& diagonal,
                              const Vectors& rhs,
                              Vectors& x,
                              const UpwindTransportSettings& settings)
            : upwind_(upwind)
            , fluxes_(fluxes)
            , weightFactor_(weightFactor)
            , diagonal_(diagonal)
            , rhs_(rhs)
            , x_(x)
            , settings_(settings)
            , numCells_(static_cast<std::size_t>(upwind.size()))
            , numVectors_(rhs.size())
            , localIndex_(numCells_, -1)
            , change_(numVectors_)
            , scale_(numVectors_)
        {
            validate_();
        }

        UpwindTransportStatistics solve()
        {
            const auto components = tarjan(upwind_);
            UpwindTransportStatistics statistics;
            statistics.numComponents = components.numComponents();
            for (const auto cells : components.components()) {
                const int size = static_cast<int>(cells.size());
                statistics.largestComponent = std::max(statistics.largestComponent, size);
                if (size == 1) {
                    solveCell_(cells[0]);
                } else {
                    ++statistics.numLoops;
                    if (size <= settings_.maxDirectSize) {
                        solveDirectly_(cells);
                    } else {
                        const int sweeps = solveByGaussSeidel_(cells);
                        statistics.maxIterations = std::max(statistics.maxIterations, sweeps);
                    }
                }
            }
            return statistics;
        }

    private:
        void validate_() const
        {
            if (fluxes_.size() != upwind_.size() || fluxes_.dataSize() != upwind_.dataSize()
                || diagonal_.size() != numCells_ || x_.size() != numVectors_) {
                throw std::invalid_argument("Opm::solveUpwindTransport: sizes do not match");
            }
            for (std::size_t q = 0; q < numVectors_; ++q) {
                if (static_cast<std::size_t>(rhs_[q].size()) != numCells_
                    || static_cast<std::size_t>(x_[q].size()) != numCells_) {
                    throw std::invalid_argument("Opm::solveUpwindTransport: wrong vector length");
                }
            }
            if (!std::ranges::all_of(diagonal_, [](Scalar d) { return d > Scalar {0}; })) {
                throw std::invalid_argument(
                    "Opm::solveUpwindTransport: the diagonal must be positive");
            }
            for (std::size_t c = 0; c < numCells_; ++c) {
                auto cell = upwind_[static_cast<int>(c)].begin();
                for (const auto flux : fluxes_[static_cast<int>(c)]) {
                    // Written so that NaN is rejected too.
                    if (!(flux >= Scalar {0}) || std::cmp_equal(*cell++, c)) {
                        throw std::invalid_argument(
                            "Opm::solveUpwindTransport: the fluxes must not be negative, "
                            "and a cell cannot be its own upwind cell");
                    }
                }
            }
        }

        /// The right-hand side of cell c for quantity q: `rhs` plus the weighted
        /// values of the upwind cells. The cells of a flow loop that is being
        /// solved directly are marked in `localIndex_`. They are unknowns of the
        /// dense system, so they are left out here. Everywhere else `localIndex_`
        /// is -1, and all upwind cells contribute.
        Scalar rightHandSide_(std::size_t q, int c) const
        {
            Scalar sum = rhs_[q][c];
            auto cell = upwind_[c].begin();
            for (const auto flux : fluxes_[c]) {
                const int u = *cell++;
                if (localIndex_[u] < 0) {
                    sum += weightFactor_ * flux * x_[q][u];
                }
            }
            return sum;
        }

        /// A cell that is not in a flow loop: one division.
        void solveCell_(int c)
        {
            for (std::size_t q = 0; q < numVectors_; ++q) {
                x_[q][c] = rightHandSide_(q, c) / diagonal_[c];
            }
        }

        /// A flow loop with a few cells: solve for all of its cells at once.
        void solveDirectly_(std::span<const int> cells)
        {
            const std::size_t m = cells.size();
            matrix_.assign(m * m, Scalar {0});
            solution_.assign(m * numVectors_, Scalar {0});
            for (std::size_t i = 0; i < m; ++i) {
                localIndex_[cells[i]] = static_cast<int>(i);
            }
            for (std::size_t i = 0; i < m; ++i) {
                const int c = cells[i];
                matrix_[i * m + i] = diagonal_[c];
                auto cell = upwind_[c].begin();
                for (const auto flux : fluxes_[c]) {
                    const int u = *cell++;
                    if (localIndex_[u] >= 0) {
                        matrix_[i * m + static_cast<std::size_t>(localIndex_[u])]
                            -= weightFactor_ * flux;
                    }
                }
                for (std::size_t q = 0; q < numVectors_; ++q) {
                    // Upwind cells outside the loop are already solved.
                    solution_[i * numVectors_ + q] = rightHandSide_(q, c);
                }
            }
            solveDenseSystem(matrix_, solution_, m, numVectors_);
            for (std::size_t i = 0; i < m; ++i) {
                localIndex_[cells[i]] = -1;
                for (std::size_t q = 0; q < numVectors_; ++q) {
                    x_[q][cells[i]] = solution_[i * numVectors_ + q];
                }
            }
        }

        /// A large flow loop: Gauss-Seidel sweeps, in the order of the cells.
        /// Returns the number of sweeps.
        int solveByGaussSeidel_(std::span<const int> cells)
        {
            for (int sweep = 1; sweep <= settings_.maxIterations; ++sweep) {
                std::fill(change_.begin(), change_.end(), Scalar {0});
                std::fill(scale_.begin(), scale_.end(), Scalar {0});
                for (const int c : cells) {
                    for (std::size_t q = 0; q < numVectors_; ++q) {
                        const Scalar updated = rightHandSide_(q, c) / diagonal_[c];
                        change_[q] = std::max(change_[q], std::abs(updated - x_[q][c]));
                        scale_[q] = std::max(scale_[q], std::abs(updated));
                        x_[q][c] = updated;
                    }
                }
                bool converged = true;
                for (std::size_t q = 0; q < numVectors_; ++q) {
                    converged = converged && change_[q] <= settings_.tolerance * scale_[q];
                }
                if (converged) {
                    return sweep;
                }
            }
            throw std::runtime_error("Opm::solveUpwindTransport: no convergence after "
                                     + std::to_string(settings_.maxIterations)
                                     + " sweeps for a flow loop of " + std::to_string(cells.size())
                                     + " cells");
        }

        const SparseTable<int>& upwind_;
        const SparseTable<Scalar>& fluxes_;
        Scalar weightFactor_;
        const std::vector<Scalar>& diagonal_;
        const Vectors& rhs_;
        Vectors& x_;
        const UpwindTransportSettings& settings_;
        std::size_t numCells_;
        std::size_t numVectors_;

        // Scratch space, reused for every flow loop
        std::vector<int> localIndex_;
        std::vector<Scalar> matrix_;
        std::vector<Scalar> solution_;
        std::vector<Scalar> change_;
        std::vector<Scalar> scale_;
    };

} // namespace detail

/// Solve the upwind transport system described in the file comment.
///
/// \param upwind        The cells each cell receives flow from.
/// \param fluxes        The flux from each of those cells. Must have the
///                      layout of `upwind`, and be non-negative.
/// \param weightFactor  Factor for the fluxes in the equations, for example
///                      the time step.
/// \param diagonal      One positive number per cell.
/// \param rhs           One vector of length `numCells` per quantity.
/// \param x             The solution, in the layout of `rhs`. The values on
///                      entry are the initial guess when a flow loop is solved
///                      by iteration; otherwise they are not used.
///
/// \return Statistics about the components that were found and solved.
///
/// \throws std::invalid_argument if the sizes do not match, if a diagonal
///         entry is not positive, if a flux is negative or NaN, or if a
///         cell is its own upwind cell.
/// \throws std::runtime_error if the iteration for a flow loop does not
///         converge, or if a flow loop has a singular matrix.
template <class Scalar, class Vectors>
UpwindTransportStatistics
solveUpwindTransport(const SparseTable<int>& upwind,
                     const SparseTable<Scalar>& fluxes,
                     const Scalar weightFactor,
                     const std::vector<Scalar>& diagonal,
                     const Vectors& rhs,
                     Vectors& x,
                     const UpwindTransportSettings& settings = {})
{
    return detail::UpwindTransportSystem<Scalar, Vectors>(
               upwind, fluxes, weightFactor, diagonal, rhs, x, settings)
        .solve();
}

} // namespace Opm

#endif // OPM_UPWINDTRANSPORTSOLVER_HEADER_INCLUDED
