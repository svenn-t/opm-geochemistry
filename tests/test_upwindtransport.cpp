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
#define BOOST_TEST_MODULE UpwindTransportTests

#include <boost/test/unit_test.hpp>

#include <opm/common/utility/SparseTable.hpp>
#include <opm/simulators/utils/UpwindGraph.hpp>
#include <opm/simulators/utils/UpwindTransportSolver.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace
{

using Vectors = std::vector<std::vector<double>>;

// A flow of `flux` from cell `from` into cell `to`.
struct Edge
{
    int from;
    int to;
    double flux;
};

struct System
{
    Opm::SparseTable<int> upwind;
    Opm::SparseTable<double> fluxes;
    std::vector<double> diagonal;
};

// The system of an implicit upwind scheme: the diagonal is the volume of each
// cell plus `factor` times the flux out of it.
System
makeSystem(const std::vector<double>& volumes, const std::vector<Edge>& edges, double factor)
{
    Opm::UpwindGraphBuilder builder(static_cast<int>(volumes.size()));
    std::vector<double> diagonal = volumes;
    for (const Edge& edge : edges) {
        builder.addConnection(edge.from, edge.to, edge.flux);
        diagonal[static_cast<std::size_t>(edge.from)] += factor * edge.flux;
    }
    return {builder.build(), builder.buildFluxes(), diagonal};
}

// Dense solve of the whole system, by Gaussian elimination with partial pivoting.
Vectors
referenceSolution(const std::vector<double>& diagonal,
                  const std::vector<Edge>& edges,
                  double factor,
                  const Vectors& rhs)
{
    const std::size_t n = diagonal.size();
    std::vector<std::vector<double>> a(n, std::vector<double>(n, 0.0));
    for (std::size_t c = 0; c < n; ++c) {
        a[c][c] = diagonal[c];
    }
    for (const Edge& edge : edges) {
        a[static_cast<std::size_t>(edge.to)][static_cast<std::size_t>(edge.from)]
            -= factor * edge.flux;
    }

    Vectors result;
    for (const auto& b : rhs) {
        auto m = a;
        auto x = b;
        for (std::size_t col = 0; col < n; ++col) {
            std::size_t pivot = col;
            for (std::size_t row = col + 1; row < n; ++row) {
                if (std::abs(m[row][col]) > std::abs(m[pivot][col])) {
                    pivot = row;
                }
            }
            std::swap(m[pivot], m[col]);
            std::swap(x[pivot], x[col]);
            for (std::size_t row = col + 1; row < n; ++row) {
                const double f = m[row][col] / m[col][col];
                for (std::size_t k = col; k < n; ++k) {
                    m[row][k] -= f * m[col][k];
                }
                x[row] -= f * x[col];
            }
        }
        for (std::size_t col = n; col-- > 0;) {
            for (std::size_t j = col + 1; j < n; ++j) {
                x[col] -= m[col][j] * x[j];
            }
            x[col] /= m[col][col];
        }
        result.push_back(x);
    }
    return result;
}

// Largest difference to the reference, or a negative number if a value is negative.
double
maxDifference(const Vectors& x, const Vectors& reference)
{
    double largest = 0.0;
    for (std::size_t q = 0; q < x.size(); ++q) {
        for (std::size_t c = 0; c < x[q].size(); ++c) {
            if (!(x[q][c] >= -1e-12)) {
                return -1.0;
            }
            largest = std::max(largest, std::abs(x[q][c] - reference[q][c]));
        }
    }
    return largest;
}

Opm::UpwindTransportSettings
alwaysDirect()
{
    Opm::UpwindTransportSettings settings;
    settings.maxDirectSize = 1000;
    return settings;
}

Opm::UpwindTransportSettings
alwaysIterate()
{
    Opm::UpwindTransportSettings settings;
    settings.maxDirectSize = 0;
    settings.tolerance = 1e-14;
    settings.maxIterations = 100000;
    return settings;
}

} // namespace

BOOST_AUTO_TEST_CASE(ChainOfCellsIsSolvedOneByOne)
{
    // 0 -> 1 -> 2 with fluxes 2 and 3 and storage 4, 5, 6 minus the outflow.
    const auto system = makeSystem({2.0, 2.0, 6.0}, {{0, 1, 2.0}, {1, 2, 3.0}}, 1.0);
    BOOST_REQUIRE_EQUAL(system.diagonal[0], 4.0);
    BOOST_REQUIRE_EQUAL(system.diagonal[1], 5.0);
    BOOST_REQUIRE_EQUAL(system.diagonal[2], 6.0);

    const Vectors rhs {{8.0, 10.0, 12.0}};
    Vectors x(1, std::vector<double>(3, 0.0));
    const auto statistics
        = Opm::solveUpwindTransport(system.upwind, system.fluxes, 1.0, system.diagonal, rhs, x);

    BOOST_CHECK_CLOSE(x[0][0], 2.0, 1e-12); // 8 / 4
    BOOST_CHECK_CLOSE(x[0][1], 2.8, 1e-12); // (10 + 2 * 2) / 5
    BOOST_CHECK_CLOSE(x[0][2], 3.4, 1e-12); // (12 + 3 * 2.8) / 6
    BOOST_CHECK_EQUAL(statistics.numComponents, 3);
    BOOST_CHECK_EQUAL(statistics.numLoops, 0);
    BOOST_CHECK_EQUAL(statistics.largestComponent, 1);
    BOOST_CHECK_EQUAL(statistics.maxIterations, 0);
}

BOOST_AUTO_TEST_CASE(TwoCellsThatFeedEachOtherAreSolvedTogether)
{
    // 3 x0 = 1 + 2 x1 and 4 x1 = 2 + 1 x0, so x1 = 0.7 and x0 = 0.8.
    const auto system = makeSystem({1.0, 1.0}, {{0, 1, 1.0}, {1, 0, 2.0}}, 1.0);
    const std::vector<double> diagonal {3.0, 4.0};
    const Vectors rhs {{1.0, 2.0}};
    for (const auto& settings : {alwaysDirect(), alwaysIterate()}) {
        Vectors x(1, std::vector<double>(2, 0.0));
        const auto statistics = Opm::solveUpwindTransport(
            system.upwind, system.fluxes, 1.0, diagonal, rhs, x, settings);
        BOOST_CHECK_CLOSE(x[0][0], 0.8, 1e-9);
        BOOST_CHECK_CLOSE(x[0][1], 0.7, 1e-9);
        BOOST_CHECK_EQUAL(statistics.numComponents, 1);
        BOOST_CHECK_EQUAL(statistics.numLoops, 1);
        BOOST_CHECK_EQUAL(statistics.largestComponent, 2);
        BOOST_CHECK_EQUAL(statistics.maxIterations > 0, settings.maxDirectSize == 0);
    }
}

BOOST_AUTO_TEST_CASE(ChainAndLoopTogether)
{
    // Cells 0 -> 1 -> 2 in a chain, and 3 and 4 feeding each other.
    const std::vector<Edge> edges {{0, 1, 1.0}, {1, 2, 1.0}, {3, 4, 1.0}, {4, 3, 1.0}};
    const auto system = makeSystem({1.0, 1.0, 1.0, 1.0, 1.0}, edges, 1.0);
    const Vectors rhs {{1.0, 1.0, 1.0, 1.0, 1.0}};
    Vectors x(1, std::vector<double>(5, 0.0));
    const auto statistics
        = Opm::solveUpwindTransport(system.upwind, system.fluxes, 1.0, system.diagonal, rhs, x);
    BOOST_CHECK_EQUAL(statistics.numComponents, 4);
    BOOST_CHECK_EQUAL(statistics.numLoops, 1);
    BOOST_CHECK_EQUAL(statistics.largestComponent, 2);
    // The solution agrees with a solve of the whole system.
    BOOST_CHECK_LT(maxDifference(x, referenceSolution(system.diagonal, edges, 1.0, rhs)), 1e-12);
}

BOOST_AUTO_TEST_CASE(RandomSystemsAgreeWithADenseSolve)
{
    // A fixed seed is intended: the test must give the same systems on every
    // run, so that a failure can be reproduced.
    std::mt19937 rng(20261009);
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    int numLoops = 0;
    int largestLoop = 0;

    for (int trial = 0; trial < 200; ++trial) {
        const int n = 1 + static_cast<int>(rng() % 30);
        const double factor = 0.1 + 3.0 * unit(rng);

        std::vector<double> volumes(n);
        for (double& volume : volumes) {
            volume = 0.2 + 1.8 * unit(rng);
        }
        std::vector<Edge> edges;
        const int numEdges = static_cast<int>(rng() % (2 * n + 1));
        for (int e = 0; e < numEdges; ++e) {
            const int from = static_cast<int>(rng() % n);
            const int to = static_cast<int>(rng() % n);
            if (from != to) {
                edges.push_back({from, to, 0.1 + 1.9 * unit(rng)});
            }
        }

        Vectors rhs(3, std::vector<double>(n));
        for (auto& b : rhs) {
            for (double& value : b) {
                value = unit(rng);
            }
        }
        const auto system = makeSystem(volumes, edges, factor);
        const auto reference = referenceSolution(system.diagonal, edges, factor, rhs);

        for (const auto& settings : {alwaysDirect(), alwaysIterate()}) {
            Vectors x(rhs.size(), std::vector<double>(n, 0.0));
            const auto statistics = Opm::solveUpwindTransport(
                system.upwind, system.fluxes, factor, system.diagonal, rhs, x, settings);
            numLoops += statistics.numLoops;
            largestLoop = std::max(largestLoop, statistics.largestComponent);
            const double difference = maxDifference(x, reference);
            BOOST_REQUIRE_MESSAGE(difference >= 0.0 && difference < 1e-8,
                                  "trial " + std::to_string(trial) + ", "
                                      + (settings.maxDirectSize == 0 ? "iterative" : "direct")
                                      + ": difference " + std::to_string(difference));
        }
    }
    // The test is only meaningful if there were flow loops, small and large.
    BOOST_CHECK_GT(numLoops, 100);
    BOOST_CHECK_GT(largestLoop, 5);
}

BOOST_AUTO_TEST_CASE(IterationThatDoesNotConvergeThrows)
{
    // Two cells that pass almost all their content to each other: the
    // Gauss-Seidel iteration needs far more than 3 sweeps.
    const std::vector<Edge> edges {{0, 1, 1.0}, {1, 0, 1.0}};
    const auto system = makeSystem({1e-3, 1e-3}, edges, 1.0);
    const Vectors rhs {{1.0, 0.0}};
    Vectors x(1, std::vector<double>(2, 0.0));

    auto settings = alwaysIterate();
    settings.maxIterations = 3;
    BOOST_CHECK_THROW(Opm::solveUpwindTransport(
                          system.upwind, system.fluxes, 1.0, system.diagonal, rhs, x, settings),
                      std::runtime_error);
}

BOOST_AUTO_TEST_CASE(InvalidInputIsRejected)
{
    const std::vector<double> volumes {1.0, 1.0};
    const Vectors rhs {{1.0, 1.0}};
    Vectors x(1, std::vector<double>(2, 0.0));

    const auto good = makeSystem(volumes, {{0, 1, 1.0}}, 1.0);
    BOOST_CHECK_NO_THROW(
        Opm::solveUpwindTransport(good.upwind, good.fluxes, 1.0, good.diagonal, rhs, x));

    // A diagonal that is not positive.
    const std::vector<double> zeroDiagonal {1.0, 0.0};
    BOOST_CHECK_THROW(
        Opm::solveUpwindTransport(good.upwind, good.fluxes, 1.0, zeroDiagonal, rhs, x),
        std::invalid_argument);

    // A cell that is its own upwind cell.
    const auto selfEdge = makeSystem(volumes, {{1, 1, 1.0}}, 1.0);
    BOOST_CHECK_THROW(
        Opm::solveUpwindTransport(selfEdge.upwind, selfEdge.fluxes, 1.0, selfEdge.diagonal, rhs, x),
        std::invalid_argument);

    // Vectors of the wrong length.
    Vectors shortX(1, std::vector<double>(1, 0.0));
    BOOST_CHECK_THROW(
        Opm::solveUpwindTransport(good.upwind, good.fluxes, 1.0, good.diagonal, rhs, shortX),
        std::invalid_argument);

    // Fluxes in another layout than the graph.
    const auto other = makeSystem(volumes, {}, 1.0);
    BOOST_CHECK_THROW(
        Opm::solveUpwindTransport(good.upwind, other.fluxes, 1.0, good.diagonal, rhs, x),
        std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(LongChainDoesNotOverflowTheStack)
{
    // x[0] = 1/2 and x[k] = (1 + x[k - 1]) / 2, so x[k] = 1 - 2^-(k+1).
    constexpr int n = 100000;
    const std::vector<double> volumes(n, 1.0);
    std::vector<Edge> edges;
    edges.reserve(n - 1);
    for (int c = 0; c + 1 < n; ++c) {
        edges.push_back({c, c + 1, 1.0});
    }
    auto system = makeSystem(volumes, edges, 1.0);
    // The diagonal is 2 for all but the last cell, which has no outflow.
    system.diagonal[n - 1] = 2.0;

    const Vectors rhs {std::vector<double>(n, 1.0)};
    Vectors x(1, std::vector<double>(n, 0.0));
    const auto statistics
        = Opm::solveUpwindTransport(system.upwind, system.fluxes, 1.0, system.diagonal, rhs, x);

    BOOST_CHECK_EQUAL(statistics.numComponents, n);
    BOOST_CHECK_CLOSE(x[0][0], 0.5, 1e-12);
    BOOST_CHECK_CLOSE(x[0][9], 1.0 - std::pow(0.5, 10), 1e-10);
    BOOST_CHECK_CLOSE(x[0][n - 1], 1.0, 1e-10);
}
