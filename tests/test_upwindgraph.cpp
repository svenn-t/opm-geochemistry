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
#define BOOST_TEST_MODULE UpwindGraphTests

#include <boost/test/unit_test.hpp>

#include <opm/common/utility/SparseTable.hpp>
#include <opm/grid/CpGrid.hpp>
#include <opm/simulators/utils/Tarjan.hpp>
#include <opm/simulators/utils/UpwindGraph.hpp>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/yaspgrid.hh>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <memory>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

namespace
{

struct MpiFixture
{
    MpiFixture()
    {
        auto& suite = boost::unit_test::framework::master_test_suite();
        Dune::MPIHelper::instance(suite.argc, suite.argv);
    }
};

} // namespace

BOOST_GLOBAL_FIXTURE(MpiFixture);

namespace
{

// Both grids are used with cells of size 1 x 1 whose centres are at
// (0.5, 0.5), (1.5, 0.5), ..., with nx cells in x and ny cells in y.
struct YaspFactory
{
    using Grid = Dune::YaspGrid<2>;

    static std::unique_ptr<Grid> make(int nx, int ny)
    {
        return std::make_unique<Grid>(
            Dune::FieldVector<double, 2> {static_cast<double>(nx), static_cast<double>(ny)},
            std::array<int, 2> {nx, ny});
    }
};

struct CpGridFactory
{
    using Grid = Dune::CpGrid;

    static std::unique_ptr<Grid> make(int nx, int ny)
    {
        auto grid = std::make_unique<Grid>();
        grid->createCartesian({nx, ny, 1}, {1.0, 1.0, 1.0});
        return grid;
    }
};

using Factories = std::tuple<YaspFactory, CpGridFactory>;

using Point = std::array<double, 2>;

struct GraphAndInfo
{
    Opm::SparseTable<int> graph;
    std::vector<Point> centers; // by cell index
    std::vector<std::pair<int, int>> fluxCalls; // arguments of every flux evaluation
};

// Build the upwind graph of a grid when the flux between two neighbouring
// cells is velocity(midpoint) . (centre of b - centre of a).
template <class Grid, class Velocity>
GraphAndInfo
upwindGraphFor(const Grid& grid, Velocity velocity)
{
    using GridView = Grid::LeafGridView;
    const GridView gridView = grid.leafGridView();
    const Dune::MultipleCodimMultipleGeomTypeMapper<GridView> mapper(gridView,
                                                                     Dune::mcmgElementLayout());

    std::vector<Point> centers(mapper.size());
    for (const auto& element : elements(gridView)) {
        const auto center = element.geometry().center();
        centers[mapper.index(element)] = {center[0], center[1]};
    }

    std::vector<std::pair<int, int>> calls;
    const auto flux = [&](int a, int b) {
        calls.emplace_back(a, b);
        const Point& ca = centers[a];
        const Point& cb = centers[b];
        const Point v = velocity(0.5 * (ca[0] + cb[0]), 0.5 * (ca[1] + cb[1]));
        return v[0] * (cb[0] - ca[0]) + v[1] * (cb[1] - ca[1]);
    };

    auto graph = Opm::makeUpwindGraph(gridView, mapper, flux);
    return {std::move(graph), std::move(centers), std::move(calls)};
}

std::vector<int>
rowOf(const Opm::SparseTable<int>& graph, int cell)
{
    const auto row = graph[cell];
    return {row.begin(), row.end()};
}

bool
isAt(const Point& p, double x, double y)
{
    return std::abs(p[0] - x) < 1e-12 && std::abs(p[1] - y) < 1e-12;
}

} // namespace

BOOST_AUTO_TEST_CASE_TEMPLATE(UniformFlowInPlusXDirection, Factory, Factories)
{
    constexpr int nx = 4;
    constexpr int ny = 3;
    const auto grid = Factory::make(nx, ny);
    const auto result = upwindGraphFor(*grid, [](double, double) { return Point {1.0, 0.0}; });
    const auto& graph = result.graph;
    const auto& centers = result.centers;

    // A cell receives from the cell to its left, and the first column
    // receives nothing. No flow in y gives no edges in that direction.
    BOOST_REQUIRE_EQUAL(graph.size(), nx * ny);
    for (int c = 0; c < nx * ny; ++c) {
        const auto row = rowOf(graph, c);
        if (centers[c][0] < 1.0) {
            BOOST_CHECK(row.empty());
        } else {
            BOOST_REQUIRE_EQUAL(row.size(), 1u);
            BOOST_CHECK(isAt(centers[row[0]], centers[c][0] - 1.0, centers[c][1]));
        }
    }

    // Every connection is evaluated exactly once, with the lower index first.
    auto calls = result.fluxCalls;
    BOOST_CHECK_EQUAL(calls.size(), static_cast<std::size_t>((nx - 1) * ny + nx * (ny - 1)));
    BOOST_CHECK(
        std::ranges::all_of(calls, [](const auto& call) { return call.first < call.second; }));
    std::ranges::sort(calls);
    BOOST_CHECK(std::ranges::adjacent_find(calls) == calls.end());

    // No flow loops: every cell is a component of its own, upstream first.
    const auto scc = Opm::tarjan(graph);
    BOOST_REQUIRE_EQUAL(scc.numComponents(), nx * ny);
    for (int c = 0; c < nx * ny; ++c) {
        for (const int upwind : rowOf(graph, c)) {
            BOOST_CHECK_LT(scc.componentOf[upwind], scc.componentOf[c]);
        }
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(UniformFlowInMinusXDirection, Factory, Factories)
{
    constexpr int nx = 4;
    constexpr int ny = 3;
    const auto grid = Factory::make(nx, ny);
    const auto result = upwindGraphFor(*grid, [](double, double) { return Point {-1.0, 0.0}; });

    for (int c = 0; c < nx * ny; ++c) {
        const auto row = rowOf(result.graph, c);
        const Point& centre = result.centers[c];
        if (centre[0] > nx - 1.0) {
            BOOST_CHECK(row.empty());
        } else {
            BOOST_REQUIRE_EQUAL(row.size(), 1u);
            BOOST_CHECK(isAt(result.centers[row[0]], centre[0] + 1.0, centre[1]));
        }
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(NoFlowGivesNoEdges, Factory, Factories)
{
    const auto grid = Factory::make(4, 3);
    const auto result = upwindGraphFor(*grid, [](double, double) { return Point {0.0, 0.0}; });

    BOOST_REQUIRE_EQUAL(result.graph.size(), 12);
    for (int c = 0; c < 12; ++c) {
        BOOST_CHECK(rowOf(result.graph, c).empty());
    }
    BOOST_CHECK_EQUAL(Opm::tarjan(result.graph).numComponents(), 12);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(CirculatingFlowIsOneLoop, Factory, Factories)
{
    // Flow around the centre (1, 1) of a grid with 2 x 2 cells. Each cell
    // receives from exactly one neighbour and delivers to another, so the
    // four cells depend on each other.
    const auto grid = Factory::make(2, 2);
    const auto result
        = upwindGraphFor(*grid, [](double x, double y) { return Point {-(y - 1.0), x - 1.0}; });

    for (int c = 0; c < 4; ++c) {
        BOOST_CHECK_EQUAL(rowOf(result.graph, c).size(), 1u);
    }
    const auto scc = Opm::tarjan(result.graph);
    BOOST_REQUIRE_EQUAL(scc.numComponents(), 1);
    BOOST_CHECK_EQUAL(scc[0].size(), 4u);
}

BOOST_AUTO_TEST_CASE(BuilderGivesTheUpwindCellsOfEachCell)
{
    Opm::UpwindGraphBuilder builder(4);
    builder.addConnection(0, 1, 2.0); // flow 0 -> 1: 0 is upwind of 1
    builder.addConnection(2, 1, 5.0); // flow 2 -> 1: 2 is upwind of 1
    builder.addConnection(1, 3, -1.0); // flow 3 -> 1: 3 is upwind of 1
    builder.addConnection(2, 3, 0.0); // no flow: no edge

    const auto graph = builder.build();
    BOOST_REQUIRE_EQUAL(graph.size(), 4);
    BOOST_CHECK(rowOf(graph, 0).empty());
    BOOST_CHECK(rowOf(graph, 2).empty());
    BOOST_CHECK(rowOf(graph, 3).empty());

    // The edges of a row are in the order the connections were added.
    const std::vector<int> expected {0, 2, 3};
    const auto row = rowOf(graph, 1);
    BOOST_CHECK_EQUAL_COLLECTIONS(row.begin(), row.end(), expected.begin(), expected.end());
}

BOOST_AUTO_TEST_CASE(BuilderGivesTheFluxOfEachEdgeInTheLayoutOfTheGraph)
{
    Opm::UpwindGraphBuilder builder(4);
    builder.addConnection(0, 1, 2.0); // 0 is upwind of 1, flux 2
    builder.addConnection(3, 2, -7.0); // 2 is upwind of 3, flux 7
    builder.addConnection(2, 1, 5.0); // 2 is upwind of 1, flux 5
    builder.addConnection(1, 3, -1.5); // 3 is upwind of 1, flux 1.5
    builder.addConnection(0, 3, 0.0); // no edge, and no flux

    const auto graph = builder.build();
    const auto fluxes = builder.buildFluxes();
    BOOST_REQUIRE_EQUAL(fluxes.size(), graph.size());
    for (int c = 0; c < graph.size(); ++c) {
        BOOST_REQUIRE_EQUAL(fluxes[c].size(), graph[c].size());
    }

    const auto fluxRow = [&](int c) {
        const auto row = fluxes[c];
        return std::vector<double> {row.begin(), row.end()};
    };
    BOOST_CHECK((rowOf(graph, 1) == std::vector<int> {0, 2, 3}));
    BOOST_CHECK((fluxRow(1) == std::vector<double> {2.0, 5.0, 1.5}));
    BOOST_CHECK((rowOf(graph, 3) == std::vector<int> {2}));
    BOOST_CHECK((fluxRow(3) == std::vector<double> {7.0}));
    BOOST_CHECK(fluxRow(0).empty());
    BOOST_CHECK(fluxRow(2).empty());
}

BOOST_AUTO_TEST_CASE(BuilderWithoutCells)
{
    const Opm::UpwindGraphBuilder builder(0);
    BOOST_CHECK_EQUAL(builder.build().size(), 0);
}

BOOST_AUTO_TEST_CASE(BuilderRejectsNaNFlux)
{
    Opm::UpwindGraphBuilder builder(3);
    BOOST_CHECK_THROW(builder.addConnection(0, 1, std::numeric_limits<double>::quiet_NaN()),
                      std::invalid_argument);
    // The rejected connection must not have left an edge behind.
    const auto graph = builder.build();
    BOOST_CHECK(rowOf(graph, 0).empty());
    BOOST_CHECK(rowOf(graph, 1).empty());
}

BOOST_AUTO_TEST_CASE(BuilderRejectsCellsOutsideTheGrid)
{
    Opm::UpwindGraphBuilder builder(3);
    BOOST_CHECK_THROW(builder.addConnection(0, 3, 1.0), std::out_of_range);
    BOOST_CHECK_THROW(builder.addConnection(-1, 1, 1.0), std::out_of_range);
    BOOST_CHECK_EQUAL(builder.build().size(), 3);
}
