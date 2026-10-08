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

#ifndef OPM_UPWINDGRAPH_HEADER_INCLUDED
#define OPM_UPWINDGRAPH_HEADER_INCLUDED

#include <opm/common/utility/SparseTable.hpp>

#include <dune/grid/common/rangegenerators.hh>

#include <cmath>
#include <cstddef>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

/// \file
///
/// The upwind graph of a flow field, as input to Opm::tarjan() (see
/// Tarjan.hpp). Row `c` of the graph lists the cells that cell `c` receives
/// flow from, i.e. there is an edge from `c` to each of its upwind
/// neighbours. Opm::tarjan() then numbers the strongly connected components
/// so that every component only depends on components with a lower number.

namespace Opm
{

/// Collects connections between cells together with the flux through them
/// and turns them into an upwind graph.
class UpwindGraphBuilder
{
public:
    /// \param numCells Number of cells, which are numbered 0, ..., numCells - 1.
    explicit UpwindGraphBuilder(int numCells)
        : numCells_(numCells)
    {
    }

    /// Record the connection between cells `a` and `b`.
    ///
    /// \param flux Flux from `a` to `b`. A positive flux makes `a` upwind of
    ///             `b`, a negative flux makes `b` upwind of `a`. Zero flux
    ///             gives no edge, since neither cell depends on the other.
    ///             Only the sign is used.
    ///
    /// Adding the same connection more than once gives duplicate edges,
    /// which Opm::tarjan() accepts.
    ///
    /// \throws std::out_of_range if a cell is not in [0, numCells).
    /// \throws std::invalid_argument if the flux is NaN, which means the
    ///         flux computation has gone wrong and the direction of the
    ///         connection is unknown.
    void addConnection(int a, int b, double flux)
    {
        checkCell(a);
        checkCell(b);
        if (std::isnan(flux)) {
            throw std::invalid_argument("Opm::UpwindGraphBuilder: NaN flux between cells "
                                        + std::to_string(a) + " and " + std::to_string(b));
        }
        if (flux > 0.0) {
            upwindOf_.emplace_back(b, a);
        } else if (flux < 0.0) {
            upwindOf_.emplace_back(a, b);
        }
    }

    /// The graph with one row per cell. The edges of a row are in the order
    /// the connections were added.
    SparseTable<int> build() const
    {
        // Counting sort of the edges by the cell they belong to.
        std::vector<int> rowStart(static_cast<std::size_t>(numCells_) + 1, 0);
        for (const auto& edge : upwindOf_) {
            ++rowStart[static_cast<std::size_t>(edge.first) + 1];
        }
        std::partial_sum(rowStart.begin(), rowStart.end(), rowStart.begin());

        std::vector<int> next(rowStart.begin(), rowStart.end() - 1);
        std::vector<int> upwindCells(upwindOf_.size());
        for (const auto& [cell, upwindCell] : upwindOf_) {
            upwindCells[static_cast<std::size_t>(next[cell]++)] = upwindCell;
        }
        return {std::move(upwindCells), std::move(rowStart)};
    }

private:
    void checkCell(int cell) const
    {
        if (cell < 0 || cell >= numCells_) {
            throw std::out_of_range("Opm::UpwindGraphBuilder: cell " + std::to_string(cell)
                                    + " is not in a grid with " + std::to_string(numCells_)
                                    + " cells");
        }
    }

    int numCells_;
    std::vector<std::pair<int, int>> upwindOf_;
};

/// Build the upwind graph of a grid from the flux between neighbouring cells.
///
/// Only the generic Dune grid interface is used (`elements`, `intersections`,
/// `neighbor()`, `outside()`), so any grid view works, for example the leaf
/// view of a `Dune::CpGrid`.
///
/// \param gridView The grid view to take the cells and their neighbours from.
/// \param mapper   An element mapper for the view. Its indices number the
///                 cells, and the rows of the graph.
/// \param flux     A function `flux(a, b)` for the flux from cell `a` to cell
///                 `b`, positive when the flow goes from `a` to `b`. It is
///                 called once for every connection between two cells, with
///                 `a < b`, so the flux of a connection is only computed once
///                 and both cells agree on its direction. If two cells share
///                 more than one intersection, it is called once for each.
///
/// Intersections without a neighbour, i.e. the boundary of the grid, give no
/// edge, so flow through them has to be treated as a source term by the
/// caller. All elements of the view are included, also overlap cells in a
/// parallel run; the caller decides how to treat those.
template <class GridView, class ElementMapper, class FluxFunction>
SparseTable<int>
makeUpwindGraph(const GridView& gridView, const ElementMapper& mapper, FluxFunction&& flux)
{
    UpwindGraphBuilder builder(static_cast<int>(mapper.size()));
    for (const auto& element : elements(gridView)) {
        const int inside = static_cast<int>(mapper.index(element));
        for (const auto& intersection : intersections(gridView, element)) {
            if (!intersection.neighbor()) {
                continue;
            }
            const int outside = static_cast<int>(mapper.index(intersection.outside()));
            // Each connection is seen from both sides; use the lower index.
            if (outside > inside) {
                builder.addConnection(inside, outside, flux(inside, outside));
            }
        }
    }
    return builder.build();
}

} // namespace Opm

#endif // OPM_UPWINDGRAPH_HEADER_INCLUDED
