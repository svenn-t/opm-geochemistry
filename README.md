# opm-geochemistry
Geochemical module tailored to porous media applications.

## Contains
The opm-geochemistry module includes:
- Interface to OPM Flow simulator
- Standalone geochemical equilibrium solver
- Standalone 1D reactive transport solver

## Build instructions
opm-geochemistry depend on:

- [opm-common](https://github.com/OPM/opm-common)
- [opm-grid](https://github.com/OPM/opm-grid)
- [opm-simulators](https://github.com/OPM/opm-simulators)

Follow the [build instructions on the OPM home page](https://opm-project.org/?page_id=36).

After the prerequisits are in place, build the opm-geochemistry module by running:
```shell
mkdir build/
cd build/
cmake ..
make
```

## Code documentation
To build the in-code documentation, [Doxygen](https://www.doxygen.nl/index.html) must be installed. Build documentation
with `make doc`. In the build folder, the code documentation will be located at `docs/doxygen` in html, latex and xml
format.

## OPM Flow interface
The executables starting with the name `flow_` are OPM Flow simulators with interface to geochemistry solver.

See [tutorials section](#tutorials) for usage.

## GeoChemX
The executable `GeoChemX` has two sub-programs `EQSOLVER` and `TRANSPORT` which invokes the equilibrium and 1D reactive
transport solvers, respectively. Both sub-programs require an input file as the only command-line argument. Example of
usage for both programs:
```shell
GeoChemX EQSOLVER <input-file>
```
```shell
GeoChemX TRANSPORT <input-file>
```

See [tutorials section](#tutorials) for usage.

## Python
A Python wrapper for `GeoChemX` can be found in [pyGeoChem](pyGeoChem). Follow this [README](pyGeoChem/README.md) to
install. Note that `GeoChemX` must be installed separately using the [build instructions](#build-instructions) before
using `pyGeoChem`.

An experimental Python binding to `EQSOLVER` using [pybind11](https://pybind11.readthedocs.io/en/stable/) is located
[here](python).

## Tutorials
Tutorials can be found here:
- [HTML version](https://opm.github.io/opm-geochemistry/)
- [Notebook/Markdown version](./docs/tutorials/README.md)

## Optional dependencies
For the `TRANSPORT` solver, it is possible to write time series files with the
[NetCDF](https://www.unidata.ucar.edu/software/netcdf/) software. NetCDF files can be read using, e.g., the
[xarray](https://docs.xarray.dev/en/stable/) Python package.

## Database
The thermodynamic database that is used in this module is the `slop07.dat` database [[1]](#1), which is distributed
under the [CC BY 4.0](database/LICENSE_slop07.txt) license.

## Tests
Unit tests for the core geochemistry code is located in [tests/tests_using_catch2](./tests/tests_using_catch2/), which
use [Catch2 v2](https://github.com/catchorg/Catch2) and [ApprovalTests](https://approvaltests.com/). To build and run
tests, use [CTest](https://cmake.org/cmake/help/book/mastering-cmake/chapter/Testing%20With%20CMake%20and%20CTest.html)
facilities (standard in most IDEs). Alternatively, set CMake option `ENABLE_TESTS=ON`, build `make catch2_tests`,
navigate to the `catch2_run_dir/` folder, and run `catch2_tests` executable.

## Papers using `GeoChemX`
- A. Hiorth, L. Cathles and M. Madland. The Impact of Pore Water Chemistry on Carbonate Surface Charge and Oil
  Wettability, Transport in Porous Media, 85, pp. 1-21, 2010,
  [10.1007/s11242-010-9543-6](https://doi.org/10.1007/s11242-010-9543-6).
- A. Hiorth, E. Jettestuen, L. Cathles and M. Madland. Precipitation, Dissolution, and Ion Exchange Processes Coupled
  With a Lattice Boltzmann Advection Diffusion Solver, Geochimica et Cosmochimica Acta, 104, pp. 99-110, 2013,
  [10.1016/j.gca.2012.11.019](https://doi.org/10.1016/j.gca.2012.11.019).
- F. Feldmann, O. Nodland, J. Sagen, B. Antonsen, T. Sira, J. L. Vinningland, R. Moe and A. Hiorth. IORSim: a
  Mathematical Workflow for Field-Scale Geochemistry Simulations in Porous Media, Transport in Porous Media, 151(9), pp.
  1781-1809, 2024, [10.1007/s11242-024-02094-9](https://doi.org/10.1007/s11242-024-02094-9).
- O. Nodland and A. Hiorth. A New Formulation of the Surface Charge/surface Potential Relationship in Electrolytes With
  Valence Less Than Three, Computational Geosciences, 28(2), pp. 289-304, 2024,
  [10.1007/s10596-023-10239-w](https://doi.org/10.1007/s10596-023-10239-w).

## License
The core geochemistry code is distributed under the MIT license, see [MIT LICENSE](LICENSE.MIT). OPM interface code is
distributed under the GPL version 3+ license, see [GPLv3+ LICENSE](LICENSE.GPL-V.3).

## Contributors
The core geochemistry solver was developed by:

- [Aksel Hiorth](https://github.com/ahiorth/)
- [Oddbjørn Nødland](https://github.com/onoedland/)
- [Espen Jettestuen](https://github.com/eje74)

## References
<a id="1">[1]</a>
GEOPIG. (2019). GEOPIG slop files [Data set]. Zenodo.
[10.5281/zenodo.2630819](https://doi.org/10.5281/zenodo.2630820)
