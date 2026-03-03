# -*- mode: cmake; tab-width: 2; indent-tabs-mode: t; truncate-lines: t; compile-command: "cmake -Wdev" -*-
# vim: set filetype=cmake autoindent tabstop=2 shiftwidth=2 noexpandtab softtabstop=2 nowrap:

###########################################################################
#                                                                         #
# Note: The bulk of the build system is located in the cmake/ directory.  #
#       This file only contains the specializations for this particular   #
#       project. Most likely you are interested in editing one of these   #
#       files instead:                                                    #
#                                                                         #
#       dune.module                Name and version number                #
#       CMakeLists_files.cmake     Path of source files                   #
#       ${project}-prereqs.cmake   Dependencies                           #
#                                                                         #
###########################################################################

# Minimum CMake version
cmake_minimum_required(VERSION 3.23)

# Mandatory call to project
project(opm-geochemistry C CXX)

# Options for OPM build
option(SIBLING_SEARCH "Search for other modules in sibling directories?" ON)
option(USE_TRACY_PROFILER "Enable tracy profiling" OFF)
option(USE_HYPRE "Use the Hypre library for linear solvers?" OFF)

# Other options
option(ENABLE_TESTS "Build geochemistry tests?" OFF)
option(BUILD_PYTHON "Attempt to build Python package?" OFF)

# Search for opm-common (copied from opm-simulators CMakeLists.txt)
if(SIBLING_SEARCH AND NOT opm-common_DIR)
  # guess the sibling dir
  get_filename_component(_leaf_dir_name ${PROJECT_BINARY_DIR} NAME)
  get_filename_component(_parent_full_dir ${PROJECT_BINARY_DIR} DIRECTORY)
  get_filename_component(_parent_dir_name ${_parent_full_dir} NAME)
  #Try if <module-name>/<build-dir> is used
  get_filename_component(_modules_dir ${_parent_full_dir} DIRECTORY)
  if(IS_DIRECTORY ${_modules_dir}/opm-common/${_leaf_dir_name})
    set(opm-common_DIR ${_modules_dir}/opm-common/${_leaf_dir_name})
  else()
    string(REPLACE ${PROJECT_NAME} opm-common _opm_common_leaf ${_leaf_dir_name})
    if(NOT _leaf_dir_name STREQUAL _opm_common_leaf
        AND IS_DIRECTORY ${_parent_full_dir}/${_opm_common_leaf})
      # We are using build directories named <prefix><module-name><postfix>
      set(opm-common_DIR ${_parent_full_dir}/${_opm_common_leaf})
    elseif(IS_DIRECTORY ${_parent_full_dir}/opm-common)
      # All modules are in a common build dir
      set(opm-common_DIR "${_parent_full_dir}/opm-common")
    endif()
  endif()
endif()
if(opm-common_DIR AND NOT IS_DIRECTORY ${opm-common_DIR})
  message(WARNING "Value ${opm-common_DIR} passed to variable"
    " opm-common_DIR is not a directory")
endif()

find_package(opm-common REQUIRED)

# OpmInit and policies
include(OpmInit)
OpmSetPolicies()

# MPI and Hypre linear solvers
if(USE_MPI)
	set(HDF5_PREFER_PARALLEL TRUE)
endif()

if(USE_HYPRE AND USE_MPI)
  find_package(HYPRE)
  if(HYPRE_FOUND)
    set(HAVE_HYPRE 1)
  else()
    message(WARNING "Hypre requested but not found. Continuing without Hypre support.")
    set(USE_HYPRE OFF)
  endif()
elseif(USE_HYPRE)
  message(WARNING "Hypre requested but MPI not activated. Continuing without Hypre support.")
  set(USE_HYPRE OFF)
  set(HYPRE_FOUND OFF)
endif()

# not the same location as most of the other projects? this hook overrides
macro (dir_hook)
endmacro (dir_hook)

# Look for the opm-tests repository; if found the variable
# HAVE_OPM_TESTS will be set to true.
include(Findopm-tests)

# list of prerequisites for this particular project; this is in a
# separate file (in cmake/Modules sub-directory) because it is shared
# with the find module
include ("${project}-prereqs")

# read the list of components from this file (in the project
# directory); it should set various lists with the names of the files
# to include
include (CMakeLists_files.cmake)

# opm-geochemistry specific dependencies
find_package(ApprovalTests CONFIG)
if (NOT ApprovalTests_FOUND)
  include(${PROJECT_SOURCE_DIR}/cmake/Modules/FindApprovalTests.cmake)
endif()
find_package(Catch2 CONFIG)
if (NOT Catch2_FOUND)
  include(${PROJECT_SOURCE_DIR}/cmake/Modules/FindCatch2.cmake)
endif()
find_package(cereal CONFIG)
if (NOT cereal_FOUND)
  include(${PROJECT_SOURCE_DIR}/cmake/Modules/Findcereal.cmake)
endif()
find_package(nlohmann_json CONFIG)
if (NOT nlohmann_json_FOUND)
  include(${PROJECT_SOURCE_DIR}/cmake/Modules/Findnlohmann_json.cmake)
endif()


# Some macro definitions
macro (config_hook)
	opm_need_version_of ("dune-common")
	opm_need_version_of ("dune-istl")
	if(dune-fem_FOUND)
		opm_need_version_of ("dune-fem")
	endif()
	opm_need_version_of ("opm-simulators")

	if(USE_TRACY_PROFILER AND Tracy_FOUND)
		set(USE_TRACY 1)
	endif()

	include(UseDamaris)
endmacro (config_hook)

macro (prereqs_hook)
endmacro (prereqs_hook)

macro (fortran_hook)
endmacro (fortran_hook)

macro (files_hook)
	if(MPI_FOUND AND HDF5_FOUND AND NOT HDF5_IS_PARALLEL)
		message(WARNING "When building parallel OPM flow we need a "
			"parallel version of hdf5, but found only a serial one. "
			"Please install a parallel hdf5 library for MPI "
			"(e.g with apt-get install libhdf5-mpi-dev) and do a clean "
			"rebuild (build after \"make clean\"). Continuing with "
			"only normal restart without hdf5 file support.")
		set(HDF5_FOUND OFF)
		unset(HAVE_HDF5)
	endif()
	include (CMakeLists_files.cmake)
endmacro (files_hook)

macro (tests_hook)
endmacro (tests_hook)

# all setup common to the OPM library modules is done here
include (OpmLibMain)

# command line completion for bash
include(OpmBashCompletion)

# HACK: fix for cjson ld errors
target_link_libraries(opmgeochemistry PRIVATE opmcommon cjson)

# Generate flow_blackoil_geochemistry
add_library(moduleVersionGeochem OBJECT opm/simulators/utils/moduleVersion.cpp)
set_property(TARGET moduleVersionGeochem PROPERTY POSITION_INDEPENDENT_CODE ON)
add_dependencies(moduleVersionGeochem opmgeochemistry)
opm_add_target_options(TARGET moduleVersionGeochem)

opm_add_test(flow_blackoil_geochemistry
             ONLY_COMPILE
             ALWAYS_ENABLE
             DEPENDS opmsimulators opmgeochemistry
             LIBRARIES opmsimulators opmgeochemistry
             SOURCES
             flow/flow_blackoil_geochemistry.cpp
             $<TARGET_OBJECTS:moduleVersionGeochem>
)

opm_add_test(flow_onephase_geochemistry
             ONLY_COMPILE
             ALWAYS_ENABLE
             DEPENDS opmsimulators opmgeochemistry
             LIBRARIES opmsimulators opmgeochemistry
             SOURCES
             flow/flow_onephase_geochemistry.cpp
             $<TARGET_OBJECTS:moduleVersionGeochem>
)

# Manually add OPM tests
opm_add_test(test_opmgeocheminterface
             DEPENDS opmsimulators opmgeochemistry
             LIBRARIES opmsimulators opmgeochemistry Boost::unit_test_framework
             SOURCES
             tests/test_opmgeocheminterface.cpp
             $<TARGET_OBJECTS:moduleVersionGeochem>
)

# Build GeoChemX application
add_executable(GeoChemX
  examples/GeoChemX.cpp
)
target_link_libraries(GeoChemX
  PRIVATE
    opmgeochemistry
)

# Build ThermoX application
add_executable(ThermoX
  examples/ThermoX.cpp
)
target_link_libraries(ThermoX
  PRIVATE
    opmgeochemistry
)

# Build tests
if(ENABLE_TESTS)
  add_subdirectory(tests)
endif()

if( NOT "${PYTHON_SEARCH_DIR}" STREQUAL "")
    message("Searching for Python interpreter in directory: ${PYTHON_SEARCH_DIR}")
    set (Python_ROOT_DIR "${PYTHON_SEARCH_DIR}")
endif()

# Build Python bindings
if(BUILD_PYTHON)
    # Check if Python is installed
    find_package (Python COMPONENTS Interpreter)
    if(Python_FOUND)
        message("Python interpreter was found at: \"${Python_EXECUTABLE}\"")
    else()
        message ("Python interpreter was not found, cannot create bindings to C++ project...")
        set(BUILD_PYTHON OFF)
    endif()

endif()

# Continue with Python build if interpreter was found
if(BUILD_PYTHON)
    # Check if pybind11 is installed, and if not download from github
    find_package(pybind11 CONFIG)
    if (NOT pybind11_FOUND)
        include(DownloadPyBind11)
    endif()

    # The relative path where site-packages directory is located for current
    # platform. Note that this is the "local" version of the Python package
    # (stored inside the CMake build directory)
    execute_process(
        COMMAND "${PYTHON_EXECUTABLE}" -c "if True:
            from distutils import sysconfig as sc
            print(sc.get_python_lib(prefix='', plat_specific=True))"
        OUTPUT_VARIABLE PYTHON_SITE_PACKAGES_DIR
        OUTPUT_STRIP_TRAILING_WHITESPACE)
    message("OUTPUT_VARIABLE PYTHON_SITE_PACKAGES_DIR: \"${PYTHON_SITE_PACKAGES_DIR}\"")

    # The directory where the Python package is built
    #
    # Note: This needs improvement. The complication below is due to the need of
    # pybind11-stubgen in importing the package, and this is done before running
    # python setup.py. Once the package is built with setuptools,
    # the directory PYTHON_PACKAGE_DIR_1 is created.
    # The directory PYTHON_PACKAGE_DIR_2 always exists..
    #
    file(TO_NATIVE_PATH ${CMAKE_BINARY_DIR}/python/package/build/lib PYTHON_PACKAGE_DIR_1)
    file(TO_NATIVE_PATH ${CMAKE_BINARY_DIR}/python/package PYTHON_PACKAGE_DIR_2)

    if(WIN32)
        set(PROJECT_PATH "${PROJECT_LIBRARY_DIR}\;$ENV{PATH}")
        set(PROJECT_PYTHONPATH "${PYTHON_PACKAGE_DIR_1}\;${PYTHON_PACKAGE_DIR_2}\;$ENV{PYTHONPATH}")
    else()
        # Surround path with '' to avoid cmake errors when ( or ) exists.
        # Example: `Program Files (x86)` when compiling in WSL.
        set(PROJECT_PATH '$ENV{PATH}')
        set(PROJECT_PYTHONPATH '${PYTHON_PACKAGE_DIR_1}:${PYTHON_PACKAGE_DIR_2}:$ENV{PYTHONPATH}')
    endif()

    message("PYTHONPATH: ${PYTHONPATH}")
    message("PROJECT_PYTHONPATH: ${PROJECT_PYTHONPATH}")
    message("PYTHON_PACKAGE_DIR_1: ${PYTHON_PACKAGE_DIR_1}")
    message("PYTHON_PACKAGE_DIR_2: ${PYTHON_PACKAGE_DIR_2}")

    # Note that this must come after setting the Python path
    add_subdirectory(python)
endif()