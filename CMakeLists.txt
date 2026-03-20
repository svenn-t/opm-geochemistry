###########################################################################
#                                                                         #
# Note: The bulk of the build system is located in the cmake/ directory.  #
#       This file only contains the specializations for this particular   #
#       project. Most likely you are interested in editing one of these   #
#       files instead:                                                    #
#                                                                         #
#       dune.module                              Name and version number  #
#       CMakeLists_files.cmake                   Path of source files     #
#       cmake/Modules/${project}-prereqs.cmake   Dependencies             #
#                                                                         #
###########################################################################

cmake_minimum_required (VERSION 3.23)
project(opm-geochemistry C CXX)

# Other options
option(SIBLING_SEARCH "Search for other modules in sibling directories?" ON)
option(BUILD_PYTHON "Attempt to build Python package?" OFF)

macro(opm-geochemistry_prereqs_hook)
  target_link_libraries(opmgeochemistry
    PUBLIC
      opmsimulators
      opmcommon
      opmgrid
  )

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

  if (BUILD_PYTHON)
    # Be backwards compatible.
    if(PYTHON_EXECUTABLE AND NOT Python3_EXECUTABLE)
      set(Python3_EXECUTABLE ${PYTHON_EXECUTABLE})
    endif()
    find_package(Python3 REQUIRED COMPONENTS Interpreter Development.Module)
    if(Python3_VERSION_MINOR LESS 3)
      # Python native namespace packages requires python >= 3.3
      message(SEND_ERROR "OPM requires python >= 3.3 but only version ${Python3_VERSION} was found")
    endif()

    find_package(pybind11 CONFIG)
    if(NOT pybind11_FOUND)
      include(DownloadPyBind11)
    endif()
  endif()
endmacro()

macro(opm-geochemistry_tests_hook)
  if(BUILD_TESTING)
    add_subdirectory(tests)
  endif()
endmacro()

macro(opm-geochemistry_targets_hook)
  # Generate flow_blackoil_geochemistry
  opm_add_library(
    TARGET
      moduleVersionGeochem
    TYPE
      OBJECT
    SOURCES
      opm/simulators/utils/moduleVersion.cpp
  )
  add_dependencies(moduleVersionGeochem opmgeochemistry)

  opm_add_test(flow_blackoil_geochemistry
               ONLY_COMPILE
               ALWAYS_ENABLE
               DEPENDS
                 opmsimulators
                 opmgeochemistry
               LIBRARIES
                 opmsimulators
                 opmgeochemistry
               SOURCES
                 flow/flow_blackoil_geochemistry.cpp
                 $<TARGET_OBJECTS:moduleVersionGeochem>
  )

  opm_add_test(flow_onephase_geochemistry
               ONLY_COMPILE
               ALWAYS_ENABLE
               DEPENDS
                 opmsimulators
                 opmgeochemistry
               LIBRARIES
                 opmsimulators
                 opmgeochemistry
               SOURCES
                 flow/flow_onephase_geochemistry.cpp
                 $<TARGET_OBJECTS:moduleVersionGeochem>
  )

  # Manually add OPM tests
  opm_add_test(test_opmgeocheminterface
               DEPENDS
                 opmsimulators
                 opmgeochemistry
               LIBRARIES
                 opmsimulators
                 opmgeochemistry
                 Boost::unit_test_framework
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

  if(BUILD_PYTHON)
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
    add_subdirectory(python)
  endif()
endmacro()

if(SIBLING_SEARCH AND NOT opm-common_DIR)
  # guess the sibling dir
  get_filename_component(_leaf_dir_name ${PROJECT_BINARY_DIR} NAME)
  get_filename_component(_parent_full_dir ${PROJECT_BINARY_DIR} DIRECTORY)
  get_filename_component(_parent_dir_name ${_parent_full_dir} NAME)
  #Try if <module-name>/<build-dir> is used
  get_filename_component(_modules_dir ${_parent_full_dir} DIRECTORY)
  if(IS_DIRECTORY ${_modules_dir}/opm-common/${_leaf_dir_name})
    set(opm-common_DIR ${_modules_dir}/opm-common/${_leaf_dir_name})
    set(opm-grid_DIR ${_modules_dir}/opm-grid/${_leaf_dir_name})
    set(opm-simulators_DIR ${_modules_dir}/opm-simulators/${_leaf_dir_name})
  else()
    string(REPLACE ${PROJECT_NAME} opm-common _opm_common_leaf ${_leaf_dir_name})
    if(NOT _leaf_dir_name STREQUAL _opm_common_leaf
        AND IS_DIRECTORY ${_parent_full_dir}/${_opm_common_leaf})
      # We are using build directories named <prefix><module-name><postfix>
      set(opm-common_DIR ${_parent_full_dir}/${_opm_common_leaf})
      string(REPLACE ${PROJECT_NAME} opm-grid _opm_grid_leaf ${_leaf_dir_name})
      set(opm-grid_DIR ${_parent_full_dir}/${_opm_grid_leaf})
      string(REPLACE ${PROJECT_NAME} opm-simulators _opm_simulators_leaf ${_leaf_dir_name})
      set(opm-simulators_DIR ${_parent_full_dir}/${_opm_grid_leaf})
    elseif(IS_DIRECTORY ${_parent_full_dir}/opm-common)
      # All modules are in a common build dir
      set(opm-common_DIR "${_parent_full_dir}/opm-common")
      set(opm-grid_DIR "${_parent_full_dir}/opm-grid")
      set(opm-simulators_DIR "${_parent_full_dir}/opm-grid")
    endif()
  endif()
endif()
if(opm-common_DIR AND NOT IS_DIRECTORY ${opm-common_DIR})
  message(WARNING "Value ${opm-common_DIR} passed to variable"
    " opm-common_DIR is not a directory")
endif()

find_package(opm-common REQUIRED)

# project information is in dune.module. Read this file and set variables.
# we cannot generate dune.module since it is read by dunecontrol before
# the build starts, so it makes sense to keep the data there then.
include (OpmInit)

# all setup common to the OPM library modules is done here
include (OpmLibMain)
