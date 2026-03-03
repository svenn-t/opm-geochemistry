// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright 2025 Equinor ASA.

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
#ifndef GEOCHEMISTRY_MODEL_HPP
#define GEOCHEMISTRY_MODEL_HPP

#include <dune/grid/common/gridenums.hh>
#include <dune/grid/common/partitionset.hh>
#include <dune/istl/bvector.hh>

#include <opm/common/OpmLog/OpmLog.hpp>

#include <opm/input/eclipse/EclipseState/Geochemistry/GenericSpeciesConfig.hpp>
#include <opm/input/eclipse/Schedule/Well/Well.hpp>
#include <opm/input/eclipse/Schedule/Well/WellTracerProperties.hpp>

#include <opm/grid/utility/ElementChunks.hpp>

#include <opm/models/parallel/threadmanager.hpp>
#include <opm/models/utils/propertysystem.hh>

#include <opm/simulators/geochemistry/OpmGeoChemInterface.hpp>
#include <opm/simulators/wells/WellTracerRate.hpp>

#include <fmt/format.h>

#include <algorithm>
#include <filesystem>
#include <memory>
#include <utility>
#include <vector>

namespace Opm {

class EclipseState;
class Well;

template <class TypeTag>
class GeochemistryModel
{
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using Grid = GetPropType<TypeTag, Properties::Grid>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;

    using CartesianIndexMapper = Dune::CartesianIndexMapper<Grid>;
    using SpeciesVector = Dune::BlockVector<Scalar>;

    enum { waterPhaseIdx = FluidSystem::waterPhaseIdx };

    static constexpr int dimWorld = Grid::dimensionworld;

public:
    /*!
    * \brief Constructor
    *
    * \param simulator Reference to simulator object
    */
    explicit GeochemistryModel(Simulator& simulator)
        : simulator_ (simulator)
        , eclState_  (simulator.vanguard().eclState())
        , cartMapper_(simulator.vanguard().cartesianIndexMapper())
        , element_chunks_(simulator.gridView(), Dune::Partitions::all, ThreadManager::maxThreads())
    {}

    /*!
    * \brief Initialize geochemistry model
    *
    * Initialize the geochemistry solver and concentration vectors
    */
    void init()
    {
        // Return if GEOCHEM is not activated
        const auto& geochem = eclState_.runspec().geochem();
        if (!geochem.enabled()) {
            return;
        }

        // Get transported species names
        const auto& species = eclState_.species();
        std::transform(species.begin(), species.end(), std::back_inserter(speciesNames_),
                       [] (const auto& item) { return item.name; } );

        // Get mineral names
        const auto& mineral = eclState_.mineral();
        if (!mineral.empty()) {
            mineralNames_.reserve(mineral.size());
            std::transform(mineral.begin(), mineral.end(), std::back_inserter(mineralNames_),
                           [] (const auto& item) { return item.name; } );
        }

        // Get ion exchange names
        const auto& ion_exchange = eclState_.ionExchange();
        if (!ion_exchange.empty()) {
            ionExNames_.reserve(ion_exchange.size());
            std::transform(ion_exchange.begin(), ion_exchange.end(), std::back_inserter(ionExNames_),
                           [] (const auto& item) { return item.name; } );
        }

        // Initialize interface to geochemistry solver
        const auto& file_name = geochem.geochem_file_name();
        std::pair<double, double> tol = std::make_pair<double, double>(geochem.mbal_tol(), geochem.ph_tol());
        bool charge_balance = geochem.charge_balance();
        int splay_tree_resolution = geochem.splay_tree_resolution();
        geoChemInterface_ = std::make_shared<OpmGeoChemInterface>();
        geoChemInterface_->initialize_from_opm_deck(file_name,
                                                    speciesNames_,
                                                    mineralNames_,
                                                    ionExNames_,
                                                    charge_balance,
                                                    tol,
                                                    splay_tree_resolution);

        // Initialize species independent vectors
        const std::size_t numGridDof = simulator_.model().numGridDof();
        pH_.resize(numGridDof, 7.0);
        sigma_.resize(numGridDof, 0.0);
        psi_.resize(numGridDof, 0.0);
        initial_equil_.resize(numGridDof, true);
        vol1_.resize(numGridDof);

        // Fill in species concentrations
        const std::size_t nSpecies = numSpecies();
        concentration_.resize(nSpecies);
        concentrationInitial_.resize(nSpecies);
        Cads_.resize(nSpecies);
        for (std::size_t speciesIdx = 0; speciesIdx < nSpecies;  ++speciesIdx) {
            const auto& single_species = species[speciesIdx];
            concentration_[speciesIdx].resize(numGridDof);
            concentrationInitial_[speciesIdx].resize(numGridDof);
            Cads_[speciesIdx].resize(numGridDof);

            // Initial concentration for species
            setInitialConcentrations_(single_species, concentration_[speciesIdx]);
        }

        // Initialize mineral concentrations
        const std::size_t nMin = geoChemInterface_->numberOfMinerals();
        if (nMin > 0) {
            Cmin_.resize(nMin);
            minWt_.resize(nMin);
            for (std::size_t minSpeciesIdx = 0; minSpeciesIdx < nMin; ++minSpeciesIdx) {
                const auto& single_mineral = mineral[minSpeciesIdx];
                Cmin_[minSpeciesIdx].resize(numGridDof);
                minWt_[minSpeciesIdx].resize(numGridDof);

                // Initial weight fraction for mineral
                setInitialConcentrations_(single_mineral, minWt_[minSpeciesIdx]);
            }
        }

        // Initialize ion exchange species
        const std::size_t nIo = geoChemInterface_->numberOfIonExchange();
        if (nIo > 0) {
            Cio_.resize(nIo);
            for (std::size_t ioIdx = 0; ioIdx < nIo; ++ioIdx) {
                const auto& single_io = ion_exchange[ioIdx];
                Cio_[ioIdx].resize(numGridDof);

                // Initial ion exchange
                setInitialConcentrations_(single_io, Cio_[ioIdx]);
            }
        }
    }

    /*!
    * \brief Calculations before a time integration
    *
    * Updates variables from the previous time step to be used in endTimeStep()
    */
    void beginTimeStep()
    {
        // Return if GEOCHEM is not activated
        if (!eclState_.runspec().geochem().enabled()) {
            return;
        }

        // Store variables from previous time step
        updateStorageCache();
    }

    /*!
    * \brief Calculations after a time integration
    *
    * \note The reactive transport step is done here!
    */
    void endTimeStep()
    {
        if (!eclState_.runspec().geochem().enabled()) {
            return;
        }

        // Reactive transport
        advanceSpeciesFieldsExplicit();
    }

    /*!
    * \brief Get number of species
    *
    * \returns Number of transported species
    */
    std::size_t numSpecies() const
    {
        return eclState_.species().size();
    }

    /*!
    * \brief Get number of minerals
    *
    * \returns Number of minerals
    */
    std::size_t numMinerals() const
    {
        return eclState_.mineral().size();
    }

    /*!
    * \brief Get number of ion exchange species
    *
    * \returns Number of ion exchange
    */
    std::size_t numIonEx() const
    {
        return eclState_.ionExchange().size();
    }

    /*!
    * \brief Get a particular species name
    *
    * \param speciesIdx Index of species
    * \returns String with queried species name
    */
    const std::string& speciesName(unsigned speciesIdx) const
    {
        return speciesNames_[speciesIdx];
    }

    /*!
    * \brief Get (reference to) a particular mineral name
    *
    * \param minSpeciesIdx Index of mineral
    * \returns String with queried mineral name
    */
    const std::string& mineralName(unsigned minSpeciesIdx) const
    {

        return mineralNames_[minSpeciesIdx];
    }

    /*!
    * \brief Get (reference to) a particular ion exchange name
    *
    * \param minSpeciesIdx Index of ion exchange
    * \returns String with queried ion exchage name
    */
    const std::string& ionExchangeName(unsigned ionExIdx) const
    {

        return ionExNames_[ionExIdx];
    }

    /*!
    * \brief Get WSPECIES for a particular species
    *
    * \param eclWell Reference to eclWell object
    * \param name Name of species
    * \param summaryState SummaryState object
    * \returns Concentration of injected species
    */
    Scalar currentWSPECIES_(const Well& eclWell, const std::string& name, const SummaryState& summaryState) const
    {
        return eclWell.getSpeciesProperties().getConcentration(WellTracerProperties::Well { eclWell.name() },
                                                               WellTracerProperties::Tracer { name },
                                                               summaryState);
    }

    /*!
    * \brief Get concentration of a species at a grid cell
    *
    * \param speciesIdx Index of species
    * \param globalDofIdx Cell index
    * \returns Concentration of species
    */
    Scalar speciesConcentration(int speciesIdx, int globalDofIdx) const
    {
        if (concentration_.empty()) {
            return 0.0;
        }

        return concentration_[speciesIdx][globalDofIdx];
    }

    /*!
    * \brief Get concentration of a mineral at a grid cell
    *
    * \param speciesIdx Index of mineral
    * \param globalDofIdx Cell index
    * \returns Concentration of mineral
    */
    Scalar mineralConcentration(int minSpeciesIdx, int globalDofIdx) const
    {
        if (Cmin_.empty()) {
            return 0.0;
        }

        return Cmin_[minSpeciesIdx][globalDofIdx];
    }

    /*!
    * \brief Get pH in a grid cell
    *
    * \param globalDofIdx Cell index
    * \returns pH
    */
    Scalar PH(int globalDofIdx) const
    {
        return pH_[globalDofIdx];
    }

    /*!
    * \brief Get all "standard" wells' species concentration rates
    *
    * \returns Container with well species concentration rates
    */
    const std::unordered_map<int, std::vector<WellTracerRate<Scalar>>>&
    getWellSpeciesRates() const
    {
        return wellSpeciesRate_;
    }

    /*!
    * \brief Get all multisegmented wells' species concentration rates
    *
    * \returns Container with well species concentration rates
    */
    const std::unordered_map<int, std::vector<MSWellTracerRate<Scalar>>>&
    getMswSpeciesRates() const
    {
        return mSwSpeciesRate_;
    }

    /*!
    * \brief Set species contration rate in a grid cell
    *
    * \param speciesIdx Index of mineral
    * \param globalDofIdx Cell index
    * \param value Species concentration for grid cell
    */
    void setSpeciesConcentration(int speciesIdx, int globalDofIdx, Scalar value)
    {
        concentration_[speciesIdx][globalDofIdx] = value;
    }

    /*!
    * \brief Serialize variables
    *
    * \param serializer Byte array conversion
    */
    template<class Serializer>
    void serializeOp(Serializer& serializer)
    {
        serializer(concentration_);
        serializer(pH_);
        serializer(wellSpeciesRate_);
        serializer(mSwSpeciesRate_);
    }


protected:
    /*!
     * \brief Set initial concentrations for a (single) species
     *
     * \param single_species Species config
     * \param concentration Initial concentration vector to be set
     */
    void setInitialConcentrations_(const GenericSpeciesConfig::SpeciesEntry& single_species,
                                   SpeciesVector& concentration)
    {
        // *BLK
        if (single_species.concentration.has_value()) {
            const auto& species_concentration = single_species.concentration.value();
            assert(species_concentration.size() == concentration.size());
            if (species_concentration.size() != static_cast<std::size_t>(cartMapper_.cartesianSize())) {
                throw std::runtime_error("Size of S/M/IBLK" + single_species.name + " is wrong!");
            }

            // Copy *BLK concentrations for each cell to concentration vector
            std::copy_n(species_concentration.begin(), species_concentration.size(), concentration.begin());
        }
        // *VDP
        else if (single_species.svdp.has_value()) {
            const auto& species_svdp = single_species.svdp.value();
            const auto& centroids = simulator_.vanguard().cellCentroids();

            // For each each grid cell, evaluate *VDP and assign to concentration vector
            std::for_each(
                concentration.begin(), concentration.end(),
                [idx = 0, &species_svdp, &centroids](auto& conc) mutable
                {
                    conc = species_svdp.evaluate("SPECIES_CONCENTRATION", centroids(idx)[2]);
                    ++idx;
                }
            );
        }
        // Zero initial condition
        else {
            OpmLog::warning(fmt::format("No S/M/IBLK or S/M/IVDP given for species {}. "
                                        "Initial values set to zero. ", single_species.name));
            std::fill(concentration.begin(), concentration.end(), 0.0);
        }
    }

    /*!
    * \brief Compute black oil equation volume term
    *
    * \param globalDofIdx Cell index
    * \param timeIdx Time step index
    * \return Black oil volume term
    */
    Scalar computeVolume_(const unsigned globalDofIdx,
                          const unsigned timeIdx) const
    {
        const auto& intQuants = simulator_.model().intensiveQuantities(globalDofIdx, timeIdx);
        const auto& fs = intQuants.fluidState();
        constexpr Scalar min_volume = 1e-10;

        return std::max(decay<Scalar>(fs.saturation(waterPhaseIdx)) *
                        decay<Scalar>(fs.invB(waterPhaseIdx)) *
                        decay<Scalar>(intQuants.porosity()),
                        min_volume);
    }

    /*!
    * \brief Compute black oil equation flux term
    *
    * \param elemCtx Reference to element context object
    * \param scvfIdx (Local) control volume index
    * \param timeIdx Time step index
    * \return Black oil flux term and boolean indicating upstream cell or not
    */
    std::pair<Scalar, bool> computeFlux_(const ElementContext& elemCtx,
                                         const unsigned scvfIdx,
                                         const unsigned timeIdx) const
    {
        const auto& stencil = elemCtx.stencil(timeIdx);
        const auto& scvf = stencil.interiorFace(scvfIdx);

        const auto& extQuants = elemCtx.extensiveQuantities(scvfIdx, timeIdx);
        const unsigned inIdx = extQuants.interiorIndex();

        unsigned upIdx = extQuants.upstreamIndex(waterPhaseIdx);
        const auto& intQuants = elemCtx.intensiveQuantities(upIdx, timeIdx);
        const auto& fs = intQuants.fluidState();
        Scalar v = decay<Scalar>(extQuants.volumeFlux(waterPhaseIdx))
                   * decay<Scalar>(fs.invB(waterPhaseIdx));

        const Scalar A = scvf.area();
        return std::pair{A * v, inIdx == upIdx};
    }

    /*!
    * \brief Update cache for reactive transport step
    *
    * \warning Everything that is need from Flow -and is not stored before endTimeStep()- must be saved here!
    */
    void updateStorageCache()
    {
        // Update concentration from previous time step
        concentrationInitial_ = concentration_;

        // Parallel loop over element chunks
        #ifdef _OPENMP
        #pragma omp parallel for
        #endif
        for (const auto& chunk : element_chunks_) {
            ElementContext elemCtx(simulator_);

            for (const auto& elem : chunk) {
                elemCtx.updatePrimaryStencil(elem);
                elemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);
                const Scalar extrusionFactor = elemCtx.intensiveQuantities(/*dofIdx=*/ 0, /*timeIdx=*/0).extrusionFactor();
                const Scalar scvVolume = elemCtx.stencil(/*timeIdx=*/0).subControlVolume(/*dofIdx=*/ 0).volume() * extrusionFactor;
                const unsigned globalDofIdx = elemCtx.globalSpaceIndex(0, /*timeIdx=*/0);

                vol1_[globalDofIdx] = computeVolume_(globalDofIdx, 0) * scvVolume;
            }
        }
    }

    /*!
    * \brief Run reactive transport solver and post-processing
    *
    * \note The actual reactive transport solver is in speciesEquationsExplicit_()
    */
    void advanceSpeciesFieldsExplicit()
    {
        // Calculate new concentration fields for each species
        speciesEquationsExplicit_();

        // Post-processing for concentration output
        constexpr Scalar tol_sat = 1e-6;
        for (std::size_t sIdx = 0; sIdx < numSpecies(); ++sIdx) {
            for (std::size_t globalDofIdx = 0; globalDofIdx < concentration_[sIdx].size(); ++globalDofIdx) {
                const auto& intQuants = simulator_.model().intensiveQuantities(globalDofIdx, 0);
                const auto& fs = intQuants.fluidState();
                const Scalar Sw = decay<Scalar>(fs.saturation(waterPhaseIdx));

                if (concentration_[sIdx][globalDofIdx] < 0.0 || Sw < tol_sat) {
                    concentration_[sIdx][globalDofIdx] = 0.0;
                }
            }
        }

        // Report produced species
        const auto& wellPtrs = simulator_.problem().wellModel().localNonshutWells();
        for (const auto& wellPtr : wellPtrs) {
            const auto& eclWell = wellPtr->wellEcl();

            // Injection rates already reported in speciesEquationWellExplicit_()
            if (!eclWell.isProducer()) {
                continue;
            }

            Scalar rateWellPos = 0.0;
            Scalar rateWellNeg = 0.0;
            const std::size_t well_index = simulator_.problem().wellModel().wellState().index(eclWell.name()).value();
            const auto& ws = simulator_.problem().wellModel().wellState().well(well_index);
            auto& speciesRate = wellSpeciesRate_[eclWell.seqIndex()];
            auto* mswSpeciesRate = eclWell.isMultiSegment() ? &mSwSpeciesRate_[eclWell.seqIndex()] : nullptr;

            for (std::size_t i = 0; i < ws.perf_data.size(); ++i) {
                const auto I = ws.perf_data.cell_index[i];
                const Scalar rate = wellPtr->volumetricSurfaceRateForConnection(I, waterPhaseIdx);

                if (rate < 0.0) {
                    for (std::size_t sIdx = 0; sIdx < numSpecies(); ++sIdx) {
                        const Scalar delta = rate * concentration_[sIdx][I];
                        speciesRate[sIdx].rate += delta;
                        if (eclWell.isMultiSegment()) {
                            (*mswSpeciesRate)[sIdx].rate[eclWell.getConnections().get(i).segment()]
                                += delta;
                        }
                    }
                }

                if (rate < 0) {
                    rateWellNeg += rate;
                } else {
                    rateWellPos += rate;
                }
            }

            // TODO: Some inconsistencies here that perhaps should be clarified.
            // The "offical" rate as reported below is occasionally significant
            // different from the sum over connections (as calculated above). Only observed
            // for small values, neglible for the rate itself, but matters when used to
            // calculate tracer concentrations.
            // const Scalar official_well_rate_total =
            //     simulator_.problem().wellModel().wellState().well(well_index).surface_rates[waterPhaseIdx];

            // const Scalar rateWellTotal = official_well_rate_total;

            // if (rateWellTotal > rateWellNeg) { // Cross flow
            //     constexpr Scalar bucketPrDay
            //         = 10.0 / (1000. * 3600. * 24.); // ... keeps (some) trouble away
            //     const Scalar factor = (rateWellTotal < -bucketPrDay) ? rateWellTotal / rateWellNeg : 0.0;
            //     for (std::size_t sIdx = 0; sIdx < numSpecies(); ++sIdx) {
            //         speciesRate[sIdx].rate *= factor;
            //     }
            // }
        }
    }

    /*!
    * \brief Reactive transport solver with explicit scheme
    */
    void speciesEquationsExplicit_()
    {
        // Clear well containers
        wellSpeciesRate_.clear();
        mSwSpeciesRate_.clear();

        // Reserve new space
        const auto& wellPtrs = simulator_.problem().wellModel().localNonshutWells();
        wellSpeciesRate_.reserve(wellPtrs.size());
        mSwSpeciesRate_.reserve(mSwSpeciesRate_.size());

        // Simulator information
        ElementContext elemCtx(simulator_);
        const Scalar dt = elemCtx.simulator().timeStepSize();

        // Initialize current species concentrations to zero
        for (std::size_t sIdx = 0; sIdx < numSpecies(); ++sIdx) {
            concentration_[sIdx] = 0.0;
        }

        // Calculate well terms
        for (const auto& wellPtr : wellPtrs) {
            speciesEquationWellExplicit_(*wellPtr, dt);
        }

        // Loop over grid blocks and calculate new concentrations for each species
        for (const auto& elem : elements(simulator_.gridView())) {
            elemCtx.updateStencil(elem);
            const std::size_t I = elemCtx.globalSpaceIndex(/*dofIdx=*/ 0, /*timeIdx=*/0);

            // Dirichlet BC
            if (elem.partitionType() != Dune::InteriorEntity) {
                for (std::size_t sIdx = 0; sIdx < numSpecies(); ++sIdx) {
                    concentration_[sIdx][I] = 0.0;
                }
                continue;
            }

            // Update block quantities
            elemCtx.updateAllIntensiveQuantities();
            elemCtx.updateAllExtensiveQuantities();

            // Volume at current time step
            const Scalar extrusionFactor =
                elemCtx.intensiveQuantities(/*dofIdx=*/ 0, /*timeIdx=*/0).extrusionFactor();
            Valgrind::CheckDefined(extrusionFactor);
            assert(isfinite(extrusionFactor));
            assert(extrusionFactor > 0.0);
            const Scalar scvVolume =
                elemCtx.stencil(/*timeIdx=*/0).subControlVolume(/*dofIdx=*/ 0).volume() * extrusionFactor;
            const Scalar vol = computeVolume_(I, 0);
            const Scalar vol0 = vol * scvVolume;

            // At simulation time == 0.0, equilibrate geochemical system
            if (initial_equil_[I]) {
                speciesEquationChemistryExplicit_(dt, I, initial_equil_[I]);
                initial_equil_[I] = false;
            }

            // Calculate volume/storage term
            speciesEquationVolumeExplicit_(I);

            // Calculate flux term
            const std::size_t numInteriorFaces = elemCtx.numInteriorFaces(/*timIdx=*/0);
            for (unsigned scvfIdx = 0; scvfIdx < numInteriorFaces; scvfIdx++) {
                // Get neighbour global index
                const auto& face = elemCtx.stencil(0).interiorFace(scvfIdx);
                const unsigned j = face.exteriorIndex();
                const unsigned J = elemCtx.globalSpaceIndex(/*dofIdx=*/ j, /*timIdx=*/0);

                // Flux term between I and J
                speciesEquationFluxExplicit_(elemCtx, scvfIdx, I, J, dt);
            }

            // Divide all terms with current volume (necessary in explicit scheme)
            for (std::size_t sIdx = 0; sIdx < numSpecies(); ++sIdx) {
                concentration_[sIdx][I] /= vol0;
            }

            // Equilibrate geochemical system
            speciesEquationChemistryExplicit_(dt, I);
        }
    }

    /*!
    * \brief Calculate single well contribution in explicit reactive transport solver
    *
    * \param well Reference to well object
    * \param dt Time step
    */
    template <class Well>
    void speciesEquationWellExplicit_(const Well& well, const Scalar& dt)
    {
        // Get simulation wells
        const auto& eclWell = well.wellEcl();

        // Reserve space for species output
        auto& speciesRate = wellSpeciesRate_[eclWell.seqIndex()];
        speciesRate.reserve(numSpecies());
        auto* mswSpeciesRate = eclWell.isMultiSegment()
            ? &mSwSpeciesRate_[eclWell.seqIndex()]
            : nullptr;
        if (mswSpeciesRate) {
            mswSpeciesRate->reserve(numSpecies());
        }

        // Init. well output to zero
        for (std::size_t sIdx = 0; sIdx < numSpecies(); ++sIdx) {
            speciesRate.emplace_back(speciesName(sIdx), 0.0);
            if (eclWell.isMultiSegment()) {
                auto& wsr = mswSpeciesRate->emplace_back(speciesName(sIdx));
                wsr.rate.reserve(eclWell.getConnections().size());
                for (std::size_t i = 0; i < eclWell.getConnections().size(); ++i) {
                    wsr.rate.emplace(eclWell.getConnections().get(i).segment(), 0.0);
                }
            }
        }

        // Get WSPECIES info
        std::vector<Scalar> wspecies(numSpecies());
        for (std::size_t sIdx = 0; sIdx < numSpecies(); ++sIdx) {
            wspecies[sIdx] = currentWSPECIES_(eclWell, speciesName(sIdx),
                                              simulator_.problem().wellModel().summaryState());
        }

        // Calculate well term
        const auto& ws = simulator_.problem().wellModel().wellState().well(well.name());
        for (std::size_t i = 0; i < ws.perf_data.size(); ++i) {
            // Get perforation rate
            const auto I = ws.perf_data.cell_index[i];
            const Scalar rate = well.volumetricSurfaceRateForConnection(I, waterPhaseIdx);

            // Injection
            if (rate > 0.0) {
                for (std::size_t sIdx = 0; sIdx < numSpecies(); ++sIdx) {
                    // Inject WSPECIES concentration
                    const Scalar inj_species_rate = rate * wspecies[sIdx];
                    concentration_[sIdx][I] += dt * inj_species_rate;

                    // Store for reporting here because WSPECIES is constant over time step
                    speciesRate[sIdx].rate += inj_species_rate;
                    if (eclWell.isMultiSegment()) {
                        (*mswSpeciesRate)[sIdx].rate[eclWell.getConnections().get(i).segment()] += inj_species_rate;
                    }
                }
            }
            // Production
            // OBS: storing well rates for reporting done in advanceSpeciesFieldsExplicit_()
            else if (rate < 0.0) {
                for (std::size_t sIdx = 0; sIdx < numSpecies(); ++sIdx) {
                    // Produce tracer concentration
                    concentration_[sIdx][I] += dt * rate * concentrationInitial_[sIdx][I];

                    // Ensure reporting of cross-flow
                    const Scalar inj_species_rate = rate * wspecies[sIdx];
                    speciesRate[sIdx].rate += inj_species_rate;
                }
            }
        }
    }

    /*!
    * \brief Volume term in explicit reactive transport solver
    *
    * \param I Grid cell index
    */
    void speciesEquationVolumeExplicit_(unsigned I)
    {
        for (std::size_t sIdx = 0; sIdx < numSpecies(); ++sIdx) {
            concentration_[sIdx][I] += vol1_[I] * concentrationInitial_[sIdx][I];
        }
    }

    /*!
    * \brief Flux term in explicit reactive transport solver
    *
    * \param elemCtx Reference to element context object
    * \param scvfIdx (Local) control volume index
    * \param I Index of focused cell
    * \param J Index of neighbor cell to I
    * \param dt Time step
    */
    void speciesEquationFluxExplicit_(const ElementContext& elemCtx,
                                      unsigned scvfIdx,
                                      unsigned I,
                                      unsigned J,
                                      const Scalar dt)
    {
        const auto& [flux, isUp] = computeFlux_(elemCtx, scvfIdx, 0);
        const int globalUpIdx = isUp ? I : J;
        for (std::size_t sIdx = 0; sIdx < numSpecies(); ++sIdx) {
            concentration_[sIdx][I] -= dt * flux * concentrationInitial_[sIdx][globalUpIdx];
        }
    }

    /*!
    * \brief Run geochemistry equilibrium solver
    *
    * \param dt Time step
    * \param I Cell index
    * \param initial_equil Optional bool indicating if this is an initial equilibrium solve
    */
    void speciesEquationChemistryExplicit_(Scalar dt,
                                           unsigned I,
                                           bool initial_equil = false)
    {
        // Ensure that dt == 0 for initial equilibration
        const double dt_geochem = initial_equil ? 0.0 : dt;

        // Initialize total and adsorbed concentrations before equilibration
        std::vector<double> Ctot;
        std::vector<double> Cads;
        const auto nAqu = geoChemInterface_->numberOfAqueousBasisSpecies();
        const auto nBasis = geoChemInterface_->numberOfBasisSpecies();
        Ctot.resize(nBasis, 0.0);
        Cads.resize(nBasis, 0.0);
        for (std::size_t k = 0; k < nAqu; ++k) {
            if (!initial_equil) {
                Cads[k] = Cads_[k][I];
            }
            Ctot[k] = initial_equil ?
                concentrationInitial_[k][I] : concentration_[k][I];
            Ctot[k] += Cads[k];
        }

        // Setup mineral concentrations before equilibration
        // NOTE: assert checks that number of internal minerals are same as in the OPM deck
        const std::size_t nMin = geoChemInterface_->numberOfMinerals();
        assert(nMin == numMinerals());

        double* log_Amin_ptr = nullptr;
        double* Cmin_ptr = nullptr;
        std::vector<double> Cmin;
        const auto& intQuants = simulator_.model().intensiveQuantities(I, 0);
        const auto& fs = intQuants.fluidState();
        const auto poro = decay<double>(intQuants.porosity());
        if (nMin > 0) {
            // Initial equilibrium solve requires calculation of initial mineral concentration
            if (initial_equil) {
                std::unordered_map<std::string, double> wt_frac;
                wt_frac.reserve(nMin);
                for (std::size_t l = 0; l < nMin; ++l) {
                    wt_frac.emplace(mineralNames_[l], minWt_[l][I]);
                }
                geoChemInterface_->calculate_initial_mineral_concentration(Cmin, poro, wt_frac);
                for (std::size_t l = 0; l < nMin; ++l) {
                    Cmin_[l][I] = Cmin[l];
                }
            }
            else {
                Cmin.resize(nMin);
                for (std::size_t l = 0; l < nMin; ++l) {
                    Cmin[l] = Cmin_[l][I];
                }
            }
            Cmin_ptr = Cmin.data();
            log_Amin_ptr = geoChemInterface_->get_log_a_mineral().data();
        }

        // Fluid properties
        const auto temp = decay<double>(fs.temperature(0));
        const auto pres = decay<double>(fs.pressure(waterPhaseIdx));

        // Phase saturations
        const double swat = decay<double>(fs.saturation(FluidSystem::waterPhaseIdx));
        double soil = 0.0;
        double sgas = 0.0;
        if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)) {
            soil = decay<double>(fs.saturation(FluidSystem::oilPhaseIdx));
        }
        if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
            sgas = decay<double>(fs.saturation(FluidSystem::gasPhaseIdx));
        }
        std::array<double, 3> mass_phase = { swat, soil, sgas };

        // Surface area
        double SA = geoChemInterface_->GetSurfaceArea();

        // Diffusion layer
        // OBS: calculated in set_surface_concentrations()!
        double frac_DL = 0.0;

        // Set surface concentration
        // NOTE: assert checks that number of internal ion exchange species are same as in the OPM deck
        const auto nIon = geoChemInterface_->numberOfIonExchange();
        assert(nIon == numIonEx());

        // Set concentrations for ion exchangers, surface complexes, and diffusive layers in Ctot
        std::unordered_map<std::string, double> Cion;
        if (nIon > 0) {
            Cion.reserve(nIon);
            for (std::size_t m = 0; m < nIon; ++m) {
                Cion.emplace(ionExNames_[m], Cio_[m][I]);
            }
        }
        geoChemInterface_->set_surface_concentrations(mass_phase[0], Ctot, frac_DL, Cion);

        // misc variables
        double pH = pH_[I];
        double sigma = sigma_[I];
        double psi = psi_[I];

        // Run equilibrium solver
        geoChemInterface_->SolveChem_I(Ctot.data(),
                                       Cads.data(),
                                       Cmin_ptr,
                                       log_Amin_ptr,
                                       temp,
                                       pres,
                                       poro,
                                       dt_geochem,
                                       SA,
                                       frac_DL,
                                       mass_phase,
                                       pH,
                                       sigma,
                                       psi);

        // Update mineral concentarations
        if (nMin > 0) {
            for (std::size_t l = 0; l < nMin; ++l) {
                Cmin_[l][I] += Cmin[l];
            }
        }

        // Update misc variables
        pH_[I] = pH;
        sigma_[I] = sigma;
        psi_[I] = psi;

        // Update species concentrations
        for (std::size_t k = 0; k < nAqu; ++k) {
            Ctot[k] -= Cads[k];
            Cads_[k][I] = Cads[k];

            if (initial_equil) {
                concentrationInitial_[k][I] = Ctot[k];
            }
            else {
                concentration_[k][I] = Ctot[k];
            }
        }
    }

private:
    std::vector<SpeciesVector> concentrationInitial_;
    std::vector<SpeciesVector> concentration_;
    std::vector<SpeciesVector> Cads_;
    std::vector<SpeciesVector> Cio_;
    std::vector<SpeciesVector> Cmin_;
    std::vector<SpeciesVector> minWt_;
    std::vector<double> pH_;
    std::vector<double> sigma_;
    std::vector<double> psi_;
    std::vector<bool> initial_equil_;
    std::vector<Scalar> vol1_;
    std::vector<std::string> speciesNames_;
    std::vector<std::string> mineralNames_;
    std::vector<std::string> ionExNames_;

    std::shared_ptr<OpmGeoChemInterface> geoChemInterface_;

    std::unordered_map<int, std::vector<WellTracerRate<Scalar>>> wellSpeciesRate_;
    std::unordered_map<int, std::vector<MSWellTracerRate<Scalar>>> mSwSpeciesRate_;

    Simulator& simulator_;
    const EclipseState& eclState_;
    const CartesianIndexMapper& cartMapper_;
    ElementChunks<GridView, Dune::Partitions::All> element_chunks_;
};  // class GeochemistryModel
} // namespace Opm

#endif