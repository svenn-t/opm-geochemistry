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
#include <opm/simulators/geochemistry/Core/SimulationSchedule.hpp>

#include <assert.h>
#include <iostream>

SimulationSchedule::SimulationSchedule()
    : currentEpochIndex_(0)
    , all_simulation_epochs_()
{
}

bool SimulationSchedule::isValid() const
{
    if(all_simulation_epochs_.size() == 0) return false;

    double previous_end_time = 0.0;
    for (const auto& epoch : all_simulation_epochs_) {
        if (epoch.start_time_ != previous_end_time) return false;
        else previous_end_time = epoch.end_time_;
    }
    return true;
}

void SimulationSchedule::resetSchedule() {
    currentEpochIndex_ = 0;
}

void SimulationSchedule::addEpoch(SimulationEpoch epoch) {
    all_simulation_epochs_.push_back(epoch);
    assert(isValid());
};

void SimulationSchedule::addEpoch(double start_time, double end_time, GeochemBoundaryCondition inlet_bc) {

    const SimulationEpoch epoch(start_time, end_time, inlet_bc);
    addEpoch(epoch);

}

void SimulationSchedule::addEpoch(double dt, GeochemBoundaryCondition inlet_bc) {
    
    if (all_simulation_epochs_.size() > 0) {
        const double start_time = all_simulation_epochs_[all_simulation_epochs_.size() - 1].end_time_;
        const double end_time = start_time + dt;
        addEpoch(start_time, end_time, inlet_bc);

    }
    else {
        addEpoch(0.0, dt, inlet_bc);
    }
}

const SimulationEpoch& SimulationSchedule::getCurrentEpoch() const
{
    // Whenever we call this function, it is assumed that we have created
    // at least one epoch...
    assert(isValid());
    assert(!simulationHasFinished());
    return all_simulation_epochs_[currentEpochIndex_];
}

double SimulationSchedule::getTimeOfCurrentEpoch() const {

    if (all_simulation_epochs_.size() > 0) {
        const auto& epoch = getCurrentEpoch();
        const double dt = epoch.end_time_ - epoch.start_time_;
        assert(dt >= 0);
        return dt;
    }
    return 0.0;
}

void SimulationSchedule::advanceToNextEpoch() {

    if (!simulationHasFinished())
    {
        ++currentEpochIndex_;
    }

}

bool SimulationSchedule::simulationHasFinished() const {
    return currentEpochIndex_ == all_simulation_epochs_.size();
}
