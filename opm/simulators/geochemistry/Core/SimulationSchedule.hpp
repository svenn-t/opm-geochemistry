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
#ifndef SIMULATION_SCHEDULE_IS_DEF_H
#define SIMULATION_SCHEDULE_IS_DEF_H

#include <vector>

struct GeochemBoundaryCondition {

    GeochemBoundaryCondition()
        : solution_index_(0)
        , flow_rate_(0.0)
    {
    }

    GeochemBoundaryCondition(int solution_index, double flow_rate)
        : solution_index_(solution_index)
        , flow_rate_(flow_rate)
    {
    }

    int solution_index_;
    double flow_rate_;

};

struct SimulationEpoch {

    SimulationEpoch()
        : start_time_(0.0)
        , end_time_(0.0)
        , inlet_bc_()
    {
    }

    SimulationEpoch(double start_time, double end_time, GeochemBoundaryCondition inlet_bc)
        : start_time_(start_time)
        , end_time_(end_time)
        , inlet_bc_(inlet_bc)
    {
    }

    double start_time_;
    double end_time_;
    GeochemBoundaryCondition inlet_bc_;

};

class SimulationSchedule {

public:

    SimulationSchedule();

    /* Checks that the schedule is in a valid state, meaning:
    * 
    *      1) The first epoch should start at t=0.
    *      2) The end time of one epoch must equal the start time of the next.
    */
    bool isValid() const;

    void resetSchedule();

    void addEpoch(SimulationEpoch epoch);
    void addEpoch(double start_time, double end_time, GeochemBoundaryCondition inlet_bc);
    void addEpoch(double dt, GeochemBoundaryCondition inlet_bc);

    const SimulationEpoch& getCurrentEpoch() const;
    double getTimeOfCurrentEpoch() const;
    void advanceToNextEpoch();

    bool simulationHasFinished() const;

private:

    std::size_t currentEpochIndex_;
    std::vector<SimulationEpoch> all_simulation_epochs_;

};

#endif
