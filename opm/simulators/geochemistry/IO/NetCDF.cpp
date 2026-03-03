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
#include <opm/simulators/geochemistry/IO/NetCDF.hpp>

/*
* Creates the file .nc-file, provided the netCDF library is present on the system and is detected.
* The Replace parameter tells netCDF to overwrite the file if it already exists.
*/
NetCDFWriter::NetCDFWriter(const std::string& file_name)
#ifdef NETCDF
: netCDF_file_(file_name, netCDF::NcFile::replace)
, dimensions_()
#endif
{
}

// Unsolved: How to handle multiple (unlimited) time-variables that should be synchronized at all times?
// (i.e., time in hours + pore volumes)
// For now, we only report regular times....
void NetCDFWriter::initialize(int nx, double* x_coordinates)
{
    nx_ = nx;

#ifdef NETCDF
    using namespace netCDF;

    // Define the dimensions
    NcDim tDim = netCDF_file_.addDim("T", NC_UNLIMITED);
    //NcDim pvDim = netCDF_file_.addDim("PV", NC_UNLIMITED);
    NcDim xDim = netCDF_file_.addDim("X", nx_);

    // Define coordinate netCDF variables
    NcVar tVar = netCDF_file_.addVar("T", ncFloat, tDim);
    //NcVar pvVar = netCDF_file_.addVar("PV", ncFloat, pvDim);
    NcVar xVar = netCDF_file_.addVar("X", ncFloat, xDim);

    float x_coord[nx_];
    for(int ix = 0; ix < nx_; ++ix){
        x_coord[ix] = x_coordinates[ix];
    }

    // Write the coordinate variable data
    xVar.putVar(x_coord);

    // Define units attributes for coordinate vars. This attaches a
    // text attribute to each of the coordinate variables, containing
    // the units. Note that we are not writing a trailing NULL, just
    // "units", because the reading program may be fortran which does
    // not use null-terminated strings. In general it is up to the
    // reading C program to ensure that it puts null-terminators on
    // strings where necessary.
    
    //tVar.putAtt("units", "hours"); // <-- xarray interprets this as a datetime object with time in nanoseconds??
    xVar.putAtt("units", "cell no."); // <-- store coordinates in cm instead?

    // Define the netCDF data variables
    dimensions_.push_back(tDim);
    //dimensions_.push_back(pvDim);
    dimensions_.push_back(xDim);

#endif
}

void NetCDFWriter::incrementTime(double time, double pvi)
{
    ++currentTimeIndex_;
    currentTime_ = time;
    currentPoreVolumes_ = pvi;
}

int NetCDFWriter::writeField(const std::string& field_name, double* field_values)
{

#ifdef NETCDF
    using namespace netCDF;
    using namespace netCDF::exceptions;

    // If variable has not already been defined in file, create it now
    if(netCDF_file_.getVar(field_name).isNull()){
        NcVar newVar = netCDF_file_.addVar(field_name, ncFloat, dimensions_);
    }

    std::vector<std::size_t> startp;
    startp.push_back(currentTimeIndex_);
    startp.push_back(0);

    std::vector<std::size_t> countp;
    countp.push_back(1);
    countp.push_back(nx_);

    try
    {
        auto timeVar = netCDF_file_.getVar("T");
        //auto pvVar = netCDF_file_.getVar("PV");
        auto currVar = netCDF_file_.getVar(field_name);

        // Question: Does it have to be a vector even if we add a single number??
        std::vector<std::size_t> t;
        t.push_back(currentTimeIndex_);
        timeVar.putVar(t, currentTime_);
        //pvVar.putVar(t, currentPoreVolumes_);

        currVar.putVar(startp, countp, &field_values[0]);
    }
    catch(NcException& e)
    {
        //e.what();
        // Return this to OS if there is a failure.
        #define NC_ERR 2
        return NC_ERR;

    }
#endif

    return 0;
}
