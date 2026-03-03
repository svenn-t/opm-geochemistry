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
#ifndef NETCDF_IS_DEF_H
#define NETCDF_IS_DEF_H

#include <iostream>
#include <string>

#ifdef NETCDF
#include <netcdf>
#endif

/*
 * Class for writing data to the NetCDF format.
 * Currently, it is only used for 1D time-dependent grid data.
 */
class NetCDFWriter{

public:

    NetCDFWriter(const std::string& file_name);

    // The .nc file is automatically closed by the destructor.
    // This frees up any internal netCDF resources associated with the file,
    // and flushes any buffers.
    ~NetCDFWriter() = default;
    
    void incrementTime(double time, double pvi=-1);
    void initialize(int nx, double* x_coordinates);
    int writeField(const std::string& field_name, double* field_values);

private:

    int nx_{0};
    int currentTimeIndex_{-1};
    double currentTime_{0.0};
    double currentPoreVolumes_{0.0};
    
#ifdef NETCDF
    netCDF::NcFile netCDF_file_;
    std::vector<netCDF::NcDim> dimensions_;
#endif

};

#endif
