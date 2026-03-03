#ifndef GEOCHEM_TESTCASE_IS_DEF
#define GEOCHEM_TESTCASE_IS_DEF

#include <memory>
#include <string>

/* Abstract base class for holding the name and contents (i.e., input "file")
*  of a geochemical simulation test case. */
struct GeoChemTestCase
{
    virtual ~GeoChemTestCase(){};
    virtual std::string name() const = 0;
    virtual std::string getInputFileAsString() const = 0;
};

using GeoChemTestCasePtr = std::unique_ptr<GeoChemTestCase>;

#endif
