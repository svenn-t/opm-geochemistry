#ifndef GEOCHEM_TESTCASE_IS_DEF
#define GEOCHEM_TESTCASE_IS_DEF

#include <memory>
#include <string>

enum class GeoChemTestType{ EqSolver, OneDTransport };

/*
* Abstract base class for holding the name and contents (i.e., input "file")
* of a geochemical simulation test case.
*/
struct GeoChemTestCase
{
    virtual ~GeoChemTestCase(){};
    virtual std::string name() const = 0;
    virtual GeoChemTestType type() const = 0;
    virtual std::string getInputFileAsString() const = 0;
};

using GeoChemTestCasePtr = std::unique_ptr<GeoChemTestCase>;

#endif
