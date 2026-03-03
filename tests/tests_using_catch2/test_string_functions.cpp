#include <catch2/catch.hpp>

#include <string>
#include <vector>

#include <opm/simulators/geochemistry/Utility/HelpfulStringMethods.hpp>

TEST_CASE("Test partitioning of vector of two elements")
{
    std::vector<std::string> strVector = {"Hello world", "Hello key world"};
    partition_vector_by_key(strVector , "KEY");

    REQUIRE(strVector[0].compare("Hello key world") == 0);
    REQUIRE(strVector[1].compare("Hello world") == 0);
}

TEST_CASE("Test partitioning of vector of three elements")
{
    std::vector<std::string> strVector = {"a", "B", "c"};
    partition_vector_by_key(strVector , "b");  // note: case-insensitive

    REQUIRE(strVector[0].compare("B") == 0);
}

TEST_CASE("Test trimming of strings")
{
    REQUIRE(trim_left("a string ") == "a string ");
    REQUIRE(trim_right("a string ") == "a string");

    REQUIRE(trim_left(" a string") == "a string");
    REQUIRE(trim_right(" a string") == " a string");

    REQUIRE(trim_left_right(" a string ") == "a string");
}

TEST_CASE("Test for strings with or without content")
{
    REQUIRE(has_content("      x      ") == true);
    REQUIRE(has_content("            ") == false);
}
