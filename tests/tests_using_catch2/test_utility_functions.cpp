#if __has_include(<catch2/catch.hpp>)
#include <catch2/catch.hpp>
#else
#include <catch2/catch_test_macros.hpp>
#endif

#include <opm/simulators/geochemistry/Utility/MiscUtilityFunctions.hpp>

#include <limits>
#include <vector>

TEST_CASE("Test fill_beginning_of_vector()")
{
    std::vector<double> src = { 1.0, 2.0, 3.0};
    std::vector<double> dest;

    fill_beginning_of_vector_from_vector(src, dest);
    REQUIRE(dest.size() == 0); // function does nothing in this case

    dest.resize(1);
    fill_beginning_of_vector_from_vector(src, dest);
    REQUIRE(dest.size() == 1);
    REQUIRE(dest[0] == src[0]);

    dest.resize(2);
    fill_beginning_of_vector_from_vector(src, dest);
    for(std::size_t i=0; i < dest.size(); ++i) REQUIRE(dest[i] == src[i]);

    dest.resize(3);
    fill_beginning_of_vector_from_vector(src, dest);
    for(std::size_t i=0; i < dest.size(); ++i) REQUIRE(dest[i] == src[i]);

    dest.resize(4);
    dest[3] = 42.0;
    fill_beginning_of_vector_from_vector(src, dest);
    for(std::size_t i=0; i < 3; ++i) REQUIRE(dest[i] == src[i]);
    REQUIRE(dest[3] == 42.0); // Only 3 elements in src vector

    std::fill(dest.begin(), dest.end(), 0.0);
    fill_beginning_of_vector_from_vector(src, dest, 1);
    REQUIRE(dest[0] == src[0]);
    for(std::size_t i=1; i < dest.size(); ++i) REQUIRE(dest[i] == 0.0);

}

TEST_CASE("Test vector_has_infinite_numbers()")
{
    auto vector_with_inf = std::vector<double>{ 1.0, 2.0, std::numeric_limits<double>::infinity() };
    auto vector_without_inf = std::vector<double>{ 1.0, 2.0, 3.0 };

    REQUIRE(vector_has_infinite_numbers(vector_with_inf));
    REQUIRE(!vector_has_infinite_numbers(vector_without_inf));
}
