#if __has_include(<catch2/catch.hpp>)
#include <catch2/catch.hpp>
#else
#include <catch2/catch_test_macros.hpp>
#endif

#include <map>
#include <string>
#include <fstream>
#include <sstream>

#include <opm/simulators/geochemistry/IO/InputReader.hpp>

TEST_CASE("Test the reading of simple key-value pairs (InputReader)")
{
    // Define keywords with default values
    std::map<std::string, std::string> simple_key_value_pairs;
    simple_key_value_pairs["add_species"] = R"("ions.txt")";
    simple_key_value_pairs["interpolate"] = "0";
    simple_key_value_pairs["equilibrate"] = "0";
    InputReader IP_(simple_key_value_pairs, false);

    REQUIRE(IP_.get_simple_keyword_value("interpolate") == "0");
    REQUIRE(IP_.get_simple_keyword_value("equilibrate") == "0");
    
}


TEST_CASE("Test the reading of a single block keyword ") {
    
    const std::string input = R"(
this part
will be ignored


START

hi
this
# a comment line, ignore...
is
a
test
* aaand another comment

STOP

hello world, ignore this
too!
    )";

    std::istringstream inputStream(input);
    
    InputReader inpReader({}, /* case_sensitive= */ false);
    inpReader.define_block_keyword("START", { "STOP" });
    inpReader.read(inputStream);
    const auto& blockContents = inpReader.get_block_keyword_content("START");
    
    REQUIRE(blockContents.size() == 5);
    REQUIRE(blockContents[0].compare("hi") == 0);
    REQUIRE(blockContents[1].compare("this") == 0);
    REQUIRE(blockContents[2].compare("is") == 0);
    REQUIRE(blockContents[3].compare("a") == 0);
    REQUIRE(blockContents[4].compare("test") == 0);
    
}

TEST_CASE("Test the reading of two block keywords, where each of them has an associated name") {
    
    const std::string input = R"(
this part
will be ignored

SOLUTION NaCl
Na 0.1
Cl 0.1
/ end

hi
this
# a comment line, ignore...
is
a

BRINE Brine 2
Ca 0.1
Cl 0.2
/ end


test
* aaand another comment

ignore this
too!
    )";

    std::istringstream inputStream(input);
    
    InputReader inpReader({}, /* case_sensitive= */ false);
    inpReader.define_block_keyword("SOLUTION", { "/ end" });
    inpReader.define_block_keyword("BRINE", { "/ end" });
    inpReader.read(inputStream);
    
    // NB: The input block keyword itself is always used as a prefix,
    //     but we also need the next word on the same line... Hence these are empty:
    REQUIRE(inpReader.get_block_keyword_content("SOLUTION").size() == 0);
    REQUIRE(inpReader.get_block_keyword_content("BRINE").size() == 0);
    
    const auto& solutionContents = inpReader.get_block_keyword_content("SOLUTION NACL");
    REQUIRE(solutionContents.size() == 2);
    REQUIRE(solutionContents[0].compare("Na 0.1") == 0);
    REQUIRE(solutionContents[1].compare("Cl 0.1") == 0);
    
    // UH OH: The "2" is not accounted for, b/c we just read the first word after discarding whitespaces:
    const auto& brineContents = inpReader.get_block_keyword_content("BRINE BRINE");
    REQUIRE(brineContents.size() == 2);
    REQUIRE(brineContents[0].compare("Ca 0.1") == 0);
    REQUIRE(brineContents[1].compare("Cl 0.2") == 0);
    
}

TEST_CASE("Test the reading of two block keywords again, where each of them has an associated name, but where the case also matters!") {
    
    const std::string input = R"(
this part
will be ignored

SOLuTION MgCl2 this is not detected because of lower case letters
Mg 0.5
Cl 1.0
/ end

SOLUTION NaCl
Na 0.1
Cl 0.1
/ end

hi

this
# a comment line, ignore...
is
a

BRINE Brine 2
Ca 0.1
Cl 0.2
/ end

test
* aaand another comment

* comment line: this is not read either, b/c of upper case characters in end marker
BRINE BaCl2
Ba 0.4
Cl 0.8
/ eND
# Note putting this part in-between some of the correct keywords would break things....

ignore this
too!
    )";

    std::istringstream inputStream(input);
    
    InputReader inpReader({}, /* case_sensitive= */ true);
    inpReader.define_block_keyword("SOLUTION", { "/ end" });
    inpReader.define_block_keyword("BRINE", { "/ end" });
    inpReader.read(inputStream);
    
    // NB: The input block keyword itself is always used as a prefix, but we also need the next word on the same line...
    // ...moreover, we need to match the case as well!
    REQUIRE(inpReader.get_block_keyword_content("SOLUTION").size() == 0);
    REQUIRE(inpReader.get_block_keyword_content("SOLUTION NACL").size() == 0);
    REQUIRE(inpReader.get_block_keyword_content("SOLUTION MGCL2").size() == 0);
    REQUIRE(inpReader.get_block_keyword_content("SOLUTION MgCl2").size() == 0);
    REQUIRE(inpReader.get_block_keyword_content("BRINE").size() == 0);
    REQUIRE(inpReader.get_block_keyword_content("BRINE BRINE").size() == 0);
    REQUIRE(inpReader.get_block_keyword_content("BRINE BACL2").size() == 0);
    REQUIRE(inpReader.get_block_keyword_content("BRINE BaCl2").size() == 0);
    
    const auto& solutionContents = inpReader.get_block_keyword_content("SOLUTION NaCl");
    REQUIRE(solutionContents.size() == 2);
    REQUIRE(solutionContents[0].compare("Na 0.1") == 0);
    REQUIRE(solutionContents[1].compare("Cl 0.1") == 0);
    
    // Again, the "2" is currently not accounted for... (see previous test)
    const auto& brineContents = inpReader.get_block_keyword_content("BRINE Brine");
    REQUIRE(brineContents.size() == 2);
    REQUIRE(brineContents[0].compare("Ca 0.1") == 0);
    REQUIRE(brineContents[1].compare("Cl 0.2") == 0);
    
}
